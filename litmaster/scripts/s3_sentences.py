"""Stage 4: read the full text (open access) of the queued papers and index every sentence that mentions the config's genes, with its
section and a stable ID (PMID:n). The agent later selects sentence IDs; it never types a quote.
Run on a compute node:  python s3_sentences.py <config.json> <workdir> [pmids_file]        (default pmids_file: queue_pmids.txt)
Show more sentences of one paper:  python s3_sentences.py <config.json> <workdir> --show <PMID>
Add sentences that lack a gene name:  python s3_sentences.py <config.json> <workdir> --grep <PMID> "<regex>"   (registered with IDs PMID:x<n>)"""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common
from common import Cfg, epmc_records, fulltext_paragraphs, sentences, clean, label_of, log

cfg = Cfg(sys.argv[1]); W = sys.argv[2]
if len(sys.argv) > 4 and sys.argv[3] == "--show":
    d = json.load(open(f"{W}/sentences.json"))[sys.argv[4]]
    top = int(sys.argv[sys.argv.index("--top") + 1]) if "--top" in sys.argv else None      # --top N: only the N best-scoring sentences (big papers have hundreds)
    print(f"{d['label']} | {d['title']} | {d['level']} | {len(d['sentences'])} gene sentences")
    for s in (sorted(d["sentences"], key=lambda s: -s["score"])[:top] if top else d["sentences"]): print(f"{s['id']} | {s['section'][:45]} | {s['text']}")
    sys.exit(0)
if len(sys.argv) > 5 and sys.argv[3] == "--grep":      # --grep <PMID> "<regex>": add sentences without a gene name, e.g. "cluster A"
    for s in common.grep_sentences(W, sys.argv[4], sys.argv[5]): print(f"{s['id']} | {s['section'][:45]} | {s['text']}")
    sys.exit(0)
resume = "--resume" in sys.argv; sys.argv = [a for a in sys.argv if a != "--resume"]      # --resume: keep papers already in sentences.json, read only the new ones
pm_file = sys.argv[3] if len(sys.argv) > 3 else f"{W}/queue_pmids.txt"
pmids = [l.split("\t")[0].strip() for l in open(pm_file) if l.strip() and not l.startswith("pmid")]
prev = json.load(open(f"{W}/sentences.json")) if resume and os.path.exists(f"{W}/sentences.json") else {}
todo = [p for p in pmids if p not in prev]; log(f"queue {len(pmids)} | already indexed (resume) {len(pmids) - len(todo)} | to read now {len(todo)}")
recs = epmc_records(todo, what="s3 records")
tier = {}
if os.path.exists(f"{W}/screen_table.tsv"):
    for l in open(f"{W}/screen_table.tsv").read().splitlines()[1:]:
        c = l.split("\t"); tier[c[0]] = c[2]
from concurrent.futures import ThreadPoolExecutor


def process(pm):
    rec = recs.get(pm)
    if not rec: return pm, None
    paras = fulltext_paragraphs(rec.get("pmcid")) if rec.get("pmcid") else None
    level = "full text" if paras else "abstract only"
    ab = clean(rec.get("abstractText"))
    body = [("abstract", ab)] + [(s, t) for s, t in (paras or []) if not s.startswith("abstract")]
    sents, n = [], 0
    for sec, txt in body:
        for st in sentences(txt):
            if cfg.gene_rx.search(st):
                sc = 2 + bool(cfg.direction_rx.search(st)) + bool(cfg.cond_rx.search(st) or cfg.cell_rx.search(st) or cfg.tissue_rx.search(st)) + (sec.split(" > ")[0].lower().startswith(("abstract", "result", "discussion")))
                n += 1; sents.append({"id": f"{pm}:{n}", "section": sec, "text": st, "score": sc})
    auth = (rec.get("authorString") or "").rstrip(".").split(", ")
    pubtypes = rec.get("pubTypeList", {}).get("pubType", []) if isinstance(rec.get("pubTypeList"), dict) else []
    return pm, {"label": label_of(rec), "title": clean(rec.get("title")), "journal": rec.get("journalTitle"), "year": rec.get("pubYear"), "doi": rec.get("doi"),
                "pmcid": rec.get("pmcid"), "level": level, "first_author": auth[0] if auth else "", "authors": auth[:3], "n_authors": len(auth),
                "tier": tier.get(pm, ""), "retracted_pubtype": any("Retracted" in p for p in pubtypes), "sentences": sents}


out = {}
with ThreadPoolExecutor(max_workers=4) as ex:          # 4 parallel fetches: reading many papers is cluster time, not agent tokens
    new = {}
    for i, (pm, d) in enumerate(ex.map(process, todo), 1):
        if d is None: log(f"{pm}: no record"); continue
        new[pm] = d; log(f"[{i}/{len(todo)}] {pm} {d['label']:<24} {d['level']:<13} gene sentences: {len(d['sentences'])}")
for pm in pmids:                       # keep queue order; refresh the tier of papers carried over
    if pm in new: out[pm] = new[pm]
    elif pm in prev: out[pm] = prev[pm]; out[pm]["tier"] = tier.get(pm, out[pm].get("tier", ""))
json.dump(out, open(f"{W}/sentences.json", "w"), indent=1)
best = lambda d: max(d["sentences"], key=lambda s: s["score"]) if d["sentences"] else None
kept = {pm: d for pm, d in out.items() if d["sentences"]}; dropped = [pm for pm, d in out.items() if not d["sentences"]]
open(f"{W}/no_gene_papers.txt", "w").write("\n".join(f"{pm}\t{out[pm]['label']}\t{out[pm]['title'][:90]}" for pm in dropped) + "\n")   # read, but no gene sentence anywhere: kept for audit, not shown
log(f"screen: {len(kept)} of {len(out)} read papers have at least one gene sentence; {len(dropped)} have none (listed in no_gene_papers.txt)")
with open(f"{W}/sentences_overview.txt", "w") as f:       # ONE line per paper WITH a gene sentence: what the agent reads first
    f.write("pmid | label | tier | level | n gene sentences | best sentence id | best sentence\n")
    for pm, d in kept.items():
        b = best(d); f.write(f"{pm} | {d['label']} | {d['tier']} | {d['level']} | {len(d['sentences'])} | {b['id']} | {b['text'][:220]}\n")
with open(f"{W}/sentences_view.txt", "w") as f:           # top 3 sentences for the first 30 papers (queue order: best tiers first); use --show PMID for the rest
    for pm, d in list(kept.items())[:30]:
        f.write(f"\n## {pm} | {d['label']} | {d['level']} | {d['title'][:100]} | {len(d['sentences'])} gene sentences\n")
        for s in sorted(d["sentences"], key=lambda s: -s["score"])[:3]:
            f.write(f"{s['id']} | {s['section'][:40]} | {s['text'][:300]}\n")
log("failures:", len(common.FAILURES), "| wrote sentences.json, sentences_overview.txt, sentences_view.txt, no_gene_papers.txt")
json.dump(common.FAILURES[:100], open(f"{W}/failures_s3.json", "w"))
