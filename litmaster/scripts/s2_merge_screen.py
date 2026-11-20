"""Stage 3: merge candidates from s1 (script) and optional agent discovery (agent_candidates/*.json), fetch title/abstract records,
flag every paper with the config's rules, and build a screening table + a full-text reading queue. Nothing is deleted:
papers without a PMID/abstract are listed as 'unassessed' so they are not mistaken for irrelevant.
Run on a compute node:  python s2_merge_screen.py <config.json> <workdir>"""
import json, sys, os, glob, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common
from common import Cfg, epmc_records, epmc_by_doi, clean, log

cfg = Cfg(sys.argv[1]); W = sys.argv[2]
cands = {}
def key(r): return r.get("pmid") or r.get("doi") or (r.get("title") or "").lower()[:80]
def add(r, tags):
    k = key(r)
    if not k: return
    e = cands.setdefault(k, {"pmid": r.get("pmid"), "doi": r.get("doi"), "title": r.get("title"), "year": r.get("year"), "found_by": set()})
    e["found_by"] |= set(tags)
    for f in ("pmid", "doi", "title", "year"):
        if not e.get(f) and r.get(f): e[f] = r[f]
for r in json.load(open(f"{W}/candidates.json")): add(r, r["found_by"])
n_agent = 0
for f in sorted(glob.glob(f"{W}/agent_candidates/*.json")):
    for r in json.load(open(f)):
        r = {k: (str(v) if v else None) for k, v in r.items() if k in ("pmid", "doi", "title", "year")}
        add(r, ["agent:" + os.path.basename(f)[:-5]]); n_agent += 1
log(f"merged: {len(cands)} unique | agent entries: {n_agent}")
# DOI-only entries: resolve a PMID in batches, then merge entries that now share a PMID
doi_only = {e["doi"]: e for e in cands.values() if not e["pmid"] and e["doi"]}
hits = common.epmc_by_dois(list(doi_only)); resolved = 0
for d, e in doi_only.items():
    r = hits.get(d.lower())
    if r and r.get("pmid"): e["pmid"] = str(r["pmid"]); resolved += 1
merged = {}
for e in cands.values():
    k = e["pmid"] or e["doi"] or (e["title"] or "").lower()[:80]
    if k in merged:
        merged[k]["found_by"] |= e["found_by"]
        for f in ("pmid", "doi", "title", "year"):
            if not merged[k].get(f) and e.get(f): merged[k][f] = e[f]
    else: merged[k] = e
log(f"DOI-only entries: {len(doi_only)} | resolved to a PMID: {resolved} | after merging same-PMID entries: {len(cands)} -> {len(merged)}"); cands = merged
recs = epmc_records([e["pmid"] for e in cands.values()], what="s2 records")
rows, unassessed = [], []
for k, e in cands.items():
    rec = recs.get(e["pmid"]) if e["pmid"] else None
    tags = sorted(e["found_by"]); srcs = sorted({t.split(".")[0] for t in tags})
    if not rec or not (rec.get("abstractText") or rec.get("title")):
        unassessed.append({"key": k, "title": e["title"], "doi": e["doi"], "found_by": tags}); continue
    ta = clean((rec.get("title") or "") + " " + (rec.get("abstractText") or ""))
    g = [n for n, rx in cfg.group_rx.items() if rx.search(ta)]
    gene = bool(cfg.gene_rx.search(ta)); tis, cel, con = (bool(rx.search(ta)) for rx in (cfg.tissue_rx, cfg.cell_rx, cfg.cond_rx))
    ft_found = any(t.startswith("FT.") for t in tags); agent = any(t.startswith("agent:") for t in tags)       # an agent search says nothing about the gene being in the full text: own tier
    tier = "B" if (gene and tis and cel and con) else "A" if (gene and tis) else "FT" if ft_found else "AG" if agent else "C"
    pubtypes = rec.get("pubTypeList", {}).get("pubType", []) if isinstance(rec.get("pubTypeList"), dict) else []
    rows.append({"pmid": e["pmid"], "year": rec.get("pubYear"), "tier": tier, "agent": agent, "n_src": len(srcs), "genes": "+".join(g) or ("family" if gene else "-"),
                 "found_by": ",".join(tags)[:70], "oa": rec.get("isOpenAccess") or "N", "retracted": any("Retracted" in p for p in pubtypes),
                 "meeting_abstract": any("Congress" in p or "Meeting" in p for p in pubtypes), "title": clean(rec.get("title"))[:110]})
order = {"B": 0, "A": 1, "FT": 2, "AG": 3, "C": 4}
rows.sort(key=lambda r: (order[r["tier"]], -r["n_src"], -(int(r["year"]) if r["year"] and r["year"].isdigit() else 0)))
with open(f"{W}/screen_table.tsv", "w") as f:
    f.write("pmid\tyear\ttier\tn_src\tgenes\toa\tretracted\tfound_by\ttitle\tdecision\treason\n")
    for r in rows: f.write(f"{r['pmid']}\t{r['year']}\t{r['tier']}\t{r['n_src']}\t{r['genes']}\t{r['oa']}\t{r['retracted']}\t{r['found_by']}\t{r['title']}\t\t\n")
# Reading queue: EVERY tier-B paper goes in first (found missing 26 of 80 in testing when B shared quotas with the other tiers);
# the remaining slots are shared by the other tiers with quotas; unused slots roll over. Tier B is capped only if it alone exceeds fulltext_queue_max.
quota = {t: v for t, v in cfg.raw.get("queue_quota", {"A": 0.45, "FT": 0.45, "C": 0.10}).items() if t != "B"}      # B is always queued in full
pool = {t: [r for r in rows if r["tier"] == t and (t != "C" or r["n_src"] >= 3)] for t in ("B", "A", "FT", "AG", "C")}
queue = pool["B"][:cfg.queue_max]; rest = cfg.queue_max - len(queue); caps = {t: int(round(quota[t] * rest)) for t in quota}
for t in quota: queue += pool[t][:caps[t]]
left = cfg.queue_max - len(queue)
for t in quota:
    if left <= 0: break
    extra = pool[t][caps[t]:][:left]; queue += extra; left -= len(extra)
inq = {r["pmid"] for r in queue}; agent_extra = [r for r in rows if r["agent"] and r["pmid"] not in inq][:int(cfg.raw.get("agent_queue_max", 400))]
queue += agent_extra          # every agent-found paper is read (cluster time only); stage 4 then keeps those that mention a gene
open(f"{W}/queue_pmids.txt", "w").write("\n".join(r["pmid"] for r in queue) + "\n")
open(f"{W}/queue_all_pmids.txt", "w").write("\n".join(r["pmid"] for t in pool for r in pool[t]) + "\n")     # everything eligible; raise fulltext_queue_max to read more
queue_stats = {t: {"eligible": len(pool[t]), "queued": sum(1 for r in queue if r["tier"] == t)} for t in pool}
json.dump(unassessed, open(f"{W}/unassessed.json", "w"), indent=1)
by_tier = collections.Counter(r["tier"] for r in rows)
by_src = collections.Counter(s for e in cands.values() for s in {t.split(".")[0] for t in e["found_by"]})
summ = {"unique_merged": len(cands), "assessed": len(rows), "unassessed": len(unassessed), "tiers": dict(by_tier), "papers_per_source": dict(by_src),
        "queue": len(queue), "agent_papers_added_to_queue": len(agent_extra), "queue_by_tier": queue_stats, "retracted_flagged": [r["pmid"] for r in rows if r["retracted"]], "meeting_abstracts": [r["pmid"] for r in rows if r["meeting_abstract"]],
        "agent_entries": n_agent, "failures": common.FAILURES[:50]}
json.dump(summ, open(f"{W}/merge_summary.json", "w"), indent=1)
log(json.dumps({k: v for k, v in summ.items() if k != "failures"}, indent=1)); log("transient fetch failures logged:", len(common.FAILURES))
