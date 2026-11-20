"""Stage 5: turn the agent's chosen sentence IDs into verified quotes. The agent supplies claims.json:
   [{"id": "C1", "pmid": "42938326", "sentence_id": "42938326:7", "tier": "T3", "claim": "short note in the agent's own words"}, ...]
This script (never the agent) produces the quote text, re-fetches the source independently, checks each quote as an exact substring,
runs the retraction checks, draws 3 random quotes for a spot-check, and writes sources_block.md + verification_notes_auto.md.
Run on a compute node:  python s4_verify.py <config.json> <workdir> [seed]"""
import json, sys, os, random, collections, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common
from common import Cfg, epmc_records, fulltext_paragraphs, clean, log

cfg = Cfg(sys.argv[1]); W = sys.argv[2]; seed = int(sys.argv[3]) if len(sys.argv) > 3 else 11
claims = json.load(open(f"{W}/claims.json")); S = json.load(open(f"{W}/sentences.json"))
rows, errors = [], []
for c in claims:
    d = S.get(c["pmid"]); s = next((x for x in d["sentences"] if x["id"] == c["sentence_id"]), None) if d else None
    if not s: errors.append(f"{c['id']}: sentence {c['sentence_id']} not found for PMID {c['pmid']}"); continue
    w = s["text"].split(); q = " ".join(w[:40]) + (" …" if len(w) > 40 else "")
    rows.append({"claim_id": c["id"], "pmid": c["pmid"], "label": d["label"], "doi": d["doi"], "journal": d["journal"], "year": d["year"], "title": d["title"],
                 "level": d["level"], "section": s["section"], "quote": q, "truncated": len(w) > 40, "tier": c.get("tier", ""), "claim": c.get("claim", ""),
                 "authors": d.get("authors") or [d.get("first_author", "")], "n_authors": d.get("n_authors", 1), "retracted_pubtype": d["retracted_pubtype"]})
if errors:
    log("CLAIM ERRORS:\n  " + "\n  ".join(errors)); sys.exit(2)
# independent re-fetch and exact-substring check (no normalisation)
pm = sorted({r["pmid"] for r in rows}); recs = epmc_records(pm, what="s4 records"); corpus = {}
for p in pm:
    rec = recs.get(p) or {}
    paras = fulltext_paragraphs(rec.get("pmcid")) if rec.get("pmcid") else None
    corpus[p] = [clean(rec.get("abstractText"))] + [t for _, t in (paras or [])]
for r in rows:
    rec = recs.get(r["pmid"]) or {}
    if not r["journal"]: r["journal"] = rec.get("journalTitle") or ((rec.get("journalInfo") or {}).get("journal") or {}).get("title") or "journal not recorded"
    q = r["quote"].replace(" …", ""); r["exact_in_fresh_source"] = any(q in t for t in corpus[r["pmid"]])
bad = [r for r in rows if not r["exact_in_fresh_source"]]
log(f"R1: {len(rows) - len(bad)} of {len(rows)} quotes found as exact substrings of a freshly fetched source")
for r in bad: log("  FAIL:", r["label"], r["quote"][:80])
# R7 spot-check: three random quotes with context
random.seed(seed); pick = random.sample(rows, min(3, len(rows))); spot = []
for r in pick:
    q = r["quote"].replace(" …", ""); ctx = ""
    for t in corpus[r["pmid"]]:
        i = t.find(q)
        if i >= 0: ctx = t[max(0, i - 160): i + len(q) + 160]; break
    spot.append({"label": r["label"], "quote": q[:80], "context": ctx})
open(f"{W}/spot_check.txt", "w").write("\n\n".join(f"{s['label']}: ...{s['context']}..." for s in spot))
# R6 retraction: PubMed publication type (from s3) + Crossref via ToolUniverse
from tooluniverse import ToolUniverse
tu = ToolUniverse(); tu.load_tools(); retr = {}
for p in pm:
    r0 = next(r for r in rows if r["pmid"] == p)
    try:
        x = tu.run({"name": "Crossref_check_retraction", "arguments": {"doi": r0["doi"]}}); d = x.get("data", {}) if isinstance(x, dict) else {}
        retr[p] = {"crossref_retracted": d.get("is_retracted"), "concern": d.get("has_expression_of_concern"), "pubtype_retracted": r0["retracted_pubtype"]}
    except Exception as ex:
        retr[p] = {"crossref_retracted": "ERROR", "error": str(ex)[:80], "pubtype_retracted": r0["retracted_pubtype"]}
flag = [p for p, v in retr.items() if v["crossref_retracted"] is True or v["pubtype_retracted"]]
log("R6: retracted or flagged among cited:", flag or "none", "| Crossref returned no value for:", [p for p, v in retr.items() if v["crossref_retracted"] in (None, "ERROR")])
# sources block (numbered by first appearance), ready to paste into the report
order = list(dict.fromkeys(r["pmid"] for r in rows)); lines = []
for n, p in enumerate(order, 1):
    rs = [r for r in rows if r["pmid"] == p]; r0 = rs[0]; a = r0["authors"]
    who = a[0] if r0["n_authors"] == 1 else f"{a[0]} & {a[1]}" if r0["n_authors"] == 2 and len(a) > 1 else f"{a[0]} et al."
    lines.append(f"{n}. {who} *{r0['journal']}* {r0['year']}. doi:{r0['doi']}. PMID {p}. {r0['level'].capitalize()}. {'/'.join(sorted({x['tier'] for x in rs if x['tier']}))}.")
    for r in rs: lines.append(f"   - Section: {r['section']}. “{r['quote']}”")
open(f"{W}/sources_block.md", "w").write("\n".join(lines) + "\n")
json.dump({"quotes": rows, "retraction": retr, "spot_check": spot}, open(f"{W}/evidence_table.json", "w"), indent=1)
# verification notes: auto-filled facts, <<FILL>> for judgments only the agent can make
M = json.load(open(f"{W}/merge_summary.json")) if os.path.exists(f"{W}/merge_summary.json") else {}
log_s1 = json.load(open(f"{W}/call_log.json")) if os.path.exists(f"{W}/call_log.json") else []
fails = collections.Counter(c["tool"] for c in log_s1 if c["error"]); codes = collections.Counter(m for c in log_s1 if c["error"] for m in re.findall(r"\b(4\d\d|5\d\d)\b|timed out", c["error"]))
contrary = [f"{c['tool']}: {c.get('query', '')[:90]}" for c in log_s1 if re.search(r"\.E$", c["tag"])]
full = sum(1 for p in order if next(r for r in rows if r["pmid"] == p)["level"] == "full text")
tools = sorted({c["tool"] for c in log_s1} | {"Crossref_check_retraction"})
ag = json.load(open(f"{W}/agent_search_log.json")) if os.path.exists(f"{W}/agent_search_log.json") else []
ag_sum = json.load(open(f"{W}/agent_search_summary.json")) if os.path.exists(f"{W}/agent_search_summary.json") else {}   # fallback when searches were not machine-logged
if not os.path.exists(f"{W}/call_log.json"): print("WARNING: no call_log.json in the workdir: stage-1 failures and contrary queries cannot be reported; copy it from the stage-1 folder")
ag_line = (f"agent searches logged: {len(ag)} | failed {sum(1 for x in ag if x['failed'])} | zero hits {sum(1 for x in ag if not x['failed'] and x['hits'] == 0)}" if ag else
           ("agent searches (from the agents' own reports, NOT machine-logged): " + "; ".join(f"{k}: {v['searches']} searches, {v['failed']} failed, {v['zero_hit']} zero-hit" for k, v in ag_sum.items()) if ag_sum else "agent searches: not recorded"))
notes = f"""## Verification Notes
Mode: <<FILL: factoid | mini-review | full>>      Read level: full text for {full} of {len(order)} cited papers; abstract only for {len(order) - full}
Counts: unique merged {M.get('unique_merged', '?')} (script + {M.get('agent_entries', 0)} agent entries) | assessed with abstract {M.get('assessed', '?')} | unassessed (no PMID/abstract) {M.get('unassessed', '?')}
        tiers {M.get('tiers', '?')} | papers per source {M.get('papers_per_source', '?')} | read queue {M.get('queue', '?')} | cited {len(order)}
Tool failures / empty sources: failed attempts per tool in stage 1: {dict(fails) or 'none'}; codes seen: {dict(codes) or 'none'}; transient fetch failures logged in stages 3-5: {len(common.FAILURES)}; {ag_line}
References checked: {len(order)} cited | found in Europe PMC with title/year | retracted by Crossref or PubMed type: {flag or 'none'} | Crossref gave no value for: {[p for p, v in retr.items() if v['crossref_retracted'] in (None, 'ERROR')] or 'none'}
Verbatim check (R1): {len(rows) - len(bad)} of {len(rows)} quotes are exact substrings of a freshly fetched source (quotes are produced by script from sentence IDs, not typed)
Spot-check (R7), 3 random quotes re-read with context (seed {seed}): {', '.join(s['label'] for s in spot)} -> see spot_check.txt: <<FILL: pass/fail after reading>>
Unsupported claims removed: <<FILL>>      Indirect evidence (quantity differs from the question): <<FILL>>
Contrary-evidence searches (R8): {len(contrary)} wording-contrary queries were run in stage 1 (examples: {contrary[:2]}) -> <<FILL: what they found; conflicts>>
Sent to outside services: only generic public queries (gene names, disease/cell terms, PMIDs, DOIs) via {', '.join(tools)}; no project data.
"""
open(f"{W}/verification_notes_auto.md", "w").write(notes)
log("wrote sources_block.md, evidence_table.json, spot_check.txt, verification_notes_auto.md")
sys.exit(1 if bad else 0)
