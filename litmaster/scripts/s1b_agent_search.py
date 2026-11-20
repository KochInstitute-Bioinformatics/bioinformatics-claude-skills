"""Stage 2 helper for discovery agents: run ONE literature search, append the hits to agent_candidates/<angle>.json, print a compact list.
Keeps tool output small (no 100 KB dumps) and writes the file for the agent (no shell redirection needed).
Usage:  python s1b_agent_search.py <workdir> <angle> "<query>" [--tool europepmc|pubmed|openalex|semanticscholar] [--limit 15] [--ft "term1,term2"]
  --ft: Europe PMC only; require these terms in the FULL TEXT (finds papers whose abstract never names the gene)."""
import json, sys, os, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import find_items, norm_item
def log_search(workdir, angle, tool, query, hits, failed, err=""):      # every agent search is logged, also failures and empty results
    f = f"{workdir}/agent_search_log.json"; L = json.load(open(f)) if os.path.exists(f) else []
    L.append({"angle": angle, "tool": tool, "query": query[:200], "hits": hits, "failed": failed, "error": err[:120]}); json.dump(L, open(f, "w"), indent=1)
ap = argparse.ArgumentParser(); ap.add_argument("workdir"); ap.add_argument("angle"); ap.add_argument("query")
ap.add_argument("--tool", default="europepmc"); ap.add_argument("--limit", type=int, default=15); ap.add_argument("--ft", default="")
a = ap.parse_args()
T = {"europepmc": ("EuropePMC_search_articles", {"query": a.query, "limit": a.limit}), "pubmed": ("PubMed_search_articles", {"query": a.query, "limit": a.limit, "include_abstract": False}),
     "openalex": ("openalex_search_works", {"search": a.query, "per_page": min(a.limit, 25)}), "semanticscholar": ("SemanticScholar_search_papers", {"query": a.query, "limit": min(a.limit, 25)})}
name, args = T[a.tool]
if a.ft and a.tool == "europepmc": args.update({"require_has_ft": True, "fulltext_terms": [t.strip() for t in a.ft.split(",") if t.strip()]})
from tooluniverse import ToolUniverse
tu = ToolUniverse(); tu.load_tools()
try:
    r = tu.run({"name": name, "arguments": args})
except Exception as e:
    log_search(a.workdir, a.angle, a.tool, a.query, 0, True, f"{type(e).__name__}: {e}"); print(f"SEARCH FAILED ({type(e).__name__}): {str(e)[:150]}"); sys.exit(1)
if isinstance(r, dict) and (r.get("error") or r.get("status") == "error"):
    log_search(a.workdir, a.angle, a.tool, a.query, 0, True, str(r.get("error") or r)); print("SEARCH FAILED:", str(r.get("error") or r)[:200], "(rate limit or outage: do not treat as 'no literature')"); sys.exit(1)
items = [norm_item(x) for x in find_items(r)]
log_search(a.workdir, a.angle, a.tool, a.query, len(items), False)
os.makedirs(f"{a.workdir}/agent_candidates", exist_ok=True); f = f"{a.workdir}/agent_candidates/{a.angle}.json"
have = json.load(open(f)) if os.path.exists(f) else []
keyf = lambda x: x.get("pmid") or x.get("doi") or x["title"].lower()[:80]
s1 = set()                                   # papers the stage-1 script already found: not "new" for the agent
if os.path.exists(f"{a.workdir}/candidates.json"):
    for x in json.load(open(f"{a.workdir}/candidates.json")): s1 |= {k for k in (x.get("pmid"), x.get("doi"), (x.get("title") or "").lower()[:80]) if k}
seen = {keyf(x) for x in have}; new = 0
print(f"{a.tool} | {a.query[:90]} | {len(items)} hits   (NEW = not found by stage 1 and not by you before; s1 = stage 1 already has it)")
for it in items:
    k = keyf(it); flag = "seen" if k in seen else ("s1  " if (k in s1 or (it["doi"] and it["doi"] in s1) or it["title"].lower()[:80] in s1) else "NEW ")
    if k not in seen: have.append(it); seen.add(k)      # every hit is recorded (also stage-1 ones, so the merge counts the extra source)
    new += flag == "NEW "
    print(f"{flag} {it['pmid'] or '-':<9}{it['year'] or '----'}  {it['title'][:95]}")
json.dump(have, open(f, "w"), indent=1); print(f"-> {new} NEW (beyond stage 1); {len(have)} total recorded in agent_candidates/{a.angle}.json")
