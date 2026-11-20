"""Stage 1: breadth search through the ToolUniverse SDK (option 2's role), generated from the config's term lists.
Run on a compute node:  <venv python> -I s1_breadth.py <config.json> <outdir>
Writes: candidates.json (with provenance), call_log.json, candidates_compact.tsv. Nothing is dropped here."""
import json, sys, time, collections, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import Cfg, find_items, norm_item, log

cfg = Cfg(sys.argv[1]); out = sys.argv[2]; os.makedirs(out, exist_ok=True)
from tooluniverse import ToolUniverse
tu = ToolUniverse(); tu.load_tools()
calls, cands = [], {}


def call(name, args, tag):
    """One tool call with up to 3 attempts on rate limits / transient errors. Every attempt is logged."""
    for attempt in range(3):
        t0 = time.time(); err = None
        try:
            r = tu.run({"name": name, "arguments": args})
            if isinstance(r, dict) and (r.get("error") or r.get("status") == "error"):
                err = str(r.get("error") or r)[:160]
        except Exception as e:
            r, err = None, f"{type(e).__name__}: {str(e)[:140]}"
        items = [] if err else find_items(r)
        calls.append({"tool": name, "tag": tag, "query": str(args.get("query") or args.get("search") or "")[:200], "attempt": attempt + 1, "sec": round(time.time() - t0, 1), "n_items": len(items), "error": err})
        log(f"{name:<30} {tag:<10} try{attempt+1} {calls[-1]['sec']:5.1f}s items={len(items):<3} {('ERR ' + err[:80]) if err else ''}")
        if not err or not any(k in err.lower() for k in ("429", "timed out", "timeout", "500", "502", "503", "504", "unavailable")):
            return items
        time.sleep(8 * (attempt + 1))
    return []


def add(items, tag):
    for it in items:
        n = norm_item(it)
        if not n["title"]:
            continue
        k = n["pmid"] or n["doi"] or n["title"].lower()[:80]
        cands.setdefault(k, {**n, "found_by": set()})["found_by"].add(tag)


def or_ep(terms):
    return "(" + " OR ".join(f'"{t}"' if re.search(r"[\s-]", t) else t for t in terms) + ")"


def or_pm(terms):
    return "(" + " OR ".join(f'"{t}"[Title/Abstract]' if re.search(r"[\s-]", t) else f"{t}[Title/Abstract]" for t in terms) + ")"


def oa_or(terms):
    return "(" + " OR ".join(f'"{t.rstrip("*")}"' if re.search(r"[\s-]", t) else t.rstrip("*") for t in terms) + ")"


C = cfg.ctx
groups = dict(cfg.groups)
if len(groups) >= 2:
    groups["BOTH"] = None            # co-occurrence of all genes
fam = {"FAMILY": cfg.family} if cfg.family else {}
templates = {   # id -> context lists (AND-ed); contrary words only in E
    "A": ["tissue", "cell"], "B": ["tissue", "method"], "C": ["tissue", "condition"], "D": ["model", "tissue", "cell"], "E": ["tissue", "cell", "contrary"],
}
for gname, terms in {**groups, **fam}.items():
    if gname == "BOTH":
        gl = [t for g in cfg.groups.values() for t in g[:2]]
        ep_g = " AND ".join(or_ep(g) for g in cfg.groups.values()); pm_g = " AND ".join(or_pm(g) for g in cfg.groups.values())
    else:
        ep_g, pm_g = or_ep(terms), or_pm(terms)
    for tid, ctxs in templates.items():
        if gname in ("BOTH", "FAMILY") and tid not in ("A", "C"):
            continue
        tag = f"{gname}.{tid}"
        if "pubmed" in cfg.sources:
            q = f"{pm_g} AND " + " AND ".join(or_pm(C[c]) for c in ctxs)
            add(call("PubMed_search_articles", {"query": q, "limit": cfg.limit, "include_abstract": False}, "PM." + tag), "PM." + tag)
        if "europepmc" in cfg.sources:
            q = f"TITLE_ABS:{ep_g} AND " + " AND ".join(f"TITLE_ABS:{or_ep(C[c])}" for c in ctxs)
            add(call("EuropePMC_search_articles", {"query": q, "limit": cfg.limit}, "EP." + tag), "EP." + tag)
    if gname not in ("BOTH",):
        # full-text-term query: the gene may appear only in the full text (this found papers the title/abstract searches missed)
        ftq = f"{or_ep(C['method'])} AND {or_ep(C['tissue'])} AND {or_ep(C['cell'])}"
        add(call("EuropePMC_search_articles", {"query": ftq, "limit": cfg.limit, "require_has_ft": True, "fulltext_terms": (terms or [])[:3]}, f"FT.{gname}"), f"FT.{gname}")
        if "openalex" in cfg.sources:
            s = f"{oa_or(terms[:3])} AND {oa_or(C['tissue'][:3])} AND {oa_or(C['cell'][:4])}"      # OpenAlex: parenthesised boolean, no wildcards
            add(call("openalex_search_works", {"search": s, "per_page": min(cfg.limit, 25)}, f"OA.{gname}"), f"OA.{gname}")
        if "semanticscholar" in cfg.sources:
            add(call("SemanticScholar_search_papers", {"query": f"{terms[0]} {C['tissue'][1].rstrip('*')} {C['cell'][0].rstrip('*')}", "limit": 25}, f"SS.{gname}"), f"SS.{gname}")
for i, q in enumerate(cfg.extra_queries):
    add(call("EuropePMC_search_articles", {"query": q, "limit": cfg.limit}, f"EP.X{i}"), f"EP.X{i}")
for pmid in cfg.seed_pmids:          # optional citation expansion (low precision in the development tests: off unless seeds are given)
    add(call("PubMed_get_cited_by", {"pmid": pmid, "limit": 20}, f"CB.{pmid}"), f"CB.{pmid}")

tab = sorted(cands.values(), key=lambda r: -len(r["found_by"]))
for r in tab:
    r["found_by"] = sorted(r["found_by"])
json.dump(tab, open(f"{out}/candidates.json", "w"), indent=1)
json.dump(calls, open(f"{out}/call_log.json", "w"), indent=1)
with open(f"{out}/candidates_compact.tsv", "w") as f:
    f.write("pmid\tyear\tn_tags\tfound_by\ttitle\n")
    for r in tab:
        f.write(f"{r['pmid'] or r['doi'] or ''}\t{r['year'] or ''}\t{len(r['found_by'])}\t{','.join(r['found_by'])[:60]}\t{r['title'][:110]}\n")
src = collections.Counter(s.split(".")[0] for r in tab for s in r["found_by"])
errs = [c for c in calls if c["error"]]
log(f"\nunique candidates: {len(tab)} | calls: {len(calls)} | failed attempts: {len(errs)} | tool time: {round(sum(c['sec'] for c in calls))} s | papers per source: {dict(src)}")
