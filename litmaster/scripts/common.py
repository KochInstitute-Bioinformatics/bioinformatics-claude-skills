"""Shared helpers for the litmaster pipeline. Python 3.10, standard library only (plus ToolUniverse in s1/s4).
All network calls retry with backoff, have timeouts and log what failed (a failed source is never silent)."""
import json, re, sys, time, urllib.parse, urllib.request, xml.etree.ElementTree as ET

UA = {"User-Agent": "litmaster/0.1 (research literature check)"}
SENT_SPLIT = re.compile(r"(?<=[.!?])\s+(?=[A-Z(\[“\"])")
FAILURES = []   # (what, error) appended by http_get / tu_call so scripts can report them


def log(*a):
    print(*a, flush=True)


def http_get(url, raw=False, tries=5, timeout=60, what=None):
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
                b = r.read()
                return b if raw else json.loads(b)
        except Exception as ex:
            FAILURES.append((what or url[:70], f"try {k+1}: {str(ex)[:60]}"))
            if k == tries - 1:
                return None
            time.sleep(4 * (k + 1))


def term_rx(t):
    """'R-spondin 3' -> regex tolerant to space/hyphen variants; case-insensitive use."""
    parts = re.split(r"[\s\-‐‑–]+", t.strip().rstrip("*"))   # trailing * marks a stem (intestin*): same regex, prefix match
    return r"[\s\-‐‑–]?".join(re.escape(p) for p in parts)


def any_rx(terms, flags=re.I, word=True):
    if not terms:
        return re.compile(r"(?!x)x")
    body = "|".join(term_rx(t) for t in terms)
    start = r"(?<![A-Za-z0-9])"                       # terms must start at a word boundary ('normal' must not match 'abnormal')
    return re.compile(start + r"(?:" + body + r")" + (r"(?![A-Za-z0-9])" if word else ""), flags)   # word=False: prefix match (stems)


class Cfg:
    """Loads config.json and exposes compiled regexes. Fields: see config.example.json."""
    def __init__(self, path):
        c = json.load(open(path)); self.raw = c
        self.name = c["name"]; self.question = c["question"]
        self.groups = c["gene_groups"]                       # {"GENE_A": [terms], "GENE_B": [terms]}
        self.family = c.get("family_terms", [])              # weaker, e.g. "R-spondin", "R-spondins"
        self.ctx = c["context_terms"]                        # tissue, cell, condition, model, method, contrary
        self.limit = int(c.get("limit", 40))
        self.sources = c.get("sources", ["pubmed", "europepmc", "openalex", "semanticscholar"])
        self.extra_queries = c.get("extra_queries", [])
        self.seed_pmids = c.get("seed_pmids", [])
        self.queue_max = int(c.get("fulltext_queue_max", 60))
        self.gene_rx = any_rx([t for g in self.groups.values() for t in g] + self.family)
        self.group_rx = {g: any_rx(t) for g, t in self.groups.items()}
        self.tissue_rx = any_rx(self.ctx["tissue"], word=False)
        self.cell_rx = any_rx(self.ctx["cell"], word=False)
        self.cond_rx = any_rx(self.ctx["condition"], word=False)
        self.direction_rx = re.compile(r"(higher|lower|increas|decreas|up-?regulat|down-?regulat|enrich|reduc|elevat|over-?express|loss|deplet|abundan|declin|expand)", re.I)


def clean(s):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", s or "")).strip()


def sentences(text):
    return [s for s in SENT_SPLIT.split(text) if s.strip()]


def epmc_records(pmids, what="epmc batch"):
    """Core records (title, abstract, authors, pubTypes, pmcid...) for PMIDs, in batches of 25."""
    out = {}
    pm = [p for p in pmids if p]
    for i in range(0, len(pm), 25):
        chunk = pm[i:i + 25]
        q = "(" + " OR ".join(f"EXT_ID:{p}" for p in chunk) + ") AND SRC:MED"
        d = http_get("https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&resultType=core&pageSize=100&query=" + urllib.parse.quote(q), what=what)
        for r in ((d or {}).get("resultList", {}).get("result", [])):
            out[str(r.get("pmid") or r.get("id"))] = r
    return out


def epmc_by_doi(doi):
    d = http_get("https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&resultType=core&query=" + urllib.parse.quote(f'DOI:"{doi}"'), what="epmc doi")
    rs = (d or {}).get("resultList", {}).get("result", [])
    return rs[0] if rs else None


def fulltext_paragraphs(pmcid):
    """[(section path, paragraph text)] from open-access full text (Europe PMC XML, then NCBI efetch); None if unavailable."""
    if not pmcid:
        return None
    for url in (f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML",
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={pmcid}&retmode=xml"):
        raw = http_get(url, raw=True, tries=3, what=f"fulltext {pmcid}")
        if not raw:
            continue
        try:
            root = ET.fromstring(raw)
        except ET.ParseError:
            continue
        parent = {c: p for p in root.iter() for c in p}
        paras = []
        for p in root.iter("p"):
            sec, n, skip = "", p, False
            while n in parent:
                n = parent[n]
                if n.tag == "sec":
                    t = n.find("title"); sec = (clean("".join(t.itertext())) if t is not None else "") + (" > " + sec if sec else "")
                if n.tag in ("ref-list",):
                    skip = True; break
                if n.tag == "abstract":
                    sec = "abstract" + (" > " + sec if sec else ""); break
            if not skip:
                paras.append((sec or "(front/other)", clean("".join(p.itertext()))))
        if paras:
            return paras
    return None


def label_of(rec):
    """'Kang 2016' from a Europe PMC record (first author surname + year)."""
    a = (rec.get("authorString") or "").split(",")[0].strip()
    sur = a.rsplit(" ", 1)[0] if a else "?"
    return f"{sur} {rec.get('pubYear') or ''}".strip()


def find_items(r):
    """first list of dicts that looks like papers inside a ToolUniverse result"""
    if isinstance(r, list) and r and isinstance(r[0], dict):
        return r
    if isinstance(r, dict):
        for k in ("data", "results", "articles", "papers", "works", "items", "result", "hits"):
            if k in r:
                x = find_items(r[k])
                if x:
                    return x
        for v in r.values():
            x = find_items(v) if isinstance(v, (dict, list)) else []
            if x:
                return x
    return []


def norm_item(it):
    g = lambda *ks: next((it[k] for k in ks if it.get(k)), None)
    doi = g("doi", "DOI"); pmid = g("pmid", "PMID", "pubmed_id", "id")
    year = g("year", "pubYear", "publication_year", "pub_year", "date", "pubdate")
    return {"title": clean(str(g("title", "display_name", "article_title") or "")), "year": str(year)[:4] if year else None,
            "pmid": str(pmid) if pmid and str(pmid).isdigit() else None,
            "doi": (str(doi).replace("https://doi.org/", "").lower() if doi else None)}


def grep_sentences(W, pmid, pattern, maxn=15):
    """Register extra sentences of one paper that match `pattern` (case-insensitive), even without a gene name in them
    (e.g. 'cluster A showed a marked decline'). IDs are PMID:x<n>; they are written into sentences.json so the agent still
    selects IDs instead of typing quotes. Returns the sentences added or already present."""
    S = json.load(open(f"{W}/sentences.json")); d = S[pmid]; rx = re.compile(pattern, re.I)
    rec = epmc_records([pmid], what="grep record").get(pmid) or {}
    paras = fulltext_paragraphs(rec.get("pmcid")) if rec.get("pmcid") else None
    body = [("abstract", clean(rec.get("abstractText")))] + [(s, t) for s, t in (paras or []) if not s.startswith("abstract")]
    have = {s["text"]: s for s in d["sentences"]}; out = []
    for sec, txt in body:
        for st in sentences(txt):
            if rx.search(st):
                if st in have: out.append(have[st]); continue
                n = sum(1 for s in d["sentences"] if ":x" in s["id"]) + 1
                new = {"id": f"{pmid}:x{n}", "section": sec, "text": st, "score": 0}; d["sentences"].append(new); have[st] = new; out.append(new)
            if len(out) >= maxn: break
    json.dump(S, open(f"{W}/sentences.json", "w"), indent=1)
    return out


def epmc_by_dois(dois, what="epmc doi batch"):
    """DOI -> Europe PMC core record, in batches of 15 (one request per batch instead of one per DOI)."""
    out = {}
    ds = [d for d in dois if d]
    for i in range(0, len(ds), 15):
        chunk = ds[i:i + 15]
        q = "(" + " OR ".join(f'DOI:"{d}"' for d in chunk) + ")"
        d = http_get("https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&resultType=core&pageSize=100&query=" + urllib.parse.quote(q), what=what)
        for r in ((d or {}).get("resultList", {}).get("result", [])):
            if r.get("doi"):
                out[r["doi"].lower()] = r
    return out
