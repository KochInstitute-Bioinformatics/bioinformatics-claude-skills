"""Stage 6: mechanical check of the finished report against evidence_table.json. Light (no network): fine on any node.
Usage: python s5_check_report.py <report.md> <evidence_table.json> [--require-all]
Fails (exit 1) if: a quoted span of 4+ words is not part of an evidence quote; a PMID in the report is not in the table;
the Verification Notes block is missing or still contains <<FILL;
Convention: text taken from a source goes in CURLY quotes “…” (checked from 3 words); straight quotes are for queries and labels (checked only from 10 words).
 (with --require-all) an evidence quote is missing from the report."""
import json, re, sys

def tn(t):   # typographic normalisation only: dash/minus variants, curly vs straight quote marks, non-breaking spaces; never letters, digits or case
    t = re.sub("[\u2010\u2011\u2012\u2013\u2014\u2212]", "-", t.replace("\u00a0", " "))
    return re.sub(r"\s+", " ", t.replace("\u201c", '"').replace("\u201d", '"').replace("\u2018", "'").replace("\u2019", "'"))
rep_raw = open(sys.argv[1], encoding="utf8").read(); rep = tn(rep_raw); E = json.load(open(sys.argv[2]))
quotes = [tn(q["quote"].replace(" …", "")).strip() for q in E["quotes"]]; pmids = {q["pmid"] for q in E["quotes"]}
def spans(text):
    out, depth, start = [], 0, 0
    for i, ch in enumerate(text):               # outermost curly-quoted spans (nested quotes are balanced)
        if ch == "“":
            if depth == 0: start = i + 1
            depth += 1
        elif ch == "”" and depth > 0:
            depth -= 1
            if depth == 0: out.append(text[start:i])
    return [("curly", x) for x in out] + [("straight", x) for x in re.findall(r'(?<![A-Za-z0-9])"([^"\n]{40,}?)"(?![A-Za-z0-9])', text)]
problems = []
for kind, s in spans(rep_raw):          # spans are found in the raw text (curly quotes mark source text), compared after normalisation
    core = tn(re.sub(r"\s*…$", "", s.strip()))
    if len(core.split()) >= (3 if kind == "curly" else 10) and not any(core.lower() in q.lower() for q in quotes):
        problems.append(f"quoted text not in the evidence table: “{core[:90]}”")
for p in sorted(set(re.findall(r"PMID[ :]*(\d{6,9})", rep))):
    if p in pmids: continue
    ctx = [l for l in rep.splitlines() if p in l]      # PMIDs on an 'excluded / retracted / do not cite' line are allowed
    if not any(re.search(r"retract|exclud|do not cite", l, re.I) for l in ctx):
        problems.append(f"PMID {p} is in the report but not in the evidence table")
if "Verification Notes" not in rep: problems.append("Verification Notes block is missing")
if "<<FILL" in rep: problems.append(f"{rep.count('<<FILL')} unfilled <<FILL>> placeholder(s) remain")
if "--require-all" in sys.argv:
    for q in quotes:
        if q not in rep: problems.append(f"evidence quote missing or altered in the report: “{q[:80]}”")
print(f"checked {len(spans(rep_raw))} quoted spans against {len(quotes)} evidence quotes")
print("PASS" if not problems else "FAIL:\n  " + "\n  ".join(problems))
sys.exit(1 if problems else 0)
