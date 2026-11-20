"""Stage 7 (optional): build a PDF report for a reader from a narrative written by the model and the verified evidence table.
The sources section is generated here from evidence_table.json, so quotes are never retyped; the narrative only refers to them by source number.
Needs WeasyPrint (not in the ToolUniverse venv; use any Python that has it installed). Light, no network.
Usage: python s6_make_pdf.py <workdir> <narrative.html> <out.pdf> --title "..." [--subtitle "..."] [--footer "..." --date YYYY-MM-DD] [--audience outside|team]
Placeholders in the narrative: {{S:PMID,PMID}} source numbers, {{N:PMID}} one number, {{SOURCES}} the sources section (see report_template/)."""
import argparse, html, json, os, re, sys

ap = argparse.ArgumentParser(); ap.add_argument("workdir"); ap.add_argument("narrative"); ap.add_argument("out")
ap.add_argument("--title", required=True); ap.add_argument("--subtitle", default=""); ap.add_argument("--footer", default="literature report"); ap.add_argument("--date", default="")
ap.add_argument("--audience", default="outside", choices=["outside", "team"]); a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); CSS = open(f"{HERE}/../report_template/report.css", encoding="utf8").read()
E = json.load(open(f"{a.workdir}/evidence_table.json")); Q = E["quotes"]
order = list(dict.fromkeys(q["pmid"] for q in Q)); num = {p: i + 1 for i, p in enumerate(order)}
esc = lambda s: html.escape(s, quote=False)
clean_title = lambda t: re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", html.unescape(html.unescape(t or ""))))      # database titles carry entity-escaped markup
tn = lambda t: re.sub(r"\s+", " ", re.sub("[‐‑‒–—−]", "-", t.replace(" ", " ")).replace("“", '"').replace("”", '"').replace("‘", "'").replace("’", "'"))

def who(r):
    au = r.get("authors") or [""]; n = r.get("n_authors", len(au))
    return au[0] if n == 1 else f"{au[0]} &amp; {esc(au[1])}" if n == 2 and len(au) > 1 else f"{esc(au[0])} et al."
def sources_html():
    out = []
    for p in order:
        rs = [q for q in Q if q["pmid"] == p]; r0 = rs[0]; grades = "/".join(sorted({q["tier"] for q in rs if q["tier"]}))
        lvl = "abstract only" if r0["level"].startswith("abstract") else "full text"
        head = (f'<div class="src"><div class="srchead"><span class="n">{num[p]}</span> <b>{who(r0)}</b> ({r0["year"]}). {esc(clean_title(r0["title"]))} '
                f'<i>{esc(r0.get("journal") or "")}</i>. doi:{esc(r0["doi"])} &middot; PMID {p} &middot; read: {lvl} &middot; grade {grades}</div>')
        qs = "".join('<div class="quote"><div class="sec">' + esc((re.sub(r"\s+", " ", q["section"])[:95]) + ("…" if len(q["section"]) > 95 else "")) +
                     '</div>&ldquo;' + esc(q["quote"]) + '&rdquo;</div>' for q in rs)
        out.append(head + qs + "</div>")
    return "\n".join(out)

nar = open(a.narrative, encoding="utf8").read(); problems, warnings = [], []
used = set()
def res_s(m):
    ps = [x.strip() for x in m.group(1).split(",") if x.strip()]
    for p in ps:
        if p not in num: problems.append(f"{{{{S:{p}}}}}: PMID not in the evidence table"); continue
        used.add(p)
    return ", ".join(str(num[p]) for p in ps if p in num)
def res_n(m):
    p = m.group(1).strip()
    if p not in num: problems.append(f"{{{{N:{p}}}}}: PMID not in the evidence table"); return "?"
    used.add(p); return str(num[p])
body = re.sub(r"\{\{S:([^}]*)\}\}", res_s, nar); body = re.sub(r"\{\{N:([^}]*)\}\}", res_n, body)
if "{{SOURCES}}" not in body: problems.append("{{SOURCES}} placeholder missing")
if re.search(r"\[[A-Za-z][^\]]{2,}\]|PMID\d|TODO|<<FILL", re.sub(r"<!--.*?-->", "", body, flags=re.S)): warnings.append("skeleton placeholder text like [Topic] or TODO may remain in the narrative")
# 1) text taken from a source is written in CURLY quotes and must be part of an evidence quote (straight quotes are for labels and queries,
#    checked only from 10 words); this is the same convention as s5_check_report.py
quotes = [tn(q["quote"].replace(" …", "")).lower() for q in Q]
visible = re.sub(r"<!--.*?-->", "", re.sub(r"<[^>]+>", " ", html.unescape(body.replace("{{SOURCES}}", ""))), flags=re.S)
spans, depth, start = [], 0, 0
for i, ch in enumerate(visible):                       # outermost curly-quoted spans (nested quotes are balanced)
    if ch == "\u201c":
        if depth == 0: start = i + 1
        depth += 1
    elif ch == "\u201d" and depth > 0:
        depth -= 1
        if depth == 0: spans.append(("curly", visible[start:i]))
spans += [("straight", x) for x in re.findall(r'(?<![A-Za-z0-9])"([^"\n]{40,}?)"(?![A-Za-z0-9])', visible)]
for kind, sp in spans:
    core = tn(sp).lower().strip(" …")
    if len(core.split()) >= (3 if kind == "curly" else 10) and not any(core in q for q in quotes): problems.append(f"quoted text not in the evidence table: “{sp[:80]}”")
# 2) every source should be referred to; outside readers should not meet 'our' / 'we'
unref = [num[p] for p in order if p not in used]
if unref: warnings.append(f"sources never referred to in the narrative: {unref}")
if a.audience == "outside":
    hit = sorted(set(re.findall(r"\b(our|ours|we|us)\b", re.sub(r"“[^”]*”", " ", visible), flags=re.I)))
    if hit: warnings.append(f"first-person words in a report for an outside reader: {hit} (write 'the current figure', 'the project data')")
body = body.replace("{{SOURCES}}", sources_html())
for w in warnings: print("WARNING:", w)
if problems:
    print("PROBLEMS (not built):"); [print("  -", p) for p in problems]; sys.exit(1)
footer = a.footer.replace('"', "'") + (" \\2022  " + a.date if a.date else "") + " \\2022  page "      # CSS escape for the bullet; callers pass plain text
css = CSS.replace("__FOOTER__", footer)
doc = (f'<html><head><meta charset="utf-8"><title>{esc(a.title)}</title><style>{css}</style></head><body><h1>{esc(a.title)}</h1>'
       + (f'<div class="sub">{esc(a.subtitle)}</div>' if a.subtitle else "") + body + "</body></html>")
open(os.path.splitext(a.out)[0] + ".html", "w", encoding="utf8").write(doc)
from weasyprint import HTML
HTML(string=doc, base_url=".").write_pdf(a.out); print("PDF written:", a.out)
# 3) preview and text checks when the libraries are available (they are in the ToolUniverse venv, not in the WeasyPrint env)
try:
    import pypdfium2, pdfplumber
    pdf = pypdfium2.PdfDocument(a.out); pv = os.path.splitext(a.out)[0] + "_preview"; os.makedirs(pv, exist_ok=True)
    for i in range(len(pdf)): pdf[i].render(scale=1.15).to_pil().save(f"{pv}/p{i+1}.png")
    t = " ".join(pg.extract_text() or "" for pg in pdfplumber.open(a.out).pages)
    print(f"{len(pdf)} pages | leftover markup: {re.findall(r'&[a-z]+;|</?[a-z]+>', t)[:3]} | black boxes: {t.count(chr(0x25a0))} | previews in {pv}/ (look at every page)")
except ImportError:
    print("Preview skipped (pypdfium2/pdfplumber not in this Python). Render the pages with the ToolUniverse venv python and LOOK at them before sending.")
