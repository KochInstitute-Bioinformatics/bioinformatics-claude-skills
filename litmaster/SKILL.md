---
name: litmaster
description: Use when a literature answer about a gene, cell type, pathway or disease question will be shown to colleagues, put in a report or email, or used to judge an analysis result (for example "what does the literature say about GENE_A and GENE_B in colon tumor stroma?"), and a missed paper, a misquoted sentence, a retracted reference or run-to-run variation of a single search would matter. Not for a quick one-fact lookup.
---

# litmaster: verified literature search (PROTOTYPE, v0.1)

**Status:** prototype (stage 7, the PDF stage, was added and tested last). Tested without a model on one known development case (one gene in colon stroma), and then once end to end with a model on two genes in the same tissue (3 discovery agents, 39 verified quotes, report passing the checker). Original no-model test: Stage 1 found 13 of the 18 papers cited in a hand-verified fact-check (the design was tuned on those papers, so this is optimistic); stages 4 to 6 reproduced 27 of 27 hand-verified quotes exactly and caught 5 of 5 deliberate report errors. **Not tested yet:** stage 2 (discovery sub-agents), a full run where a model screens, picks sentences and writes the report, and any question the design was not tuned on. Say so when you report results.

## Core idea
Scripts do retrieval and quoting; the agent does judgment. **The agent never types a quote.** It picks sentence IDs; a script cuts the quote out of the source text, re-checks it against a freshly fetched copy, and checks retractions. A final script fails the report if any quoted text is not in the evidence table.

## Failures this prevents (all measured in a development comparison of 3 search setups, 2026-11)
| Observed failure | Mechanism here |
|---|---|
| One search or one agent session missed key papers; identical runs cited 3 to 10 papers | Stage 1 script (wide, synonyms, full-text-term query) **plus** several independent agent angles, pooled |
| Key paper names the gene only in its full text | Full-text-term query (`FT.*` tags); `--grep` for sentences with no gene name (e.g. "cluster A declined") |
| Family-level wording missed (a paper says "R-spondin 2-4") | `family_terms` in the config count as gene hits |
| Quotes typed wrong or from a different abstract version | Quotes come from sentence IDs; exact-substring check against a fresh fetch |
| 100 KB tool outputs nobody could read; MCP server failed to connect in 4 of 6 sessions | Results go to files; agents see a compact list; no MCP dependency (SDK / scripts only) |
| A rate-limited or failed tool looked like "no papers" | Every call and retry is logged; failures are printed in the Verification Notes |
| A retracted paper almost cited | Crossref check + PubMed publication type on every cited paper |
| No contrary search; first answer was one-sided | Contrary-wording queries are generated in stage 1 and a contrary agent angle is required |

## Stages
Set `S=<scripts dir>`, `W=<local workdir>` (**not under /net or /home when running headless Claude sessions**), `CFG=<config.json>` (copy `config.example.json` and edit genes, context terms and the question). Run every stage on a compute node: `sbatch --export=ALL,LITMASTER_SCRIPTS=$S $S/run_stage.sbatch <script> <args>` (create `logs/` first; show the user the sbatch command before submitting).

| # | Stage | Command (arguments after the script name) | Output in W |
|---|---|---|---|
| 1 | Breadth search (script) | `s1_breadth.py $CFG $W` | candidates.json, call_log.json |
| 2 | Discovery agents (optional, **ask first: uses plan allowance**) | sub-agents from `agent_brief.md`, each calling `s1b_agent_search.py $W <angle> "<query>"` | agent_candidates/*.json |
| 3 | Merge and screen | `s2_merge_screen.py $CFG $W` | screen_table.tsv, queue_pmids.txt (all tier B first; the remaining slots shared A 45%, FT 45%, C 10%; then EVERY agent-found paper up to `agent_queue_max`, default 400), queue_all_pmids.txt, merge_summary.json |
| 4 | Read full text, index sentences (4 parallel fetches; allow about 40 min per 500 papers) | `s3_sentences.py $CFG $W` | sentences.json, sentences_overview.txt (1 line per paper that has a gene sentence), no_gene_papers.txt (read, no gene sentence anywhere), sentences_view.txt |
| 5 | Verify the chosen quotes | `s4_verify.py $CFG $W` | evidence_table.json, sources_block.md, verification_notes_auto.md |
| 6 | Check the report | `python s5_check_report.py report.md $W/evidence_table.json --require-all` | PASS / FAIL |
| 7 | PDF for a reader (optional) | `python s6_make_pdf.py $W narrative.html report.pdf --title ... --footer ...` (WeasyPrint Python) | report.pdf, report.html, page previews |

### Stage 3 to 4: screening (your judgment)
Read `screen_table.tsv`. Tier B = gene + tissue + cell + condition words in title/abstract; A = gene + tissue; FT = found by a full-text-term query; AG = found only by an agent search (an agent search says nothing about the gene being in the paper, so AG papers are all read and then screened by stage 4); C = other. **The tiers are a filter, not relevance**: a paper can be relevant without the gene in its abstract. Fill the `decision` and `reason` columns (include / exclude / unclear) for papers you read, and put every paper you cite in the queue. Papers with no PMID or abstract are in `unassessed.json`: report them as unassessed, never as irrelevant.

### Stage 4 to 5: choosing evidence
Read `sentences_overview.txt` first (one line per paper: tier, read level, number of gene sentences, best sentence). Open the papers that matter with `--show <PMID> [--top N]` (gene sentences with sections; `--top 5` for big papers). **Do not rate a paper from its best sentence alone**: in testing, Kobayashi 2021 looked merely related from one sentence and was a direct tumor-vs-normal paper once opened. `sentences_view.txt` has the top 3 sentences of the first 30 papers. Reading many papers costs cluster time, not agent tokens, so keep the queue large (`fulltext_queue_max`, default 300; every agent-found paper is read on top of that). After changing the cap or adding agent papers, re-run stage 3, then stage 4 with `--resume` so only new papers are read. The overview is ordered by tier (B first); stop reading it when the lines stop being relevant. For more of one paper: `s3_sentences.py $CFG $W --show <PMID>`. If the key sentence has no gene name: `s3_sentences.py $CFG $W --grep <PMID> "<regex>"` registers it with an ID. Write `W/claims.json`:
`[{"id":"C1","pmid":"...","sentence_id":"<PMID>:<n>","tier":"T3","claim":"short note in your own words"}]`
Tiers (the grades used in ToolUniverse's literature skill): T1 mechanistic, T2 functional, T3 association, T4 mention or secondary source.
Rules for choosing (verification addendum R1 to R9): one quote per claim; say what quantity the paper measured and mark it INDIRECT when it differs from the question (expression vs fraction of cells, bulk vs single-cell, injury vs cancer, human vs mouse); report a finding no stronger than the authors; `not stated` and `unclear` are valid answers; never fill a gap from memory.

### Stage 6: the report
Mark every phrase taken from a source in CURLY quotes “…”; use straight quotes only for search queries and labels (the checker verifies curly-quoted text of 3+ words). Write `report.md`: Question; Answer (graded, with the opposite or conflicting evidence stated); **paste `sources_block.md` unchanged**; paste `verification_notes_auto.md` and replace every `<<FILL ...>>` (mode, spot-check result after reading `spot_check.txt`, unsupported claims removed, indirect evidence, what the contrary searches found); Limitations. **Audit your own claims before the checker:** the checker only proves that quoted text matches the evidence table; it cannot tell whether a paraphrase is supported. For every claim in the Answer, confirm the selected quote says it (in the first full pass, 3 claims went beyond their quotes: a sample size, a mutation list from an unselected sentence, an inference about fusions). Then run the checker. If it FAILS, fix the report (never edit the evidence table to make it pass).

### Stage 7: PDF report for a reader (optional)
When the user wants a document to hand to a person, read `pdf_report.md` first. In short: write a narrative from `report_template/narrative_skeleton.html` that refers to sources only by placeholder; the script generates the sources section from the evidence table and refuses to build if quoted text is not in the table. Write for the reader (no "our" or "we"), state nothing about the reader's own data that the user has not confirmed, audit every paraphrase against its quote, and **look at every rendered page** before handing over. For a technical report about the method, every number must be recomputed from a file, not recalled.

## Always
- Queries contain only gene names, disease and cell terms, PMIDs and DOIs. No unpublished data, sample identifiers or file contents go to any service.
- State the tier and the read level (full text vs abstract) for each source.
- If stage 1 or an agent reports failures (HTTP 429/5xx, timeouts), say so; an empty result from a failed search is not evidence of absence.
- A single run is one sample. If the answer matters, run stage 2 with at least two agent angles and compare what each found.

## Workdir rules
Keep `call_log.json` (stage 1) and `agent_search_log.json` (stage 2) in the workdir that stage 5 runs in; if you copy a workdir, copy them too, otherwise the Verification Notes under-report tool failures (stage 5 prints a WARNING when `call_log.json` is missing). Agent searches are logged by `s1b_agent_search.py`.

## Known gaps (v0.1)
Screening is by keyword flags plus your reading; no automatic relevance model. Citation expansion is off by default (5% on-topic in the development test). Tested on two questions (one gene; two genes, same tissue), so treat it as untuned for other topics. Sentence splitting glues figure captions and abstract section labels onto sentences (quotes stay exact but can be ugly). Crossref retraction lookups return no value intermittently (reported, not hidden).
