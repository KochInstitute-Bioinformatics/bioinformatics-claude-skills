# `litmaster`: verified literature search (prototype v0.1)

A Claude Code folder skill for literature questions whose answer will be shown to colleagues or used to judge an analysis result. Scripts do the retrieval and the quoting; Claude does the screening and judgment. Claude never types a quote: it picks sentence IDs, and a script cuts the quote from the source text, re-checks it against a freshly fetched copy and checks for retractions. A final script fails the report if any quoted text is missing from the evidence table.

**Status: prototype.** It was tested on two development questions only (one gene, then two genes, in the same tissue), so treat it as untuned for other topics. A model-driven pass was run once; the discovery sub-agent count has not been tuned. `SKILL.md` lists what is and is not tested.

## Install

Folder skill, like `/handoff`:

```bash
mkdir -p ~/.claude/skills
cp -r litmaster ~/.claude/skills/
```

## Prerequisites

- A Python with the [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SDK installed (set `LITMASTER_PY` to it).
- Network access to PubMed, Europe PMC, OpenAlex, Crossref (and optionally Semantic Scholar, which rate-limits). The scripts call these public services directly or through ToolUniverse; their terms of use apply.
- Slurm for the provided `scripts/run_stage.sbatch` (it is the only scheduler-specific file; pass your partition with `sbatch -p`). Other schedulers: run the stage scripts on a compute node by hand.
- Optional, for the PDF stage: WeasyPrint in some Python.

## Use

Copy `config.example.json`, set your genes, synonyms, context terms and question, then follow `SKILL.md` (stages 1 to 7). Review everything it writes before sharing it.

## Contents

`SKILL.md` (the skill), `agent_brief.md` (brief for discovery sub-agents), `pdf_report.md` (optional stage 7), `scripts/` (stages 1 to 6, PDF build), `report_template/`, `config.example.json`. The development tests and data are not included.

## Credits

- The search tools come from [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) (mims-harvard, Apache-2.0). The T1 to T4 evidence grades (mechanistic, functional, association, mention) follow the scheme in its literature skill.
- The verification rules (a supporting quote for every claim, "not stated" as a valid answer, visible counts, random spot-checks) draw on the habits taught in the "Prompts to Protocols" workshop (Computation Questions Module 4, 13 November 2026, presented by Alexander Greer).
- The scripts and the rest of the wording are original.
