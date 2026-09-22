# `/handoff`: Structured Session Reports

A Claude Code skill that writes a durable, schema-validated report at the end of a working session, so the next session (or the next person) can pick up where this one stopped.

Instead of a free-form Markdown note, Claude fills in a structured form of typed entries (what happened, decisions, risks, blockers, open questions, "resume here" notes, user corrections and preferences). A bundled Python script validates the form, writes it as JSON under `<project>/.claude/reports/`, and adds a one-line pointer to the project's `CLAUDE.md` so every new session sees it.

---

## Installation

Unlike the single-file skills in this repo, `/handoff` is a folder skill: it ships a validator script and reference files alongside `SKILL.md`, and it must live at `~/.claude/skills/handoff/`.

```bash
mkdir -p ~/.claude/skills
cp -r handoff ~/.claude/skills/
```

Then invoke it in Claude Code:

```
/handoff
```

Claude also picks it up from plain requests such as "write a handoff" or "write up this session".

---

## Prerequisites

| Requirement | Notes |
|-------------|-------|
| [uv](https://docs.astral.sh/uv/) | Runs the validator; it installs `pydantic` on the fly from the script's inline dependency block |
| Python ≥ 3.10 | Fetched by uv if not present |
| Linux or macOS | The `CLAUDE.md` update uses a POSIX file lock (`fcntl`) |

---

## What the skill does

1. **Gathers the session**: files touched, commands run, commits, decisions and why, risks, blockers, open questions, and whether each claim was actually tested.
2. **Splits it into typed entries**: one `session_narrative` plus one entry per decision, risk, blocker, open question or handoff note. Each entry has a summary, details, and optional guidance for the next agent, caveats, and evidence pointers (`file:line`, commit SHA, command, URL).
3. **Writes a form** at `<project>/.claude/.handoff-form.json`.
4. **Runs the validator**:
   ```bash
   uv run ~/.claude/skills/handoff/scripts/handoff_report.py --form <project>/.claude/.handoff-form.json
   ```
   which checks the form, writes `<project>/.claude/reports/YYYY-MM-DD-<slug>.json`, and inserts a pointer at the top of the `## Session reports / handoffs` section of the project `CLAUDE.md` (creating it if missing).
5. **Reports back** the report path and the index line.

On the next session start, Claude reads `CLAUDE.md`, follows the newest pointer, and reads that report in full.

---

## Safety

- Never overwrites an existing report (exit 3); pick a more specific slug instead.
- Refuses to use your home directory as the project root, or to write into your global `~/.claude/CLAUDE.md`, unless you pass `--root` or `--claude-md` explicitly.
- Concurrent `/handoff` runs cannot corrupt `CLAUDE.md` (file lock plus atomic rename).

Validator flags: `--root`, `--claude-md`, `--no-index`, `--force`, `--keep-form`. Exit codes: `0` ok, `2` validation or usage error, `3` report already exists, `4` could not lock `CLAUDE.md` (the report was still written).

---

## Files

| Path | Purpose |
|------|---------|
| `SKILL.md` | The skill definition Claude follows |
| `scripts/handoff_report.py` | Pydantic validator, report writer and `CLAUDE.md` indexer |
| `references/annotation-form.md` | Every form field, enum value and default |
| `examples/example-form.json` | A complete multi-entry form to copy |

To try the validator without Claude, run it on the example in a throwaway folder (it writes `.claude/reports/` and `.claude/CLAUDE.md` there):

```bash
mkdir -p /tmp/handoff-try && cd /tmp/handoff-try
cp ~/.claude/skills/handoff/examples/example-form.json form.json
uv run ~/.claude/skills/handoff/scripts/handoff_report.py --form form.json
```
