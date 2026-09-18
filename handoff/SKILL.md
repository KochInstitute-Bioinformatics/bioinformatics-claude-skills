---
name: handoff
description: This skill should be used when the user types "/handoff", or asks to "write a handoff", "create a session report", "log a handoff note", "record what happened this session", "hand off to the next session", "write up this session", or at the end of substantive work that needs a durable cross-session report. It captures the session as typed structured annotations (session narrative plus any decisions, risks, blockers, open questions, handoff notes, user preferences/corrections) in the Structured Reporting System (SRS) annotation format, validates them with a bundled pydantic script, writes a JSON report under <project>/.claude/reports/, and atomically adds a one-line pointer to the project CLAUDE.md index.
version: 0.1.0
---

# Handoff — structured session reports

## Purpose

Replace free-form Markdown session reports with **deterministic, schema-validated structured
reports**. `/handoff` is a present-day stand-in for the future Structured Reporting System
(SRS) `/narrative` command: it produces an *SRS-shaped* artifact — an *annotation file* of
typed, evidence-bearing entries — instead of prose. The annotation *envelope* matches SRS
field-for-field; evidence refs are a documented relaxation SRS would migrate on ingest (see
"Evidence" below). Filling out a structured form (not writing Markdown) and validating it with
a real schema makes reports consistent now and far cheaper for SRS to ingest later.

The flow: **gather the session → compose a form (typed annotations) → run the validator
script → it writes the JSON report and atomically indexes it in CLAUDE.md.** The agent never
hand-writes the report file; the script generates it from the validated form.

## When to use

Use at the end of substantive work (a build, postmortem, design decision, workflow change,
incident, or handoff) or whenever the user invokes `/handoff`. Skip for trivial one-off
edits. This is the project's reporting mechanism — it supersedes the older "write a timestamped
Markdown report + add a pointer" instruction.

## Core model (read this before composing)

- **One report = one annotation file = many typed annotations.** Always include at least one
  `session_narrative`. Add a separate annotation for each distinct decision, risk, blocker,
  open question, handoff note, user correction, or user preference. SRS builds its
  decision/risk/open-question views *only* from these typed entries — so decompose, don't
  cram everything into one narrative.
- **Evidence is lightweight here.** Each annotation *may* carry `evidence_refs`, but they are
  optional and use a simple `ref` string (`src/foo.rs:42-50`, a commit SHA, a command, a URL)
  — **not** transcript byte-offsets or hashes. This is the one real divergence from the SRS
  schema (SRS evidence is transcript-anchored and `ref` is not an SRS field); SRS ingestion
  re-derives transcript anchors. See `references/annotation-form.md` for the exact divergences.
- **User-sourced entries are authored by the user.** For `user_correction`, `user_preference`,
  or a user-stated `decision`, set `created_by: user` (it defaults to `agent`).
- **Verifiable claims only.** Put only what was observed or can be pointed to (file:line,
  command output, a verbatim quote) in `details`. Mark inferences as inferences. Set
  `verification_status` honestly (`verified` only when actually confirmed, e.g. tests passed)
  and `must_verify_before_use: true` for anything risky-but-unconfirmed.
- **Reports vs memory.** Reports = dated narrative tied to an event (this skill). Atemporal
  facts (profile, preferences, recurring feedback) go to `~/.claude/projects/.../memory/`, not
  here.

## Workflow

### 1. Gather the session

Collect what changed and why: files/symbols touched, commands run, commits, decisions and
their rationale, risks, blockers, open questions, and what the next session needs. Note the
verification status of each claim (was it actually run/tested, or asserted?).

### 2. Decompose into typed annotations

Pick an `annotation_type` per entry from this table (full guidance in
`references/annotation-form.md`):

| Type | Use for |
|---|---|
| `session_narrative` | the overall "what happened this session" (always include one) |
| `decision` / `inferred_decision` | a choice made (user-stated) / inferred from actions |
| `design_pivot` | a change in direction/approach |
| `risk` / `blocker` | a hazard / something stopping progress |
| `open_question` | unresolved question for the next session |
| `handoff_note` | "resume here" / next-step pointer |
| `user_correction` / `user_preference` | a correction or stated preference from the user |
| `hypothesis` / `implementation_note` / `retrieval_hint` / `warning` / `other` | as named |

Each annotation needs only `annotation_type` + `blocks.summary` + `blocks.details`.
Everything else has sensible defaults (see the field reference). Add `blocks.agent_guidance`
("what the next agent should do") and `blocks.caveats` where useful.

### 3. Write the form

Use the Write tool to create the form at `<project>/.claude/.handoff-form.json`. Minimum
viable form:

```json
{
  "report_meta": {
    "slug": "short-kebab-slug",
    "gist": "One-line summary for the CLAUDE.md index."
  },
  "annotations": [
    {"annotation_type": "session_narrative",
     "blocks": {"summary": "…", "details": "…", "agent_guidance": "…"}}
  ]
}
```

- `slug` is auto-normalized to kebab-case; the final filename is `YYYY-MM-DD-<slug>.json`.
- `gist` becomes the index pointer text. Keep it to one or two sentences.
- `report_meta` also accepts optional `title`, `date` (YYYY-MM-DD, defaults to today), and
  `session_id`.
- See `examples/example-form.json` for a complete multi-annotation example, and
  `references/annotation-form.md` for every field, enum value, and the SRS forward-compat
  mapping.

### 4. Run the validator

```bash
uv run ~/.claude/skills/handoff/scripts/handoff_report.py --form <project>/.claude/.handoff-form.json
```

By default the script: auto-detects the project root from the cwd, validates the form with
pydantic, writes `<root>/.claude/reports/YYYY-MM-DD-<slug>.json`, atomically inserts the
pointer at the **top** of the `## Session reports / handoffs` list in the project CLAUDE.md
(creating the section/file if absent), and deletes the input form on success. The CLAUDE.md
update is guarded by an advisory `flock` + atomic rename (POSIX-only), so concurrent **/handoff
runs** cannot corrupt it — a human editing CLAUDE.md by hand is not serialized against it. The
script prints the resolved absolute report path and CLAUDE.md so a wrong-tree write is visible.

**Safety:** the script refuses (exit 2) to auto-target your home directory as the project root
or to auto-index into your global `~/.claude/CLAUDE.md`. Run it from inside a project, or pass
`--root` / `--claude-md` explicitly if targeting an unusual location is genuinely intended.

Useful flags:

| Flag | Effect |
|---|---|
| `--claude-md <path>` | index into a specific CLAUDE.md (e.g. a different scope) instead of the auto-detected project one |
| `--no-index` | write the report only; touch no CLAUDE.md |
| `--root <dir>` | set the project root explicitly |
| `--force` | overwrite an existing report of the same name (default: refuse, exit 3) |
| `--keep-form` | keep the input form file after success |

Exit codes: `0` ok · `2` validation/usage error · `3` report already exists · `4` could not
acquire the CLAUDE.md lock.

### 5. Confirm and report back

On success the script prints the report path, the CLAUDE.md it indexed into, and the exact
`INDEX_LINE:`. Relay the report path and index line to the user.

- **Exit 2** (validation/usage or a safety refusal): read the message, fix the form (or pass
  `--root`/`--claude-md`), and rerun.
- **Exit 3** (report of this date+slug already exists): pick a more specific slug — do not
  blindly `--force`, and never silently overwrite a prior report. `--force` also applies even
  with `--no-index`.
- **Exit 4** (could not lock CLAUDE.md): the report *was* written; only indexing failed. Add
  the printed `INDEX_LINE` to the top of `## Session reports / handoffs` yourself, **or** re-run
  with `--force` (this recovery case is the one exception to the no-`--force` rule).

## Onboarding behavior this enables

Reports are JSON. On session start, after reading CLAUDE.md, **read the most recent report(s)
in full** (the JSON `blocks.summary` / `blocks.details` / `blocks.agent_guidance` of each
annotation) — do not skim. The index in CLAUDE.md is only pointers; the report is the
cross-session context.

## Additional resources

- **`references/annotation-form.md`** — complete field-by-field form reference: every enum
  value with guidance, the evidence relaxation, defaults, slug/gist conventions, the output
  JSON shape, and the mapping to the SRS `annotation.schema.json`.
- **`examples/example-form.json`** — a complete, realistic multi-annotation form to copy and
  adapt.
- **`scripts/handoff_report.py`** — the pydantic validator + atomic report writer + atomic
  CLAUDE.md indexer (PEP 723 inline deps; run with `uv run`).
