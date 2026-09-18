# Handoff form — field reference

This is the complete reference for the `/handoff` form. The authoritative validator is
`scripts/handoff_report.py` (pydantic). The form mirrors the Structured Reporting System (SRS)
annotation schema (`structured-reporting-tool/_unpacked/schema_draft/annotation.schema.json`)
with exactly one relaxation: evidence is optional and lightweight (see "Evidence" below).

## Top-level form shape

```jsonc
{
  "report_meta": { ... },      // required: metadata for filename + index pointer
  "annotations": [ { ... } ]   // required: >= 1 typed annotation
}
```

The validator writes a report whose `annotation_file` block is the SRS-faithful payload, and
wraps it with a `report_meta` block. The agent only writes the *form*; the script produces the
report.

## `report_meta`

| Field | Required | Default | Notes |
|---|---|---|---|
| `slug` | yes | — | normalized to kebab-case; filename is `YYYY-MM-DD-<slug>.json` |
| `gist` | yes | — | one/two sentences; becomes the CLAUDE.md pointer text |
| `title` | no | = `gist` | human title stored in the report |
| `date` | no | today | `YYYY-MM-DD`; prepended to the filename and used in the pointer |
| `session_id` | no | `handoff-<date>-<slug>` | SRS session id for the annotation file |

The index pointer is rendered as:
`- <date> — <gist>. See \`.claude/reports/<date>-<slug>.json\`.`
and inserted at the **top** of the `## Session reports / handoffs` list (newest first).

## `annotations[]` — one entry per distinct thing

Minimum per annotation: `annotation_type`, `blocks.summary`, `blocks.details`. All other
fields default.

| Field | Required | Default | Notes |
|---|---|---|---|
| `annotation_type` | yes | — | enum, see below |
| `blocks.summary` | yes | — | one-line headline |
| `blocks.details` | yes | — | the substance; verifiable claims only |
| `blocks.agent_guidance` | no | null | what the next agent should do |
| `blocks.caveats` | no | null | warnings / things not to trust |
| `interpretation_source` | no | `mixed` | enum, see below |
| `created_by` | no | `agent` | actor string, or `{"actor_type": "..."}` |
| `created_at` | no | generation time | ISO 8601 if provided |
| `verification_status` | no | `unverified` | enum, see below — be honest |
| `must_verify_before_use` | no | `false` | `true` for risky/unconfirmed entries |
| `evidence_refs` | no | `[]` | lightweight pointers, see below |
| `supersedes_annotation_ids` | no | `[]` | `ANN-…` ids this entry replaces |
| `claims` | no | `[]` | atomic subject/predicate/object triples |
| `annotation_id` | no | `ANN-<n>` | auto-assigned in order if omitted |

### `annotation_type` (enum)

`session_narrative` (always include one), `user_correction`, `design_pivot`, `decision`,
`inferred_decision`, `risk`, `blocker`, `open_question`, `hypothesis`, `implementation_note`,
`retrieval_hint`, `warning`, `user_preference`, `handoff_note`, `other`.

Guidance: split the session into one `session_narrative` plus one annotation per decision /
risk / blocker / open question / handoff note. Use `decision` when the user chose it,
`inferred_decision` when you inferred it from actions. Use `user_correction` /
`user_preference` to capture what the user told you, with `interpretation_source: user_stated`.

### `interpretation_source` (enum)

`user_stated` (the user said it), `agent_suggested` (you proposed/observed it), `tool_derived`
(from tool/command output), `mixed` (a blend — the common default), `other`.

### `verification_status` (enum)

`verified` (the annotation's content is confirmed — tests passed, output observed, or it is a
fact the user directly stated, e.g. a `decision`/`user_preference` with
`interpretation_source: user_stated`), `unverified` (default — an agent claim not yet
confirmed), `partially_verified`, `superseded` (replaced by a later entry), `contradicted`
(later evidence contradicts it). Do not mark an agent-derived claim `verified` unless it truly
was confirmed.

### `created_by` actor types

`user`, `agent`, `tool`, `system`, `unknown`. Accepts a bare string (`"user"`) or the SRS
nested form (`{"actor_type": "user"}`); both serialize to the nested form on output. Defaults
to `agent` — **override to `user` for `user_correction`, `user_preference`, and user-stated
`decision` entries**, otherwise the user's own words are misattributed to the agent.

## Evidence (`evidence_refs[]`) — the one relaxation

SRS requires every annotation to carry >= 1 evidence ref with a full transcript anchor (line +
byte offsets + a SHA-256 of the source record). Those values are produced deterministically by
the future SRS Rust parser and cannot be hand-authored reliably, so here evidence is
**optional and lightweight**:

```jsonc
{
  "evidence_type": "file",                // file | command | commit | transcript | artifact | url | other
  "ref": "src/foo.rs:42-50",              // required: the lightweight pointer
  "excerpt": "optional small copied text",
  "anchor": { ... }                        // optional full SRS TranscriptAnchor, if you have it
}
```

- `evidence_id` auto-assigns to `EV-<annotationIndex>.<n>` if omitted.
- Prefer concrete refs: `path:line`, a commit SHA, a command line, a URL, an artifact path.

**This shape is SRS-*shaped*, not SRS-valid — it would not pass the strict SRS `EvidenceRef`
as-is.** The concrete divergences a future SRS ingest must handle:

| Here (handoff) | SRS `EvidenceRef` (`annotation.schema.json:111-136`) |
|---|---|
| `evidence_type` ∈ 7-value set | `const: "transcript"` only |
| `ref` string (the pointer) | no `ref` field; `additionalProperties: false` → would be dropped |
| `excerpt` is a bare string | `excerpt` is an object `{text, excerpt_hash, hash_algorithm}` |
| `anchor` optional | `anchor` required (full byte+sha256 `TranscriptAnchor`) |
| `evidence_refs` may be empty | `minItems: 1` per annotation |

So SRS ingestion re-derives transcript anchors from the source transcript; the `ref`/string-
`excerpt` are convenience fields, not a drop-in upgrade path. The annotation *envelope* around
evidence does match SRS field-for-field.

## `claims[]` (optional, advanced)

Atomic factual triples that downstream tooling can index:

```jsonc
{ "subject": "report", "predicate": "format_is", "object": "json",
  "verification_status": "verified" }
```

`claim_id` auto-assigns to `CLM-<annotationIndex>.<n>`. Keep `subject`/`predicate`/`object`
short and atomic. The report stores only the triple (no human-readable claim text); a consumer
can render prose from the three fields.

## Output report shape (what the script writes)

```jsonc
{
  "report_meta": { "schema_version": "handoff/v1", "tool": "handoff", "date", "slug",
                   "title", "gist", "generated_at", "report_path", "indexed_in" },
  "annotation_file": {                     // <- SRS-shaped AnnotationFile payload
    "schema_version": "handoff/v1",
    "session_id": "...",
    "annotations": [ /* SRS-shaped annotations with auto IDs + nested created_by */ ]
  }
}
```

When SRS exists, `annotation_file` is the part it ingests (after the evidence migration noted
above); `report_meta` is handoff-tool bookkeeping.

## Conventions recap

- Filename: `<project>/.claude/reports/YYYY-MM-DD-<slug>.json`. The `slug` is normalized to
  kebab-case, so two visually-distinct slugs can collapse to the same filename (and then collide
  at exit 3); the printed report path shows the normalized form.
- Index: top of `## Session reports / handoffs` in the project CLAUDE.md, newest first. The
  index write uses an advisory `fcntl.flock` (POSIX-only) + atomic rename, so it is safe against
  other concurrent `/handoff` runs — but not against a human editing CLAUDE.md by hand.
- Safety: the script refuses (exit 2) to auto-target `$HOME` as root or your global
  `~/.claude/CLAUDE.md`; pass `--root`/`--claude-md` to override intentionally.
- Never overwrite an existing report (the script refuses with exit 3 unless `--force`).
- Onboarding: read the most recent report(s) **in full**, not skimmed.
