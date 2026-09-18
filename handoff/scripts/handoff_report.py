# /// script
# requires-python = ">=3.10"
# dependencies = ["pydantic>=2,<3"]
# ///
"""
handoff_report.py — validate a /handoff form and emit a structured session report.

Reads a "form" JSON (report_meta + a list of typed annotations), validates it with
pydantic against a relaxed, SRS-shaped variant of the Structured Reporting
System (SRS) annotation schema, writes the validated report to
``<project>/.claude/reports/YYYY-MM-DD-<slug>.json`` (atomic; never clobbers an existing
report), and inserts a one-line pointer at the TOP of the "## Session reports / handoffs"
list in the project CLAUDE.md (atomic + advisory ``fcntl.flock`` so two concurrent runs of
this script cannot corrupt the file; POSIX-only, advisory — a non-handoff editor is not
serialized against it).

The on-disk report keeps an SRS-shaped ``AnnotationFile`` payload under the
``annotation_file`` key, wrapped with a thin ``report_meta`` block carrying date / slug /
gist used for the index pointer.

Relaxations versus the SRS annotation schema (a migration step is required before this
validates against the strict SRS schema — it is SRS-*shaped*, not byte-identical):
  - ``evidence_refs`` is OPTIONAL and may be empty (SRS requires >= 1 per annotation).
  - each evidence entry is a lightweight pointer (``evidence_type`` in a wider set; a
    ``ref`` string such as ``src/foo.rs:42-50`` / commit SHA / command / URL; an optional
    string ``excerpt``) instead of SRS's transcript-only ref (``evidence_type`` const
    ``transcript``, a required ``anchor`` with byte offsets + SHA-256, and an object
    ``excerpt``). SRS ingestion must re-derive transcript anchors; the ``ref``/string
    ``excerpt`` fields do not exist in the SRS ``EvidenceRef`` and would be dropped.
The annotation *envelope* (annotation_type, interpretation_source, created_by.actor_type,
created_at, blocks, verification_status, must_verify_before_use, supersedes_annotation_ids,
claims) does match SRS field-for-field.

Exit codes: 0 ok · 2 validation/usage error · 3 report already exists (no --force) ·
4 could not acquire the CLAUDE.md lock.

Run via: uv run handoff_report.py --form FORM.json [flags]
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
import re
import sys
import tempfile
import time
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import Literal, NoReturn, Optional

from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator

SCHEMA_VERSION = "handoff/v1"
SECTION_HEADER = "## Session reports / handoffs"
DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")
ID_RE = {
    "annotation": re.compile(r"^ANN-[A-Za-z0-9][A-Za-z0-9_.-]*$"),
    "evidence": re.compile(r"^EV-[A-Za-z0-9][A-Za-z0-9_.-]*$"),
    "claim": re.compile(r"^CLM-[A-Za-z0-9][A-Za-z0-9_.-]*$"),
}

# --- enums (mirror the SRS annotation schema; EvidenceType is a relaxed superset) ---------
ActorType = Literal["user", "agent", "tool", "system", "unknown"]
VerificationStatus = Literal[
    "verified", "unverified", "partially_verified", "superseded", "contradicted"
]
AnnotationType = Literal[
    "session_narrative", "user_correction", "design_pivot", "decision",
    "inferred_decision", "risk", "blocker", "open_question", "hypothesis",
    "implementation_note", "retrieval_hint", "warning", "user_preference",
    "handoff_note", "other",
]
InterpretationSource = Literal[
    "user_stated", "agent_suggested", "tool_derived", "mixed", "other"
]
EvidenceType = Literal[
    "transcript", "file", "command", "commit", "artifact", "url", "other"
]


# --- models -------------------------------------------------------------------------------
class TranscriptAnchor(BaseModel):
    """Full SRS transcript anchor — optional here, accepted for forward-compatibility."""

    model_config = ConfigDict(extra="forbid")
    transcript_path: str
    line_start: int = Field(ge=1)
    line_end: int = Field(ge=1)
    byte_start: int = Field(ge=0)
    byte_end: int = Field(ge=0)
    source_record_hash: str = Field(min_length=1)
    hash_algorithm: Literal["sha256"] = "sha256"


class EvidenceRef(BaseModel):
    model_config = ConfigDict(extra="forbid")
    evidence_id: Optional[str] = None  # auto-assigned EV-<ann>.<n> if absent
    evidence_type: EvidenceType = "file"
    ref: str = Field(min_length=1)  # lightweight pointer
    excerpt: Optional[str] = None  # small copied text (optional)
    anchor: Optional[TranscriptAnchor] = None  # full SRS anchor (optional)

    @field_validator("evidence_id")
    @classmethod
    def _check_evidence_id(cls, v: Optional[str]) -> Optional[str]:
        if v is not None and not ID_RE["evidence"].match(v):
            raise ValueError("evidence_id must match ^EV-[A-Za-z0-9][A-Za-z0-9_.-]*$")
        return v


class AnnotationBlocks(BaseModel):
    model_config = ConfigDict(extra="forbid")
    summary: str = Field(min_length=1)
    details: str = Field(min_length=1)
    agent_guidance: Optional[str] = None
    caveats: Optional[str] = None


class Claim(BaseModel):
    """Relaxed SRS claim: an atomic subject-predicate-object triple."""

    model_config = ConfigDict(extra="forbid")
    claim_id: Optional[str] = None  # auto-assigned CLM-<ann>.<n> if absent
    subject: str = Field(min_length=1)
    predicate: str = Field(min_length=1)
    object: str = Field(min_length=1)
    verification_status: VerificationStatus = "unverified"

    @field_validator("claim_id")
    @classmethod
    def _check_claim_id(cls, v: Optional[str]) -> Optional[str]:
        if v is not None and not ID_RE["claim"].match(v):
            raise ValueError("claim_id must match ^CLM-[A-Za-z0-9][A-Za-z0-9_.-]*$")
        return v


class Annotation(BaseModel):
    model_config = ConfigDict(extra="forbid")
    annotation_id: Optional[str] = None  # auto-assigned ANN-<n> if absent
    annotation_type: AnnotationType
    interpretation_source: InterpretationSource = "mixed"
    created_by: ActorType = "agent"  # serialized as {"actor_type": ...} on output
    created_at: Optional[str] = None  # ISO 8601; defaults to generation time
    blocks: AnnotationBlocks
    evidence_refs: list[EvidenceRef] = Field(default_factory=list)
    verification_status: VerificationStatus = "unverified"
    must_verify_before_use: bool = False
    supersedes_annotation_ids: list[str] = Field(default_factory=list)
    claims: list[Claim] = Field(default_factory=list)

    @field_validator("annotation_id")
    @classmethod
    def _check_annotation_id(cls, v: Optional[str]) -> Optional[str]:
        if v is not None and not ID_RE["annotation"].match(v):
            raise ValueError("annotation_id must match ^ANN-[A-Za-z0-9][A-Za-z0-9_.-]*$")
        return v

    @field_validator("created_by", mode="before")
    @classmethod
    def _normalize_created_by(cls, v):
        # Accept either a bare actor string or the SRS nested {"actor_type": ...} form.
        if isinstance(v, dict):
            return v.get("actor_type", "agent")
        return v


class ReportMeta(BaseModel):
    model_config = ConfigDict(extra="forbid")
    slug: str
    gist: str = Field(min_length=1)
    title: Optional[str] = None
    date: Optional[str] = None  # YYYY-MM-DD; defaults to today
    session_id: Optional[str] = None

    @field_validator("slug")
    @classmethod
    def _normalize_slug(cls, v: str) -> str:
        s = re.sub(r"[^a-z0-9]+", "-", v.strip().lower())
        s = re.sub(r"-+", "-", s).strip("-")
        if not s:
            raise ValueError("slug is empty after normalization")
        return s

    @field_validator("date")
    @classmethod
    def _check_date(cls, v: Optional[str]) -> Optional[str]:
        if v is not None and not DATE_RE.match(v):
            raise ValueError("date must be YYYY-MM-DD")
        return v


class HandoffForm(BaseModel):
    model_config = ConfigDict(extra="forbid")
    report_meta: ReportMeta
    annotations: list[Annotation] = Field(min_length=1)

    @field_validator("annotations")
    @classmethod
    def _unique_ids(cls, v: list[Annotation]) -> list[Annotation]:
        seen: set[str] = set()
        for a in v:
            if a.annotation_id is None:
                continue
            if a.annotation_id in seen:
                raise ValueError(f"duplicate annotation_id: {a.annotation_id}")
            seen.add(a.annotation_id)
        return v


# --- helpers ------------------------------------------------------------------------------
def die(code: int, msg: str) -> NoReturn:
    print(f"handoff_report: {msg}", file=sys.stderr)
    sys.exit(code)


def today() -> str:
    return datetime.now().strftime("%Y-%m-%d")


def now_iso() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def find_root(start: Path) -> Path:
    cur = start.resolve()
    for p in [cur, *cur.parents]:
        if (p / ".git").exists() or (p / ".claude").is_dir():
            return p
    return cur


def lock_path_for(target: Path) -> Path:
    """A stable lock file in the system temp dir, keyed by the target's absolute path.

    Kept out of the repository so it never shows up in git status, while concurrent writers
    targeting the same CLAUDE.md still share one lock (same absolute path -> same hash).
    """
    digest = hashlib.sha256(str(target.resolve()).encode("utf-8")).hexdigest()[:16]
    return Path(tempfile.gettempdir()) / f"handoff-claude-md-{digest}.lock"


def pick_claude_md(root: Path, override: Optional[str]) -> Path:
    if override:
        return Path(override).expanduser().resolve()
    nested = root / ".claude" / "CLAUDE.md"
    flat = root / "CLAUDE.md"
    if nested.exists():
        return nested
    if flat.exists():
        return flat
    return nested  # default creation location


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix=".tmp-handoff-", suffix=path.suffix)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(text)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


@contextmanager
def file_lock(lock_path: Path, timeout: float = 15.0):
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(str(lock_path), os.O_CREAT | os.O_RDWR, 0o644)
    start = time.monotonic()
    try:
        while True:
            try:
                fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except OSError:
                if time.monotonic() - start > timeout:
                    raise TimeoutError(f"could not acquire lock {lock_path} within {timeout}s")
                time.sleep(0.1)
        yield
    finally:
        try:
            fcntl.flock(fd, fcntl.LOCK_UN)
        finally:
            os.close(fd)


def insert_index_line(content: str, new_line: str) -> str:
    """Insert ``new_line`` at the top of the bullet list under SECTION_HEADER.

    Newest-first ordering. Creates the section (appended) if it does not exist.
    """
    had_trailing_nl = content.endswith("\n")
    lines = content.splitlines()

    header_idx = None
    for i, ln in enumerate(lines):
        if ln.strip().lower().startswith("## session reports"):
            header_idx = i
            break

    if header_idx is None:
        prefix = content if content.endswith("\n") or content == "" else content + "\n"
        return f"{prefix}\n{SECTION_HEADER}\n\n{new_line}\n"

    sec_end = len(lines)
    for j in range(header_idx + 1, len(lines)):
        if lines[j].startswith("## "):
            sec_end = j
            break

    first_bullet = None
    for j in range(header_idx + 1, sec_end):
        if lines[j].lstrip().startswith("- "):
            first_bullet = j
            break

    if first_bullet is not None:
        lines.insert(first_bullet, new_line)
    else:
        ins = sec_end
        while ins - 1 > header_idx and lines[ins - 1].strip() == "":
            ins -= 1
        block = []
        if ins - 1 >= header_idx and lines[ins - 1].strip() != "":
            block.append("")
        block.append(new_line)
        lines[ins:ins] = block

    result = "\n".join(lines)
    if had_trailing_nl:
        result += "\n"
    return result


def clean_gist(gist: str) -> str:
    return gist.strip().rstrip(".").strip()


def build_output(form: HandoffForm, date: str, generated_at: str, rel_report: str,
                 claude_md_rel: Optional[str]) -> dict:
    annotations_out = []
    for ai, a in enumerate(form.annotations, start=1):
        aid = a.annotation_id or f"ANN-{ai}"
        evidence_out = []
        for ei, ev in enumerate(a.evidence_refs, start=1):
            entry = {
                "evidence_id": ev.evidence_id or f"EV-{ai}.{ei}",
                "evidence_type": ev.evidence_type,
                "ref": ev.ref,
            }
            if ev.excerpt is not None:
                entry["excerpt"] = ev.excerpt
            if ev.anchor is not None:
                entry["anchor"] = ev.anchor.model_dump()
            evidence_out.append(entry)
        claims_out = []
        for ci, c in enumerate(a.claims, start=1):
            claims_out.append({
                "claim_id": c.claim_id or f"CLM-{ai}.{ci}",
                "subject": c.subject,
                "predicate": c.predicate,
                "object": c.object,
                "verification_status": c.verification_status,
            })
        annotations_out.append({
            "annotation_id": aid,
            "annotation_type": a.annotation_type,
            "interpretation_source": a.interpretation_source,
            "created_by": {"actor_type": a.created_by},
            "created_at": a.created_at or generated_at,
            "blocks": {
                "summary": a.blocks.summary,
                "details": a.blocks.details,
                "agent_guidance": a.blocks.agent_guidance,
                "caveats": a.blocks.caveats,
            },
            "evidence_refs": evidence_out,
            "verification_status": a.verification_status,
            "must_verify_before_use": a.must_verify_before_use,
            "supersedes_annotation_ids": a.supersedes_annotation_ids,
            "claims": claims_out,
        })

    return {
        "report_meta": {
            "schema_version": SCHEMA_VERSION,
            "tool": "handoff",
            "date": date,
            "slug": form.report_meta.slug,
            "title": form.report_meta.title or form.report_meta.gist,
            "gist": form.report_meta.gist,
            "generated_at": generated_at,
            "report_path": rel_report,
            "indexed_in": claude_md_rel,
        },
        "annotation_file": {
            "schema_version": SCHEMA_VERSION,
            "session_id": form.report_meta.session_id or f"handoff-{date}-{form.report_meta.slug}",
            "annotations": annotations_out,
        },
    }


def read_form(form_arg: str) -> str:
    if form_arg == "-":
        return sys.stdin.read()
    p = Path(form_arg).expanduser()
    if not p.is_file():
        die(2, f"form file not found: {form_arg}")
    return p.read_text(encoding="utf-8")


# --- main ---------------------------------------------------------------------------------
def parse_args(argv=None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Validate a /handoff form and write a structured report.")
    ap.add_argument("--form", required=True, help="path to the form JSON, or '-' for stdin")
    ap.add_argument("--root", default=None, help="project root (default: auto-detect from cwd)")
    ap.add_argument("--claude-md", default=None, help="explicit CLAUDE.md to index into")
    ap.add_argument("--no-index", action="store_true", help="write the report but do not touch any CLAUDE.md")
    ap.add_argument("--force", action="store_true", help="overwrite an existing report with the same name")
    ap.add_argument("--keep-form", action="store_true", help="do not delete the input form file on success")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)

    raw = read_form(args.form)
    try:
        data = json.loads(raw)
    except json.JSONDecodeError as e:
        die(2, f"form is not valid JSON: {e}")
    try:
        form = HandoffForm.model_validate(data)
    except ValidationError as e:
        die(2, "form failed schema validation:\n" + str(e))

    root = Path(args.root).expanduser().resolve() if args.root else find_root(Path.cwd())
    home = Path.home().resolve()
    if args.root is None and root == home:
        die(2, "auto-detected project root is your home directory ($HOME); refusing to write\n"
               "  a session report there. cd into a project, or pass --root explicitly.")
    date = form.report_meta.date or today()
    slug = form.report_meta.slug
    filename = f"{date}-{slug}.json"
    report_path = root / ".claude" / "reports" / filename
    rel_report = f".claude/reports/{filename}"

    if report_path.exists() and not args.force:
        die(3, f"report already exists: {report_path}\n  choose a different slug or pass --force")

    claude_md = None
    claude_md_rel = None
    if not args.no_index:
        claude_md = pick_claude_md(root, args.claude_md)
        if args.claude_md is None and claude_md.resolve() == (home / ".claude" / "CLAUDE.md").resolve():
            die(2, "auto-detected CLAUDE.md is your global ~/.claude/CLAUDE.md; refusing to index\n"
                   "  session reports into your global instructions. Pass --claude-md explicitly if\n"
                   "  that is truly intended, or pass --no-index.")
        try:
            claude_md_rel = os.path.relpath(claude_md, root)
        except ValueError:
            claude_md_rel = str(claude_md)

    generated_at = now_iso()
    out = build_output(form, date, generated_at, rel_report, claude_md_rel)
    atomic_write(report_path, json.dumps(out, indent=2, ensure_ascii=False) + "\n")

    index_line = None
    if not args.no_index:
        assert claude_md is not None  # set above whenever indexing is enabled
        index_line = f"- {date} — {clean_gist(form.report_meta.gist)}. See `{rel_report}`."
        lock = lock_path_for(claude_md)
        try:
            with file_lock(lock):
                existing = claude_md.read_text(encoding="utf-8") if claude_md.exists() else "# CLAUDE.md\n"
                atomic_write(claude_md, insert_index_line(existing, index_line))
        except TimeoutError as e:
            print(f"handoff_report: {e}", file=sys.stderr)
            print(f"  the report WAS written to {report_path}, but indexing did not complete.", file=sys.stderr)
            print(f"  recover WITHOUT overwriting: add this line at the top of the", file=sys.stderr)
            print(f"  '## Session reports / handoffs' list in {claude_md}:", file=sys.stderr)
            print(f"      {index_line}", file=sys.stderr)
            print(f"  (or re-run with --force to retry indexing; it rewrites the identical report.)", file=sys.stderr)
            sys.exit(4)

    if args.form != "-" and not args.keep_form:
        try:
            os.unlink(args.form)
        except OSError:
            pass

    print(f"OK: wrote {report_path} ({len(form.annotations)} annotation(s))")
    if index_line:
        print(f"OK: indexed in {claude_md}")
        print(f"INDEX_LINE: {index_line}")
    else:
        print("OK: indexing skipped (--no-index)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
