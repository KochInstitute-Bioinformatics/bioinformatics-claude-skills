# nfcore-rnavar-setup: all-in params file — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the defects found by the first real run of `nfcore-rnavar-setup` (numeric CLI params rejected as strings under Nextflow 26.04; empty reports/trace on re-run; `params.max_*` warning; wrong output-dir names; STAR build on tiny genomes; no custom-reference option) by moving every pipeline parameter into a generated params file.

**Architecture:** The skill (`nfcore-rnavar-setup/nfcore-rnavar-setup.md`) generates `{PARAMS_YAML}` containing every rnavar parameter (paths/strings quoted, numbers/booleans unquoted) and a launch script whose `nextflow run` line carries only `-params-file`. The static checker (`nfcore-rnavar-setup/tests/check_skill.sh`) gains two rules that would have caught the defect: rnavar parameters must not appear as `--flag` on the launch command, and every key in the skill's params-file template must be a real 1.3.0 schema parameter.

**Tech Stack:** Markdown skill file, bash + grep + curl checker (no python on this host), nf-core/rnavar 1.3.0, Nextflow 26.04.6, YAML.

**Spec:** the approved in-chat design of 2026-09-28 (recorded here) extending `docs/superpowers/specs/2026-09-28-nfcore-rnavar-setup-design.md`. Evidence: the first live run in `/net/bmc-lab3/data/bcc/yannvrb/rnavar_test/` (Nextflow rejected `--read_length 151`: "Value is [string] but should be [number]"; fixed there by `-params-file`).

## Global Constraints

- Every key in the params file and every `--flag` in prose must be a parameter in the rnavar 1.3.0 `nextflow_schema.json` (or a documented allowlisted non-rnavar tool flag). Never emit `annotation_cache`, `gencode`, `strandedness`.
- The launch command is exactly `nextflow run nf-core/rnavar -r {VERSION} -c nextflow.config -profile slurm,singularity -params-file {PARAMS_YAML}` — no rnavar `--flag` on it.
- Params-file typing: paths and strings double-quoted; numbers and booleans unquoted (`read_length: 151`, `skip_baserecalibration: true`, `star_twopass: false`).
- Keys equal to the pipeline default are omitted. `read_length` is always present; `seq_platform: "illumina"` is present.
- Known sites: either the four `dbsnp`, `dbsnp_tbi`, `known_indels`, `known_indels_tbi` keys or `skip_baserecalibration: true`; the contig guard must enforce that the VCF contig names equal the FASTA's before the pipeline starts.
- `nextflow.config` is written only if absent, never overwritten. Existing `{PARAMS_YAML}` / launch script / samplesheet are never overwritten silently (ask: 1. overwrite · 2. choose another filename).
- `resourceLimits` in the config are literal values (cpus 16, memory '64 GB', time '24h'); no `params.max_*`.
- Helper jobs never run on the login node; no embedded unverified URLs; all finite-choice questions numbered.
- The skill file stays self-contained (no references to other skills' steps).
- Git: show `git status` and `git diff --stat` before each commit; never `git push` without explicit user approval. Commit messages end with the attribution line from the session system-reminder (`Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`).

## Review Focus

1. **Empty optional blocks in the params file** (no annotation, no variant options, BQSR skipped): the YAML stays valid; no blank/dangling keys. Owned by Task 1.
2. **String that looks numeric** (`multiqc_title` such as `260928` would be parsed as a number): titles are always double-quoted. Owned by Task 1.
3. **Small genome STAR build**: `--genomeSAindexNbases` never above 14 and never below 4; resources scale down for small genomes. Owned by Task 2.
4. **Custom reference with a non-Ensembl name** (helper naming, `{ENS_VERSION}` unbound): no unbound placeholder in helper names. Owned by Task 2.
5. **Re-running after a failed launch** reproduces reports/trace (`overwrite = true`). Owned by Task 3.
6. **Contig guard**: empty contig on either side fails; first-contig read cannot SIGPIPE; params keys point at post-guard filenames. Owned by Task 2.

## File Structure

| File | Responsibility |
|---|---|
| `nfcore-rnavar-setup/nfcore-rnavar-setup.md` | The skill prompt. Modified in Tasks 1–3. |
| `nfcore-rnavar-setup/tests/check_skill.sh` | Static checker. Modified in Tasks 1–3 (new rules in Task 1). |
| `nfcore-rnavar-setup/README.md` | User README. Modified in Task 3. |

Paths relative to `/net/bmc-lab3/data/bcc/yannvrb/Github_repos/bioinformatics-claude-skills`.

---

### Task 1: Params-file core (Steps 5–9, 11) and the two new checker rules

**Files:**
- Modify: `nfcore-rnavar-setup/tests/check_skill.sh`
- Modify: `nfcore-rnavar-setup/nfcore-rnavar-setup.md` (Steps 6 `--star_index` sentence, 7, 8, 9, 11, Notes)

**Interfaces:**
- Consumes: placeholders bound in Steps 1–6 (`{VERSION}`, `{SAMPLESHEET_CSV}`, `{FASTA_PATH}`, `{GTF_PATH}`, `{STAR_INDEX}`, `{READ_LENGTH}`, `{MULTIQC_TITLE}`, `{OUTDIR}`, `{USER_EMAIL}`, `{CONDA_ENV}`).
- Produces: `{PARAMS_YAML}` (= the samplesheet file name with `_samplesheet.csv` replaced by `_params.yaml`, same directory), and the fragments `{KNOWN_SITES_PARAMS}`, `{VARIANT_PARAMS}`, `{ANNOTATION_PARAMS}` (each: zero or more complete YAML `key: value` lines, no trailing backslashes). Tasks 2–3 rely on these names.

- [ ] **Step 1: Add the two new checker rules first (test-first)**

In `check_skill.sh`, before the final `[ $fail ...` line, add a block `# --- Params-file rules` containing:

```bash
# (1) The launch command must carry no rnavar --flag (params file only)
launch=$(grep 'nextflow run nf-core/rnavar' "$SKILL" | head -1)
if [ -z "$launch" ]; then echo "FAIL: no 'nextflow run nf-core/rnavar' launch line found"; fail=1
elif echo "$launch" | grep -qE ' --[A-Za-z_]'; then echo "FAIL: launch line carries a --flag (must use -params-file only): $launch"; fail=1
fi
need "-params-file {PARAMS_YAML}"

# (2) Every key in the params-file template (fenced yaml block) must be a schema parameter
yaml_keys=$(awk '/^```yaml/{f=1;next} /^```/{f=0} f' "$SKILL" | grep -oE '^[a-z_0-9]+:' | tr -d ':' | sort -u)
[ -n "$yaml_keys" ] || { echo "FAIL: no fenced yaml params template found"; fail=1; }
for k in $yaml_keys; do
  echo "$schema_names" | grep -qx -- "$k" || { echo "FAIL: params-file key not in rnavar schema: $k"; fail=1; }
done
# --- end Params-file rules
```
Run the checker. Expected: FAIL (no yaml block, launch line still has flags).

- [ ] **Step 2: Rewrite the skill text**

  - **Step 6** (`--star_index '{STAR_INDEX}' is always passed ...` sentence): say the `star_index` key is always written to the params file so rnavar never rebuilds the index.
  - **Step 7:** option 1 emits the keys `dbsnp`, `dbsnp_tbi`, `known_indels`, `known_indels_tbi` (double-quoted paths) instead of `--dbsnp ...` flags; option 2 emits `skip_baserecalibration: true`. Keep everything else in Step 7 (schema `.vcf.gz`/`.tbi` requirement, run-time URL resolution, mandatory contig-name check and its wording) byte-identical; only change how the result is stored: "Store the result as `{KNOWN_SITES_PARAMS}`: the four key lines for option 1, or the single line `skip_baserecalibration: true` for option 2."
  - **Step 8:** every emitted option becomes a params-file key: `remove_duplicates: true`, `star_twopass: false`, `gatk_hc_call_conf: <int>`, `gatk_vf_qd_filter: <number>`, `gatk_vf_fs_filter: <number>`, `gatk_vf_window_size: <int>`, `gatk_vf_cluster_size: <int>`, `skip_variantfiltration: true`, `generate_gvcf: true`, `bam_csi_index: true`, `save_align_intermeds: true`. Collect as `{VARIANT_PARAMS}`. Keep "omit values equal to the default".
  - **Step 9:** `tools: "snpeff"` / `"vep"` / `"merge"`, `snpeff_cache`, `vep_cache`, `snpeff_db`, `vep_genome`, `vep_species`, `vep_cache_version` as keys (quoted strings). Keep the no-silent-`download_cache` rule (write it as "`download_cache: true` only if the user explicitly chooses it after being warned"); keep the dash-less `annotation_cache` note. Collect as `{ANNOTATION_PARAMS}`.
  - **Step 11** (rename to "Generate the params file and the submission script"): write `{PARAMS_YAML}` (first check whether it exists; if so ask 1. overwrite · 2. choose another filename) from this template inside a fenced ```yaml block:

    ```yaml
    # nf-core/rnavar {VERSION} parameters — generated by /nfcore-rnavar-setup
    input: "{SAMPLESHEET_CSV}"
    outdir: "{OUTDIR}"
    multiqc_title: "{MULTIQC_TITLE}"
    fasta: "{FASTA_PATH}"
    gtf: "{GTF_PATH}"
    star_index: "{STAR_INDEX}"
    read_length: {READ_LENGTH}
    seq_platform: "illumina"
    {KNOWN_SITES_PARAMS}
    {VARIANT_PARAMS}
    {ANNOTATION_PARAMS}
    ```
    State the typing rules (paths and strings double-quoted — `multiqc_title` always quoted so a numeric-looking title stays a string; numbers and booleans unquoted) and that an empty fragment placeholder line is deleted entirely. State **why**: under Nextflow 26.04, values passed as `--flag value` reach the schema validator as strings (observed: `--read_length 151` rejected as "[string] but should be [number]"); a params file keeps YAML types. Then the submission script template: the existing `#SBATCH` header, module/conda lines, and

    ```
    nextflow run nf-core/rnavar -r {VERSION} -c nextflow.config -profile slurm,singularity -params-file {PARAMS_YAML}
    ```
    Keep the "check whether the script exists and ask" sentence, the "Script written / To submit" block, and the helper-order sentence (`sbatch --dependency=afterok:<helper_jobid> ...`). Remove the old "Each `{..._LINES}` placeholder ..." paragraph.
  - **Notes:** replace "Omit flags that equal the pipeline default" with the params-file wording; add "Never pass rnavar parameters on the `nextflow run` command line — always the params file"; keep the "never emit" list.
  - Search the whole skill for remaining `{KNOWN_SITES_LINES}`, `{VARIANT_LINES}`, `{ANNOTATION_LINES}` and `\` continuation wording and update to the `_PARAMS` names.

- [ ] **Step 3: Update the older checker lines that asserted removed text**

Replace (do not delete without an equivalent) every `need`/`forbid` that asserts the old CLI form: `need "--read_length {READ_LENGTH}"`, `need "--star_index '{STAR_INDEX}'"`, `need "already ends in"`, `need "never leave a bare"`, and the three `forbid '{X_LINES}\'` lines. Keep `need "nextflow run nf-core/rnavar -r {VERSION}"` and `forbid "{VERSION_ENS}"`. New equivalents: `need "read_length: {READ_LENGTH}"`, `need 'star_index: "{STAR_INDEX}"'`, `need "skip_baserecalibration: true"`, `need "{KNOWN_SITES_PARAMS}"`, `need "{VARIANT_PARAMS}"`, `need "{ANNOTATION_PARAMS}"`, a `need` for the quoting-rule phrase you actually write, `forbid "{KNOWN_SITES_LINES}"`, `forbid "{VARIANT_LINES}"`, `forbid "{ANNOTATION_LINES}"`. Keep prose coverage: for each option moved to a key, keep a `need` for the key name (e.g. `need "gatk_hc_call_conf"`).

- [ ] **Step 4: Run the checker; expected PASS**

Run: `bash nfcore-rnavar-setup/tests/check_skill.sh nfcore-rnavar-setup/nfcore-rnavar-setup.md`
Expected: `PASS`. Also test the two new rules on a temp copy under /tmp: add `--read_length 151` to the launch line and confirm `FAIL: launch line carries a --flag`; add a bogus key `bogus_key: 1` to the yaml block and confirm `FAIL: params-file key not in rnavar schema: bogus_key`; delete the copy.

- [ ] **Step 5: Commit**

```bash
git status && git diff --stat
git add nfcore-rnavar-setup
git commit -m "nfcore-rnavar-setup: pass all pipeline parameters via a params file; add checker rules"
```

---

### Task 2: Custom reference option (Step 6) and helper scripts (Step 12)

**Files:**
- Modify: `nfcore-rnavar-setup/tests/check_skill.sh`
- Modify: `nfcore-rnavar-setup/nfcore-rnavar-setup.md` (Step 6, Step 7 contig paragraph where it mentions the helper, Step 12)

**Interfaces:**
- Consumes: `{PARAMS_YAML}`, `{KNOWN_SITES_PARAMS}` (Task 1), `{GENOME_DIR}`, `{FASTA_PATH}`, `{GTF_PATH}`, `{SJDB_OVERHANG}`, `{ASSEMBLY}`, `{ENS_VERSION}` (existing).
- Produces: `{REF_TAG}` (helper-name tag: `{ASSEMBLY}_ens{ENS_VERSION}` for an Ensembl reference, or `custom_{WD_NAME}` for a custom one) used in helper file names; `{GENOME_LENGTH}` and `{SA_INDEX_NBASES}` in the STAR helper.

- [ ] **Step 1: Write failing checker lines**

Add block `# --- Task 2 (custom reference + helpers)` with `need` lines for the distinctive phrases you will write, at least: `need "Custom reference"`, `need "{REF_TAG}"`, `need "genomeSAindexNbases"`, `need "{SA_INDEX_NBASES}"`, `need "tabix -l"`, a `need` for the non-empty guard text you write, `need "post-guard"`, and `forbid "build_star_index_rnavar_{ASSEMBLY}_ens{ENS_VERSION}.sh"` (old fixed name). Run: expected FAIL.

- [ ] **Step 2: Step 6 — custom reference choice**

Right after the organism/base-directory questions in Step 6, add a numbered question: "1. Ensembl release in the standard folder (default) · 2. Custom reference — I already have a FASTA and GTF". Option 2: ask for the FASTA and GTF paths and the directory that will hold indexes (`{GENOME_DIR}`), skip the Ensembl download/version logic, set `{REF_TAG}` = `custom_{WD_NAME}`, still apply the per-read-length STAR index rule. For option 1, `{REF_TAG}` = `{ASSEMBLY}_ens{ENS_VERSION}`. The existing "Other organisms" bullet becomes "Other organisms: use option 2". Report the GTF source as before.

- [ ] **Step 3: Step 12 — helpers**

  - Helper file names use the tag: `build_star_index_rnavar_{REF_TAG}.sh`, `prepare_known_sites_{REF_TAG}.sh`, `prepare_annotation_cache_{TOOL}.sh`.
  - STAR helper: compute the genome length (`grep -v '^>' {FASTA_PATH} | tr -d '\n' | wc -c`) as `{GENOME_LENGTH}` and `{SA_INDEX_NBASES}` = `min(14, floor(log2({GENOME_LENGTH})/2 - 1))`, clamped to a minimum of 4; pass `--genomeSAindexNbases {SA_INDEX_NBASES}` after `--sjdbOverhang` (for a 40 kb genome this is 6; for GRCh38 it is 14). Scale the helper's SBATCH resources by genome size: genomes over 1 Gb `-n 8 --mem=64G -t 4:00:00`; genomes under 100 Mb `-n 4 --mem=8G -t 00:30:00`; in between `-n 8 --mem=32G -t 2:00:00` (these respect the HPC defaults of at most 64 G and 4 h).
  - Contig guard (in `prepare_known_sites_{REF_TAG}.sh`): read the first contig with `tabix -l FILE | head -n1` (no SIGPIPE), read the FASTA's first header as before, **fail with exit 1 if either contig is empty**, keep the rename/`exit 1` logic, and add a paragraph headed "post-guard filenames": the `dbsnp`/`known_indels` keys (and their `.tbi`) in `{PARAMS_YAML}` must point at the files that pass the guard (the renamed `NEW.vcf.gz` when a rename happened), so the wizard decides the final filenames before writing the params file (writing the `NEW` names up front when a rename is expected).
  - Step 7's last bullet: when VCFs exist and the wizard finds a mismatch, the wizard also generates `prepare_known_sites_{REF_TAG}.sh` (option (ii) is no longer "rely on the helper" — it generates it).
  - The Step 12 opening sentence with `-n 8 --mem=64G -t 8:00:00` is replaced by "resources scaled as below".

- [ ] **Step 4: Run the checker; expected PASS**

Run: `bash nfcore-rnavar-setup/tests/check_skill.sh nfcore-rnavar-setup/nfcore-rnavar-setup.md` → `PASS`.

- [ ] **Step 5: Commit**

```bash
git status && git diff --stat
git add nfcore-rnavar-setup
git commit -m "nfcore-rnavar-setup: custom reference option; small-genome STAR helper; robust contig guard"
```

---

### Task 3: Config, hand-off, docs (Steps 10, 13, Notes, README)

**Files:**
- Modify: `nfcore-rnavar-setup/tests/check_skill.sh`
- Modify: `nfcore-rnavar-setup/nfcore-rnavar-setup.md` (Step 10, Step 13)
- Modify: `nfcore-rnavar-setup/README.md`

**Interfaces:** consumes Task 1–2 names; produces nothing later tasks need.

- [ ] **Step 1: Write failing checker lines**

Add block `# --- Task 3`: `need "overwrite = true"` and a count check `[ "$(grep -c 'overwrite = true' "$SKILL")" -ge 4 ] || { echo "FAIL: overwrite = true must appear for timeline, report, trace and dag"; fail=1; }`; `forbid "params.max_"`, `forbid "max_cpus"`, `need "reports/multiqc"`, `forbid "  multiqc/           MultiQC report"`, `need "annotation/"`. Run: expected FAIL.

- [ ] **Step 2: Step 10 config**

Change the `nextflow.config` template: delete the `params { max_cpus ... }` block and set `process { resourceLimits = [ cpus: 16, memory: '64 GB', time: '24h' ] }` with literal values; add `overwrite = true` to `timeline`, `report`, `trace` and `dag`. Add one sentence explaining why (a rerun after a failed launch otherwise leaves an empty trace and no HTML reports because the first launch already created the files; and `params.max_*` are not rnavar parameters and trigger a schema warning). Keep the selector text and the "compare with the trace after the first run" sentence.

- [ ] **Step 3: Step 13 hand-off**

Replace the directory block with the verified 1.3.0 layout: `variant_calling/` (`{SAMPLE}.haplotypecaller.filtered.vcf.gz`, raw `{SAMPLE}.haplotypecaller.vcf.gz`, gVCFs with `generate_gvcf: true`), `preprocessing/` (`{SAMPLE}.md.bam`, and `{SAMPLE}.recal.bam` unless `skip_baserecalibration: true`), `annotation/` (only when `tools` is set), `reports/` (FastQC, samtools, Picard, STAR and `reports/multiqc/{MULTIQC_TITLE}_multiqc_report.html`), `pipeline_info/`. Keep the WebFetch verification sentence.

- [ ] **Step 4: README**

Update `nfcore-rnavar-setup/README.md`: outputs table adds `{date}_{project}_params.yaml`; describe that all parameters go through the params file and why (one sentence); add the custom-reference option to the wizard-steps table; update the Validation status section to: the static checker passes; **one end-to-end run on the nf-core rnavar test data was completed on the cluster** (paired-end, custom reference, known sites, no annotation) — annotation, BAM input, gVCF and the contig-guard rename branch remain untested; prerequisites mention Nextflow >= 24.04 (validated with 26.04.6); remove any statement that the `withName` selectors are unverified and say they were verified against real tasks for STAR_ALIGN, GATK4_SPLITNCIGARREADS, GATK4_BASERECALIBRATOR, GATK4_HAPLOTYPECALLER.

- [ ] **Step 5: Run the checker; expected PASS; commit**

```bash
bash nfcore-rnavar-setup/tests/check_skill.sh nfcore-rnavar-setup/nfcore-rnavar-setup.md
git status && git diff --stat
git add nfcore-rnavar-setup
git commit -m "nfcore-rnavar-setup: config overwrite/resourceLimits, corrected hand-off dirs, README"
```

---

### Task 4: Acceptance on the test dataset (controller-run, not a subagent)

Run by the controller in this session because it needs the user's email and conda env (already provided: `yannvrb@mit.edu`, `nf-env`) and cluster submission.

- [ ] **Step 1:** Install the branch's skill to `~/.claude/commands/nfcore-rnavar-setup.md` (after showing the diff against the installed copy and keeping a `.bak`).
- [ ] **Step 2:** In a fresh directory `/net/bmc-lab3/data/bcc/yannvrb/rnavar_test2/` (copy the test data), follow the updated skill by hand with the same answers as the first run, choosing the "custom reference" option; confirm the generated launch line has only `-params-file` and `{PARAMS_YAML}` has `read_length: 151` unquoted.
- [ ] **Step 3:** Submit the STAR helper then the pipeline with `--dependency=afterok`; expect `Pipeline completed successfully`, a **non-empty** `pipeline_info/execution_trace.txt`, and HTML report/timeline rendered.
- [ ] **Step 4:** Contig-guard rename test: run the generated `prepare_known_sites` helper against a FASTA copy whose contig is renamed `22` (VCFs stay `chr22`), and confirm it either renames and passes or exits 1 with the specified message.
- [ ] **Step 5:** Record results; do not push without the user's approval.

---

## Self-review notes

- Spec coverage: every item of the approved design maps to a task — params file (T1), custom reference, small-genome helper and contig-guard hardening (T2), config overwrite/resourceLimits, hand-off dirs, README (T3), verification incl. the contig-guard rename test (T4). Checker: two new rules (T1). The four deferred Minor items from the last review are covered in T2: SIGPIPE (`tabix -l`), both-empty guard, post-guard filenames, wizard generates the helper on mismatch.
- Consistency: placeholders defined once — `{PARAMS_YAML}`, `{KNOWN_SITES_PARAMS}`, `{VARIANT_PARAMS}`, `{ANNOTATION_PARAMS}` (T1); `{REF_TAG}`, `{GENOME_LENGTH}`, `{SA_INDEX_NBASES}` (T2).
- Deliberately not changed: annotation, BAM-input and gVCF branches beyond their key names; the other deferred Minors (wording, checker hardening of the shared /tmp schema path).
