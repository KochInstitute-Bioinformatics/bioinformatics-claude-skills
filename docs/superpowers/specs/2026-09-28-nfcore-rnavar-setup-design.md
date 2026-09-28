# Design: `nfcore-rnavar-setup` skill

Date: 2026-09-28 · Status: draft for review · Sub-project 1 of 3 (see Roadmap)

## Purpose

An interactive Claude Code skill that sets up and submits **nf-core/rnavar** (GATK best-practice RNA-seq short-variant calling: 2-pass STAR → MarkDuplicates → SplitNCigarReads → BQSR → HaplotypeCaller → soft filtering → optional annotation) on the SLURM + Singularity cluster. It follows the wrap-an-nf-core-pipeline pattern already validated by `nfcore-rnaseq-setup` and `nfcore-scrnaseq-setup`, and keeps the same step order, question style and conventions.

Success criteria:
- A user with FASTQs and a genome build gets a correct `nextflow.config`, samplesheet CSV and `nf-core_rnavar_{version}.sh` that submit cleanly with `sbatch`.
- The skill never lets the pipeline fail late for a foreseeable reason (missing known sites, wrong `--read_length`, annotation caches that cannot be downloaded from compute nodes).
- The run's outputs (filtered VCFs, recalibrated BAMs) are reported in a hand-off note usable by sub-project 2 (`ase-pipeline`).

Non-goals: annotation-driven biology (kept optional), RNA editing detection, allele-specific expression, and any downstream R analysis.

## Roadmap (context)

| # | Sub-project | Depends on | Status |
|---|---|---|---|
| 1 | `nfcore-rnavar-setup` (this spec) | — | design |
| 2 | `ase-pipeline` (F1 N-masked mode + outbred WASP mode; genotypes from external VCF or rnavar VCF; beta-binomial testing) | 1 (optional for F1 mode) | separate spec |
| 3 | `nfcore-rnasplice-setup` (rMATS, SUPPA2, LeafCutter, DEXSeq/DRIMSeq, IsoformSwitchAnalyzeR; overlaps `bulk-rnaseq-pipeline` Rmd 04 and must be positioned against it) | — | separate spec |
| — | RNA editing detection (needs matched DNA/parental genotypes to exclude SNPs) | 1 | parked |

## Wizard steps

Unchanged from `nfcore-rnaseq-setup`: Step 0 working directory / `{WD_NAME}` / `{TODAY_YYMMDD}`; Step 1 email; Step 2 conda env (asked separately from email); Step 3 latest version via `WebFetch` of `https://api.github.com/repos/nf-core/rnavar/releases/latest` (never `gh`, which is not installed); numbered-option questions; MultiQC title; output directory; sample-name sanitisation and collision warning.

Specific to rnavar:

**Step 4 — FASTQ scan and samplesheet.** Same scan and paired-end detection. Columns: `sample,fastq_1,fastq_2` (no `strandedness`). Also accept pre-aligned input (`bam`/`bai` or `cram`/`crai`), one file type per sample — mixing FASTQ and BAM for the same sample errors in the pipeline, so the skill checks for it. Sample names shared across rows are merged as lanes (`CAT_FASTQ`), so the duplicate-name warning from the rnaseq skill applies unchanged. Single-end is allowed by the pipeline but flagged as unusual for variant calling.

**Step 5 — Read length.** Detect from the FASTQs with the same `awk` one-liner as the rnaseq skill and always pass `--read_length {N}` (schema default 150 would silently give a wrong `sjdbOverhang` for other lengths). Report "Detected read length N → sjdbOverhang N−1".

**Step 6 — Organism and genome.** Same folder convention `{genome_base}/{organism}/{assembly}_ens{version}/`, reusing an existing FASTA/GTF. The STAR index is **built separately for rnavar** (`index/star_rnavar_sjdb{N-1}/`) because `sjdbOverhang` is fixed at build time and the rnaseq skill's index may have a different one. Passing `--star_index` is optional; if omitted rnavar builds it, but the skill prefers a pre-built index (with `--save_reference` off) to avoid rebuilding per run. The GTF source (Ensembl vs GENCODE) is still detected and reported, but unlike rnaseq **no flag is emitted**: the 1.3.0 schema has no `gencode` parameter, so emitting it would fail parameter validation.

**Step 7 — Known sites (BQSR).** Required unless BQSR is skipped: the pipeline passes `dbsnp`/`known_indels` straight to BaseRecalibrator with no automatic fallback, so missing files make the run fail late. The skill therefore asks explicitly:
1. Human: use the GATK resource-bundle dbSNP and known-indels VCFs (with `.tbi`) — pass `--dbsnp`, `--dbsnp_tbi`, `--known_indels`, `--known_indels_tbi`.
2. Mouse: use Mouse Genomes Project variants, downloaded and tabix-indexed by a generated helper script.
3. Skip base recalibration: add `--skip_baserecalibration` and state the trade-off (slightly less accurate base qualities; the appropriate choice for organisms without curated variant sets).
The skill must resolve resource URLs by fetching the current Ensembl/GATK locations at run time and must not embed URLs it has not verified.

**Step 8 — Variant options.** Numbered choices, default-valued flags omitted from the script:
- Duplicates: keep marked (default) or `--remove_duplicates`.
- `--star_twopass` is true by default; only emit `--star_twopass false` if the user disables it.
- HaplotypeCaller confidence `--gatk_hc_call_conf` (default 20) and soft filters `--gatk_vf_qd_filter` (2), `--gatk_vf_fs_filter` (30), `--gatk_vf_window_size` (35), `--gatk_vf_cluster_size` (3): keep defaults unless the user asks; `--skip_variantfiltration` optional.
- `--generate_gvcf` when the user plans joint calling across samples.
- `--bam_csi_index` only for genomes with chromosomes > 512 Mb (and it disables filtration — state that).

**Step 9 — Optional annotation.** Ask: none / SnpEff / VEP / merge, via `--tools`. If chosen, annotation caches must exist on disk before submission (`--snpeff_cache`, `--vep_cache`, plus `--snpeff_db`, `--vep_genome`, `--vep_species`, `--vep_cache_version`). `--download_cache` needs internet from compute nodes and can stall; the skill states this and offers a login-node-free pre-download job instead. Note: `--annotation_cache` mentioned on the rnavar usage page is **not** in the 1.3.0 schema and must not be emitted.

**Step 10 — `nextflow.config`.** Generated only if absent (never overwritten). Same SLURM/Singularity template as the rnaseq skill with `resourceLimits` capped at 64 GB / 24 h. rnavar's `base.config` defines label-based resources only (`process_medium` 6 CPU/36 GB/8 h, `process_high` 12 CPU/72 GB/16 h), so the config sets `withName` overrides using regex selectors (`'.*:FASTQ_ALIGN_STAR:.*'`, `'.*:GATK4_HAPLOTYPECALLER'`, `'.*:GATK4_BASERECALIBRATOR'`) rather than a hard-coded workflow prefix, with STAR at ≤ 64 GB. Exact selector strings are confirmed against a real run's trace file during implementation testing.

**Step 11 — Submission script.** `nf-core_rnavar_{version}.sh`, same header and conda/Singularity block as the rnaseq skill:

```
nextflow run nf-core/rnavar -r {VERSION} -c nextflow.config -profile slurm,singularity \
  --input {SAMPLESHEET_CSV} --fasta {FASTA} --gtf {GTF} --star_index {STAR_INDEX} \
  --read_length {N} --seq_platform illumina {KNOWN_SITES_LINES} {VARIANT_LINES} \
  {ANNOTATION_LINES} --multiqc_title {TITLE} --outdir {OUTDIR}
```
Flags equal to pipeline defaults are omitted.

**Step 12 — Helper scripts (only when resources are missing):** genome + STAR index build (rnaseq skill's Step 14b logic with `sjdbOverhang = N−1`), known-sites download/index, annotation-cache pre-download.

**Step 13 — Hand-off note.** After writing scripts, print where the outputs will land (filtered VCFs, recalibrated BAMs under `{OUTDIR}`) and state that they are the inputs expected by `ase-pipeline`.

## Verified facts (nf-core/rnavar 1.3.0, 2026-06-03)

Checked against `nextflow_schema.json`, `conf/base.config` and `workflows/rnavar.nf` at tag 1.3.0:
- Required parameters: `input`, `outdir`, `aligner` (default `star`), `seq_platform` (default `illumina`).
- No automatic skipping of BQSR when known sites are missing (workflow errors).
- Duplicate sample names are merged as lanes via `CAT_FASTQ`.
- Annotation is selected with `--tools` (`snpeff`, `vep`, `bcfann`, `merge`) and `--skip_tools`.
- Module names present in the workflow include `FASTQ_ALIGN_STAR`, `BAM_MARKDUPLICATES_PICARD`, `GATK4_BASERECALIBRATOR`, `GATK4_HAPLOTYPECALLER`, `GATK4_VARIANTFILTRATION`; subworkflows `SPLITNCIGAR`, `RECALIBRATE`, `VCF_ANNOTATE_ALL`.

## Error handling and notes for the assistant

- Never run heavy computation on the login node; all work through `sbatch`.
- Raw FASTQ/BAM files are read-only.
- Present every finite-choice question as a numbered list.
- Use `WebFetch` for all GitHub API calls.
- For non-mouse/non-human organisms, ask the user for FASTA/GTF/known-sites paths and offer `--skip_baserecalibration`.

## Testing and acceptance

The skill is a prompt, so acceptance is behavioural:
1. Run the skill against a small real dataset (paired-end human or mouse) and confirm the generated samplesheet, config and script pass `nextflow run ... -preview` (or an equivalent parse check) on a compute node.
2. Run the pipeline's own `-profile test,singularity` on a compute node to confirm the environment, Singularity image pulls, and config resource selectors match real process names (read them from the trace file).
3. Confirm the two failure paths behave as designed: missing known sites without `--skip_baserecalibration` is prevented by the wizard, and `--read_length` always matches the detected value.
4. Confirm the hand-off note lists real, existing output paths after a successful run.
