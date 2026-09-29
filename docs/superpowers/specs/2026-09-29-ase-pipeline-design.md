# Design: `ase-pipeline` skill (allele-specific expression)

Date: 2026-09-29 · Status: draft for review · Sub-project 2 of the RNA-seq variant/splicing roadmap (see `2026-09-28-nfcore-rnavar-setup-design.md`)

## Purpose

An interactive Claude Code skill that generates a complete **allele-specific expression (ASE)** analysis for bulk RNA-seq on the SLURM + Singularity cluster: reference preparation, allele-aware alignment, allele counting, and downstream statistics as R Markdown reports. It supports two experimental designs behind one wizard:

- **F1 cross mode** (for example B6 × AJ mice): every SNP that differs between the parental strains is heterozygous by construction; mapping bias is handled with a **third-allele masked reference**.
- **Outbred / human mode**: heterozygous SNPs come from a per-individual genotype VCF; mapping bias is handled with **STAR's WASP re-mapping filter**.

Success criteria:
- A user with FASTQs, a reference and a genotype (or parental-difference) VCF gets runnable sbatch scripts and Rmds that finish on the cluster and produce per-SNP and per-gene imbalance tables, without hand edits.
- The report never silently presents biased ratios: a reference-bias diagnostic is always produced and flagged when the mean reference-allele fraction is far from 0.5.
- On synthetic ground-truth data, planted imbalance is recovered at FDR 0.05, and null genes show near-nominal false positives.

Non-goals: RNA editing detection (needs matched DNA to exclude SNPs; parked in the roadmap), variant calling itself (that is `nfcore-rnavar-setup`), and any nf-core wrapper (no nf-core pipeline exists for ASE).

## Decisions already made (with the user, 2026-09-29)

| Decision | Choice |
|---|---|
| Skill shape | Standalone `ase-pipeline` skill; not an rnavar option and not a `bulk-rnaseq-pipeline` module (that skill starts from Salmon `quant.sf`, ASE needs BAMs and genotypes) |
| Modes | Both F1 cross and outbred/human, one wizard with a mode fork (like paired-end vs 3′ DGE in `bulk-rnaseq-pipeline`) |
| Genotype source (outbred) | Matched WGS/array VCF when available, or the rnavar VCF; the wizard says plainly that RNA-derived genotypes are biased (see Risks) |
| F1 SNP source | A user-supplied parental-difference VCF is the core input (any organism); a mouse helper extracts it from the Mouse Genomes Project multi-strain VCF |
| Analyses in v1 | All four: per-sample imbalance, reciprocal F1, differential ASE, read-backed phasing (phASER) — as one spec, staged plan |
| Upstream engine | Generated SLURM scripts and array jobs (no Nextflow/Snakemake); tools from Singularity biocontainers (STAR, GATK, samtools, bcftools) |

## Architecture

The skill is one Markdown wizard file (`ase-pipeline/ase-pipeline.md`), written in the style of `bulk-rnaseq-pipeline` (numbered steps, numbered-option questions, mode fork, "Notes for the assistant"). It generates:

1. Upstream sbatch scripts: one reference-prep job, then a per-sample array job. STAR runs from the cached 2.7.10b container for both index building and alignment.
2. Downstream Rmd files (self-contained, no `source()` of helper files; `cache = FALSE`; Bioconductor packages loaded before tidyverse), each with its own sbatch script.
3. A summary report page linking all outputs (same pattern as `bulk-rnaseq-pipeline` Step 15).

Repository layout: `ase-pipeline/ase-pipeline.md`, `ase-pipeline/README.md`, `ase-pipeline/tests/check_skill.sh`, `ase-pipeline/tests/fixtures/` (recorded `--help` flag lists for STAR and `gatk ASEReadCounter`), plus a row in the root README and an installed copy in `~/.claude/commands/`.

## Wizard steps

0. Working directory, `{WD_NAME}`, `{TODAY_YYMMDD}`, results directory `results/{YYMMDD}_{WD_NAME}` (same convention as `bulk-rnaseq-pipeline`; raw FASTQ/BAM/VCF inputs are read-only).
1. Email. 2. A conda environment is **not** needed (no Nextflow); the phASER stage creates its own environment.
3. Mode: **F1 cross** or **outbred/human** (numbered, with the difference explained).
4. Organism and reference: FASTA and GTF (custom paths or the standard folder convention used by the other skills), read length detection (mode of read lengths; `sjdbOverhang = read length − 1`, a separate index per read length, same rule as `nfcore-rnavar-setup`).
5. Sample sheet `{WD_NAME}_samples.csv` (CSV, not xlsx, because no R runs on the login node), scaffolded from the FASTQs. Columns: `sample, fastq_1, fastq_2, condition`, plus `individual` (outbred) or `cross_direction` and the two strain names (F1). Validation: no dashes/spaces in names; at least one replicate per condition warned.
6. Genotype source: (F1) the parental-difference VCF or the mouse helper; (outbred) per-individual VCF path(s) from external WGS/array data, or the rnavar `variant_calling/` VCFs.
7. Analysis menu (numbered, multi-select) with enforced constraints: per-sample imbalance always on; reciprocal F1 only in F1 mode and only when both cross directions are present; differential ASE only with at least two conditions with replicates; phASER only in outbred mode.
8. Site-filter and significance constants, defaults shown and editable: `MIN_DEPTH` 10 total reads per SNP, `FDR_SIG` 0.05, `ABS_DEV_SIG` 0.1 (|alt fraction − 0.5|).
9. SLURM resources (defaults at most 64 G and 4 h; STAR alignment for human is the largest request).
10. Project name and author for the Rmd headers.
11–13. Generate the upstream scripts, the Rmds and their sbatch scripts, then the summary report.

## Upstream — F1 cross mode

**Reference prep job** (`prep_f1_reference.sh`):
1. Subset the parental-difference VCF to **biallelic SNPs** where the two strains genotype differently (indels and multi-allelic sites are excluded from counting and documented as a limitation).
2. Build the **third-allele masked FASTA**: at each SNP position write a base that is neither strain's allele (deterministic rule: the first of A, C, G, T that is not one of the two alleles). Record the masked positions in a BED. Both parental alleles then mismatch the reference equally, so neither strain is favoured in mapping.
3. Build a STAR index (`sjdbOverhang = read length − 1`, `--genomeSAindexNbases` from the genome length, formula as in `nfcore-rnavar-setup`).
4. The reference-prep job also writes `f1_het_sites.vcf.gz` (single sample `F1`, genotype `0/1` at every parental SNP, bgzipped and tabix-indexed) because ASEReadCounter counts at heterozygous sites; the parental-difference VCF stays the source of truth. Task 2 confirmed that a sites-only VCF gives 0 rows, so this file is required. The reference also needs a `.fai` (`samtools faidx`) and a `.dict` (`gatk CreateSequenceDictionary`).

**Per-sample array job** (`align_count_f1.sh`):
1. STAR alignment to the masked reference with `--outSAMattrRGline ID:{sample} SM:{sample} PL:ILLUMINA` (without a read group ASEReadCounter's default read-group filter silently drops every read), coordinate-sorted BAM, then samtools index and duplicate marking.
2. **GATK ASEReadCounter** (GATK 4.4.0.0 container) at the parental SNP sites, using `f1_het_sites.vcf.gz` (bgzipped, tabix-indexed, heterozygous genotype column: required, verified) with the original REF/ALT alleles from the parental VCF, independent of the masked FASTA base: required arguments `--input` and `--variant`; the script also passes `--reference`, `--output`, `--min-base-quality` (default 0), `--min-mapping-quality` and `--min-depth-of-non-filtered-base` from the wizard's site-filter constants, and leaves `--count-overlap-reads-handling` at its default `COUNT_FRAGMENTS_REQUIRE_SAME_BASE`. ASEReadCounter's default read filters include `NotDuplicateReadFilter`, so duplicates are **marked, not removed**, and are skipped by the tool itself. Output: the standard ASEReadCounter table.
3. Alleles are mapped to strains through the VCF (REF strain versus ALT strain), so every table row has strain-A and strain-B counts.

**Mouse helper** (`extract_mgp_parental_vcf.sh`): extract the two named strains from the Mouse Genomes Project multi-strain VCF (about 22 GB, `.csi` index) with `bcftools view -s` and keep sites where the genotypes differ. Because the source file is very large, the helper streams by region/chromosome and never stores the full VCF twice. It runs as a per-chromosome array job followed by `bcftools concat` (about 80 s fixed remote overhead and about 80 s per 10 Mb, so a serial whole-genome run would take about 6 hours). MGP contigs are `1`, `2`, ... (no `chr` prefix), so the FASTA must use the same names. The release (`REL-2112-v8-SNPs_Indels`, files `mgp_REL2021_snps.vcf.gz` and `mgp_REL2021_indels.vcf.gz`) is on **GRCm39** with C57BL/6J as the reference, contains 52 strains whose VCF sample names use underscores (`A_J`, `CAST_EiJ`, `C57BL_6NJ`, `WSB_EiJ`, ...), encodes the reference allele as `0/0`, marks sites `PASS` or `LowQual`, and records per-genotype confidence in the `FI` tag (1 = high confidence): the helper keeps only genotypes with `FI=1` for both strains (verified from the release README, 2026-09-29).

## Upstream — outbred / human mode

**Genotype prep** (`prep_genotypes.sh`): for each individual produce a **single-sample** VCF (STAR reads genotypes from the 10th column, so multi-sample files are not used) restricted to biallelic heterozygous SNPs: `bcftools view -s SAMPLE -g het -v snps -m2 -M2`, plus quality filters (external WGS/array: `GQ` and `DP` thresholds; rnavar-derived: `FILTER=PASS`, `DP`, `QUAL`). The wizard states the RNA-derived genotype caveat.

**Per-sample array job** (`align_wasp_count.sh`): STAR with `--varVCFfile <het VCF> --waspOutputMode SAMtag` (the VCF must be the heterozygous-only one: STAR does not ignore homozygous genotypes, verified) and `--outSAMattributes` including `NH HI AS nM vA vG vW`; keep alignments with `vW:i:1` and alignments carrying no `vW` tag (reads not overlapping any SNP); index and mark duplicates; run **ASEReadCounter** at the individual's heterozygous sites. Reads with `vW` values 2–7 failed WASP filtering and are excluded. The script also writes a WASP removal table per sample (alignments by `vW` code and by `vA` allele class), shown in Rmd 01.

## Downstream Rmd modules

All Rmds are self-contained, use the shared constants (`MIN_DEPTH`, `FDR_SIG`, `ABS_DEV_SIG`) for tables and figures so counts can never disagree, and write xlsx tables plus a checkpoint RDS.

**Rmd 01 — Import and QC.** Read the ASEReadCounter tables; apply site filters (total depth at least `MIN_DEPTH`; exclude sites with a high fraction of low-quality or other-base depth or improper pairs); map alleles to strains (F1) or REF/ALT (outbred); **reference-bias diagnostic**: distribution and mean of the reference-allele fraction across heterozygous sites per sample, flagged when the deviation exceeds max(0.03, 3 × SE), SE = sqrt(0.25 / total reads); per-sample coverage and site counts.

**Rmd 02 — Per-sample allelic imbalance (always).** Per-SNP two-sided beta-binomial test against p = 0.5, overdispersion (ρ) estimated per sample by maximum likelihood (implemented with base R `optim`, so this stage needs no extra packages), BH adjustment per sample, significance = `FDR_SIG` and `ABS_DEV_SIG`. Gene-level results: **F1 mode** — SNPs are phased by strain, so gene counts are exact sums of the two strains' reads across the gene's SNPs, tested with the same beta-binomial model; **outbred without phasing** — SNPs cannot be pooled, so the report gives per-SNP tests and a gene-level combined p-value (Cauchy combination of the SNP p-values), explicitly labelled "unphased, no direction"; **outbred with phASER** — the gene-level haplotypic counts replace pooled SNP counts.

**Rmd 03 — Reciprocal F1 (F1 mode, both cross directions present).** Per gene, a beta-binomial model `logit(p) = β0 + β1·d` where p is the strain-A allele fraction and d = +1 when strain A is the maternal allele and −1 when it is the paternal allele. **β0 is the strain effect (cis-regulatory divergence); β1 is the parent-of-origin effect (imprinting).** Requires replicates in each direction.

**Rmd 04 — Differential ASE (at least two conditions with replicates).** `logit(p) = β0 + β1·condition`, with an individual term for outbred samples; per-SNP tests pair the same SNP across conditions within an individual; F1 mode also tests per gene on strain-summed counts.

The beta-binomial GLMs in Rmd 03 and 04 use `aod::betabin` (fixed effects with a dispersion parameter), and `lme4::glmer` for the individual random effect in Rmd 04; both are **already in the `bulkrnaseq` image** (aod 1.3.3, lme4 2.0.6, verified 2026-09-29), so no image rebuild is needed. `VGAM`, `glmmTMB`, `MBASED` and `VariantAnnotation` are absent and are not used. The wizard's package check (Step 2 style of `bulk-rnaseq-pipeline`) verifies the six mandatory packages (`aod`, `lme4`, `openxlsx`, `tidyverse`, `GenomicRanges`, `rtracklayer`) and stops with a rebuild request if any is missing.

**Stage 3 — phASER (outbred only, gated).** A helper creates a conda environment from phASER's current `environment.yml` (Python 3.9 or later, conda-forge and bioconda packages), runs phASER per sample on the alignment BAM and the genotype VCF, then runs the gene-level haplotypic expression step with a gene-annotation BED. The stage is kept only if the environment installs and a test run succeeds during Stage 3; otherwise it is dropped and documented.

## Task 0 — verification gate (before any Stage 1 implementation)

**Already verified (2026-09-29), fixtures to be copied into `ase-pipeline/tests/fixtures/`** from `/net/bmc-lab3/data/bcc/yannvrb/rnavar_test2/ase_fixtures/`:
- `gatk ASEReadCounter --help` from the cached GATK 4.4.0.0 container (`gatk_ASEReadCounter_help.txt`): required `--input`/`--variant`; optional `--reference`, `--output`, `--min-base-quality`, `--min-mapping-quality`, `--min-depth-of-non-filtered-base`, `--count-overlap-reads-handling`, `--output-format`, `--read-filter`, `--disable-read-filter`, `--max-depth-per-sample`.
- STAR 2.7.10b (`star_version.txt`, `star_help.txt`): `--varVCFfile`, `--waspOutputMode`, and the `vA`/`vG`/`vW` attributes exist (`vW`: 1 = passes WASP filtering; 2–7 = fails).
- R packages in the `bulkrnaseq` image (`r_packages.txt`): `aod`, `lme4`, `openxlsx`, `tidyverse`, `apeglm`, `GenomicRanges`, `rtracklayer` present; `VGAM`, `glmmTMB`, `MBASED`, `VariantAnnotation` absent.
- Mouse Genomes Project release: GRCm39, 52 strains, sample-name format, `0/0` reference encoding, `PASS`/`LowQual`, `FI` genotype-confidence tag (see the mouse helper above).

**Items 1-3 (STAR WASP behaviour, `.csi` region query, ASEReadCounter columns) were completed on 2026-09-29, with results in `ase-pipeline/tests/fixtures/verification.md`. Still open before Stage 1 implementation of the masking step; the plan does not proceed on assumptions:**
1. DONE: STAR behaviour with WASP, from one small test run in the cached 2.7.10b container: whether homozygous or unphased genotypes in column 10 are ignored, whether only the first sample of a multi-sample VCF is used, which `--outSAMattributes` and `--outSAMtype` settings are required for `vW`, and the interaction with multi-mappers.
2. DONE: Whether the Mouse Genomes Project `.csi` index supports region queries with `bcftools view -r` from the cluster, and how long extracting two strains for one chromosome takes (sizes the helper's sbatch request).
3. DONE: The full ASEReadCounter output-table column names, from a test run on the small synthetic data (used by Rmd 01).
4. The third-allele masking script: the user is recovering scripts from the earlier published analysis; if found they inform the implementation, otherwise the deterministic rule above is used.

## Testing and acceptance

Static checker (`ase-pipeline/tests/check_skill.sh`), in the same style as `nfcore-rnavar-setup/tests/check_skill.sh`: required and forbidden strings, and every flag of STAR and `gatk ASEReadCounter` used in the skill's script templates must appear in the recorded `--help` fixtures. It cannot catch runtime typing problems, so each stage also has a cluster acceptance run.

**Synthetic ground-truth datasets (both modes):** paired-end reads simulated from the two haplotypes on a small region with planted imbalance in chosen genes (for example 70:30) and 50:50 elsewhere; the outbred version also plants reference-biased mapping that the WASP filter must remove. **Real-data smoke tests:** a small public human RNA-seq restricted to chr22 with matching genotypes, and a public F1 mouse dataset if one is small enough; the user's own B6 × AJ F1 data is the preferred final check.

Acceptance criteria (measurable):
- Planted-imbalance genes are recovered at FDR 0.05; false positives among null genes are close to the nominal rate.
- On the null genes only, the mean reference-allele fraction is within max(0.03, 3 × SE) of 0.5 both without and with WASP (measured +0.02 unfiltered and −0.03 filtered on the synthetic data; no improvement promised); the masked-versus-unmasked comparison is reported honestly (an earlier analysis found little difference from B6-only mapping, so no improvement is promised).
- Under a simulated null the beta-binomial p-values are approximately uniform (calibration check run in the container).
- Every Rmd builds from the real checkpoint files without hand edits.

## Staging

- **Stage 1:** wizard, F1 and outbred upstream, Rmd 01 and 02, checker, README, synthetic acceptance run.
- **Stage 2:** Rmd 03 (reciprocal F1) and Rmd 04 (differential ASE), with a synthetic reciprocal design.
- **Stage 3:** phASER, gated as above.

Each stage follows the process used for `nfcore-rnavar-setup`: plan, subagent-driven implementation with per-task review, final whole-branch review, a live cluster acceptance run, and the user's approval before any merge or push.

## Risks and known limitations

- **RNA-derived genotypes are circular:** heterozygous sites with strong imbalance may be called homozygous, so genotypes called from the same RNA-seq bias results toward balance and miss lowly expressed sites. External genotypes are preferred; the wizard states this whenever the rnavar VCF is chosen.
- **Indels and multi-allelic sites** are excluded from counting; masking covers only biallelic SNPs.
- **Multi-mapping and paralogous regions** can give spurious ratios; low-mappability filtering is optional and documented, not automatic.
- **Very large inputs:** the Mouse Genomes Project VCF is about 22 GB, so the helper must stream by region and needs a long enough sbatch limit within the 4 h default or an explicit user choice.
- **phASER** is an older tool whose maintenance status is unclear; the stage is gated.
- Downstream Rmds 03 and 04 require `aod` and `lme4` in the container (present today); the wizard verifies them and stops with an image-rebuild request if either disappears, rather than silently substituting a different model.
