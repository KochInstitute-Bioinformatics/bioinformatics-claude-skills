# ase-pipeline — allele-specific expression (F1 cross and outbred/human)

You are helping the user set up a complete allele-specific expression (ASE) analysis for bulk RNA-seq on an HPC cluster using SLURM and Singularity: reference preparation, allele-aware alignment (STAR), allele counting (GATK ASEReadCounter), and downstream statistics as R Markdown reports. Walk through each step below in order, asking the user only what you need and performing automated steps silently.

Two experimental designs share one wizard:

- **F1 cross** (for example B6 × AJ mice): every SNP that differs between the two parental strains is heterozygous by construction. Mapping bias is handled with a third-allele masked reference.
- **Outbred / human**: heterozygous SNPs come from a per-individual genotype VCF. Mapping bias is handled with STAR's WASP re-mapping filter.

**Honesty rule.** Allelic ratios are never reported without the reference-bias diagnostic (the distribution and mean of the reference-allele fraction across heterozygous sites, per sample). If the diagnostic is flagged, say so next to every ratio table and figure. No improvement from masking or WASP is promised; the report shows the comparison as measured.

**HPC rule.** Never run R, STAR, GATK or any heavy job on the login node. Everything that is not a trivial shell command goes through `sbatch` (debug with `salloc` or `srun`). Default requests are at most 64 G and 4 h unless the user asks otherwise. Raw FASTQ, BAM and VCF inputs are read-only.

---

## Step 0 — Establish working directory

Before anything else, run `pwd` to record the current working directory (`{CWD}`). All scripts, the sample sheet, the Rmd files and the results are written here. Inform the user: "Working directory: `{CWD}`"

Also derive:
- `{WD_NAME}` = basename of `{CWD}`
- `{TODAY_YYMMDD}` = today's date formatted as `YYMMDD`
- `{RESULTS_DIR}` = `{CWD}/results/{TODAY_YYMMDD}_{WD_NAME}` (create it with `mkdir -p` when the first output is written)

---

## Step 1 — Email address

Ask the user: "What email address should SLURM use for job notifications (completion/failure)?" Store it as `{USER_EMAIL}`.

Do **not** pre-suggest or pre-fill any email address. Wait for the answer before Step 2.

---

## Step 2 — Environment check (no R on the login node)

No conda environment is needed for this skill. Check the three things below; run them silently and report a one-line result for each.

**(a) Singularity.** All tools except R run from cached Singularity biocontainers. Load the module with a checked `module add` and confirm the command exists:

```bash
module add singularity/3.10.4 || { echo "ERROR: cannot load singularity" >&2; exit 1; }
command -v singularity || { echo "ERROR: singularity not on PATH" >&2; exit 1; }
```

**(b) Cached containers.** Set `{SIF_DIR}` = `${NXF_SINGULARITY_CACHEDIR:-$HOME/.singularity/cache}`. Container files are named `depot.galaxyproject.org-singularity-<name>-<tag>.img`. Locate each one below and store the full path in the variable shown:

| Variable | File name in `{SIF_DIR}` | Download URL |
|---|---|---|
| `{STAR_SIF}` | `depot.galaxyproject.org-singularity-star-2.7.10b--h9ee0642_0.img` | `https://depot.galaxyproject.org/singularity/star:2.7.10b--h9ee0642_0` |
| `{GATK_SIF}` | `depot.galaxyproject.org-singularity-gatk4-4.4.0.0--py36hdfd78af_0.img` | `https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0` |
| `{BCFTOOLS_SIF}` | `depot.galaxyproject.org-singularity-bcftools-1.20--h8b25389_0.img` | `https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0` |
| `{SAMTOOLS_SIF}` | `depot.galaxyproject.org-singularity-samtools-1.21--h50ea8bc_0.img` | `https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0` |
| `{PICARD_SIF}` | `depot.galaxyproject.org-singularity-picard-3.1.1--hdfd78af_0.img` | `https://depot.galaxyproject.org/singularity/picard:3.1.1--hdfd78af_0` |

All five URLs were verified (HTTP 200) on 2026-09-29 (`ase-pipeline/tests/fixtures/verified_urls.txt`). For any file that is missing, the generated helper script downloads it (these images are large, so never download on the login node in the foreground):

```bash
[ -s "$SIF" ] || { wget -c -O "$SIF.part" "URL" && mv "$SIF.part" "$SIF"; } || { echo "ERROR: could not download $SIF" >&2; exit 1; }
```

Tell the user which containers are already cached and which the helper will download.

**(c) R packages.** R is only available in the `bulkrnaseq` image, `/net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif` (`{R_SIF}`). Do not start R on the login node. Write `check_r_packages.sh` and submit it as a tiny `sbatch` job (`-n 1 --mem=4G -t 0:10:00`):

```bash
#!/bin/bash
#SBATCH -J ase_rpkgs
#SBATCH -p bcc
#SBATCH -n 1 --mem=4G -t 0:10:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={USER_EMAIL}
#SBATCH -o check_r_packages_%j.out
module add singularity/3.10.4 || { echo "ERROR: cannot load singularity" >&2; exit 1; }
singularity exec --bind {CWD} {R_SIF} \
  Rscript -e 'for (p in c("aod","lme4","openxlsx","tidyverse","GenomicRanges","rtracklayer")) cat(p, requireNamespace(p, quietly = TRUE), "\n")'
```

Mandatory packages: `aod`, `lme4`, `openxlsx`, `tidyverse`, `GenomicRanges`, `rtracklayer`. Submit it with `sbatch -p bcc check_r_packages.sh` (every job is submitted with `-p bcc`), wait for the job to finish, then read the result from the job output file `check_r_packages_<jobid>.out`. If any mandatory package prints `FALSE`, **stop** and ask the user to rebuild the `bulkrnaseq` image; do not substitute another model or package. (`VGAM`, `glmmTMB`, `MBASED` and `VariantAnnotation` are absent from the image and must never be required.)

---

## Step 3 — Mode

Ask the user (numbered options):

"Which experimental design is this?
1. F1 cross
2. Outbred / human"

- **1. F1 cross** — the animals are F1 offspring of two inbred parental strains. Every SNP that differs between the strains is heterozygous by construction, so genotypes come from a parental-difference VCF (or the mouse helper), not from the animals. Reads are mapped to a third-allele masked reference so neither strain is favoured in mapping.
- **2. Outbred / human** — each individual has its own genotype. Heterozygous SNPs come from a per-individual VCF (external WGS/array data, or the rnavar VCF with a stated caveat). Reads are mapped to the normal reference and STAR's WASP filter removes reads whose mapping depends on the allele they carry.

Store the answer as `{MODE}` = `f1` or `outbred`. Every later step that differs between the modes branches on `{MODE}`.

---

## Step 4 — Reference

Ask (numbered): "1. Ensembl-style reference in the standard folder (default) · 2. Custom reference: I already have a FASTA and GTF". Never download or decompress anything on the login node in the foreground; the prep job does it.

- **Option 1 (standard folder):** ask for the organism and the base directory that holds genomes (`{genome_base}`). Use the folder convention `{genome_base}/{organism}/{assembly}_ens{version}/` (for example `{genome_base}/mouse/mm39_ens{version}/` or `{genome_base}/human/hg38_ens{version}/`). If `{genome_base}/{organism}/` already has `{assembly}_ens{N}` folders, use the highest N unless the user asks otherwise. Set `{GENOME_DIR}` to that folder and `{FASTA_PATH}` / `{GTF_PATH}` to the primary-assembly FASTA and the GTF inside it. If they are missing, the prep job downloads them (URLs resolved at run time from the current Ensembl release listing and shown to the user first; never typed from memory).
- **Option 2 (custom):** ask for the FASTA path, the GTF path and the directory that will hold the indexes (`{GENOME_DIR}`). If a file ends in `.gz`, the prep job decompresses it with `gunzip -c file.gz > {GENOME_DIR}/<name>` and `{FASTA_PATH}` / `{GTF_PATH}` point at the decompressed copies; the originals are never modified.

The reference also needs a `.fai` and a sequence `.dict`; the prep job creates them (do not create them here). In F1 mode with the mouse helper the FASTA must be GRCm39 with Ensembl-style contig names (`1`, `2`, ..., no `chr` prefix), matching the Mouse Genomes Project VCF.

The read length, `{SJDB_OVERHANG}` and `{STAR_INDEX}` depend on the FASTQs, so they are determined at the end of Step 5, once the FASTQ directory is chosen.


---

## Step 5 — Sample sheet

Scan `{CWD}` for FASTQs (`find {CWD} -name "*.fastq.gz" -o -name "*.fq.gz" | head -30`), list the unique directories and ask which to use (or for another path). Use paths relative to `{CWD}` when the input is inside it, otherwise absolute.

**Paired-end detection:** `_R1_`/`_R2_`, `_1.fastq.gz`/`_2.fastq.gz`, or `_1_sequence`/`_2_sequence` mean paired-end. Otherwise single-end: warn "Single-end data works but gives fewer usable allele-specific reads; continue?". Pair each R1 with its R2 by substitution; warn and leave `fastq_2` empty if R2 is missing.

**Name sanitisation (always, before showing the user):** sample name = filename up to `_S\d+`, `_R1` or `_1_sequence`; replace `-`, spaces, `/`, `(`, `)` and other special characters with `_` and note substitutions in the preview. After sanitisation, names must be unique; on a collision warn: "Name collision '{NAME}': each row is aligned as its own sample, so duplicate names would overwrite each other's files. Provide distinct names." and go to custom naming.

Scaffold `{SAMPLES_CSV}` = `{CWD}/{WD_NAME}_samples.csv` (CSV, not xlsx) with the columns:

```
sample,fastq_1,fastq_2,condition,cross_direction,individual
```

- `condition`: ask the user (plain-language description of the groups, then map samples to it).
- **F1:** ask once for the two strain names, `{STRAIN_A}` (the reference strain) and `{STRAIN_B}`, and, per sample, `cross_direction` (for example `AxB` or `BxA`, maternal strain first). Leave `individual` equal to `sample` (each F1 animal is its own individual).
- **Outbred:** `individual` per row (which person/animal the sample comes from; several samples may share one); leave `cross_direction` as `NA`.

Show the full table, ask "Does this look correct?", and write the file only after confirmation (if it exists, ask: overwrite or choose another filename).

**Validation rules:** no dashes or spaces in `sample`, `condition`, `cross_direction` or `individual`; sample names unique; each condition should have replicates (defined below; warn, do not stop, when it has only one sample); reciprocal analysis needs both cross directions present (checked again in Step 7).

**Replicates (definition used here and in Step 7):** replicates = at least 2 samples in the condition. A condition with a single sample triggers the warning above.

### Step 5, part 2 — Read length and STAR index

**Read length (always detect it).** STAR's `sjdbOverhang` is fixed when the index is built. Define `{FASTQ_FILE}` as the first `fastq_1` of each of up to 5 distinct samples in `{SAMPLES_CSV}`. The `zcat | head` read below is bounded (1000 reads) and intentionally light, so it is acceptable on the login node:

```bash
zcat {FASTQ_FILE} | awk 'NR%4==2 {print length($0)}' | head -n 1000 | sort -n | uniq -c | sort -rn | head -3
```

The mode of the lengths is `{READ_LENGTH}`. If lengths differ between samples, report the distribution and warn that the index is tuned to the most common length. Set `{SJDB_OVERHANG}` = `{READ_LENGTH}` − 1 and tell the user: "Detected read length {READ_LENGTH} bp -> sjdbOverhang {SJDB_OVERHANG}."

**STAR index: one index per read length and per mode.** Set `{STAR_INDEX}` to:

- F1: `{GENOME_DIR}/index/star_ase_masked_sjdb{SJDB_OVERHANG}/` (built from the third-allele masked genome, so it never collides with an unmasked index of the same read length);
- outbred: `{GENOME_DIR}/index/star_ase_sjdb{SJDB_OVERHANG}/` (built from the normal reference).

If `SA`, `Genome` and `sjdbList.out.tab` already exist and are non-empty there, reuse it; otherwise the prep job builds it at that path (STAR 2.7.10b from `{STAR_SIF}`, on a compute node). The prep job computes the genome length and `genomeSAindexNbases` = min(14, max(4, floor(log2(L)/2 - 1))) with L = genome length (for example 300,000 bp gives 8; 40,001 gives 6; 3.1e9 gives 14; 2 kb gives 4) itself, in shell at run time, because the FASTA may not exist yet when the wizard writes the script; the wizard never substitutes that value.

---

## Step 6 — Genotype source

The allelic analysis needs heterozygous SNP positions.

**F1 mode.** Ask (numbered): "1. I have a parental-difference VCF · 2. Build it from the Mouse Genomes Project (mouse helper)".

1. `{PARENTAL_VCF}`: a biallelic SNP VCF (bgzipped and tabix-indexed) where the REF allele is `{STRAIN_A}` (the reference strain) and the ALT allele is `{STRAIN_B}`. Check with `bcftools` (from `{BCFTOOLS_SIF}`, in a compute job, never on the login node) that contig names match the FASTA.
2. The mouse helper queries the Mouse Genomes Project (release REL-2112-v8, GRCm39, contigs `1`, `2`, ... without `chr`) for the two strains as a per-chromosome array job and writes `{PARENTAL_VCF}` with the same convention. It is a compute job; do not run it on the login node.

Either way the prep job also writes `f1_het_sites.vcf.gz`, a bgzipped and tabix-indexed het-sites VCF (single sample `F1`, genotype `0/1` at every parental SNP), beside the parental-difference VCF. ASEReadCounter needs this genotype column and an indexed VCF; a sites-only VCF gives zero counts.

**Outbred mode.** One VCF per `individual`, `{GENOTYPE_VCFS}` (a mapping individual -> path). Ask (numbered): "1. Matched WGS or array genotypes (preferred) · 2. The rnavar `variant_calling/` filtered VCFs".

- Option 2 requires this caveat, stated verbatim to the user: "RNA-derived genotypes are circular: heterozygous sites with strong imbalance may be called homozygous, so genotypes called from the same RNA-seq bias results toward balance and miss lowly expressed sites. External genotypes are preferred." Continue only after the user confirms.
- The prep job turns each individual's VCF into a single-sample VCF restricted to biallelic heterozygous SNPs (`bcftools view -s SAMPLE -g het -v snps -m2 -M2`, plus quality filters). The per-individual VCF given to STAR must be heterozygous-only: STAR does not ignore homozygous genotypes and uses only the first sample column.

Input VCFs are read-only; the derived files go to `{RESULTS_DIR}`.

---

## Step 7 — Analysis menu

Ask (numbered, multi-select, for example "1,3"):

1. **Per-sample allelic imbalance** (always on, cannot be deselected): reference-bias diagnostic plus per-gene and per-site allelic ratios for every sample.
2. **Reciprocal F1 analysis** (parent-of-origin versus strain effect): offered only when `{MODE}` = `f1` and both `cross_direction` values are present in `{SAMPLES_CSV}`.
3. **Differential ASE between conditions**: offered only when at least two conditions each have replicates (replicates = at least 2 samples in the condition, as defined in Step 5).
4. **phASER haplotype phasing**: offered only when `{MODE}` = `outbred`.

Show only the options whose preconditions hold, and say why any other is hidden. Store the selection as `{ANALYSES}`. Stage 1 implements only the always-on analysis: if the user chooses reciprocal, differential or phASER, say "available in a later stage" and continue with the per-sample analysis.

---

## Step 8 — Constants

Show these defaults and let the user edit any of them:

| Constant | Default | Meaning |
|---|---|---|
| `MIN_DEPTH` | 10 | minimum total reads at a site or gene for testing |
| `FDR_SIG` | 0.05 | FDR threshold for significance |
| `ABS_DEV_SIG` | 0.1 | minimum absolute deviation of the reference fraction from 0.5 to call a gene imbalanced |
| `BIAS_TOL` | 0.03 | reference-bias tolerance: a sample is flagged when its reference-fraction deviation from 0.5 exceeds `max(BIAS_TOL, 3 x SE)`, SE = sqrt(0.25 / total reads) |

ASEReadCounter defaults: `--min-mapping-quality` 10 and `--min-base-quality` 10 (the tool's own defaults are 0). Record all values; they are written into the scripts and the Rmd parameters.

---

## Step 9 — Resources

The per-sample steps (STAR alignment, duplicate marking, ASEReadCounter) run as one SLURM array job:

```bash
#SBATCH --array=1-{ARRAY_N}
```

with `{ARRAY_N}` = the number of rows in `{SAMPLES_CSV}`. Set `{MEM}` and `{TIME}` for the array job:

- human or other large genome (over 1 Gb): `-n 8 --mem=48G -t 4:00:00`;
- small genome (under 100 Mb): `-n 4 --mem=8G -t 1:00:00`;
- in between: `-n 8 --mem=32G -t 2:00:00`.

The prep jobs (index build, reference and VCF preparation) are scaled by genome size (approximate it from `stat -c %s FASTA` or an existing `.fai`, never by reading the FASTA on the login node): over 1 Gb `-n 8 --mem=64G -t 4:00:00`; under 100 Mb `-n 4 --mem=8G -t 0:30:00`; in between `-n 8 --mem=32G -t 2:00:00`. Never request more than 64 G or 4 h unless the user asks. Every job is submitted with `sbatch -p bcc`.
