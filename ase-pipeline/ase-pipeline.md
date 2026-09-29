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
#SBATCH -n 1 --mem=4G -t 0:10:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={USER_EMAIL}
#SBATCH -o check_r_packages_%j.out
module add singularity/3.10.4 || { echo "ERROR: cannot load singularity" >&2; exit 1; }
singularity exec --bind {CWD} /net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif \
  Rscript -e 'for (p in c("aod","lme4","openxlsx","tidyverse","GenomicRanges","rtracklayer")) cat(p, requireNamespace(p, quietly = TRUE), "\n")'
```

Mandatory packages: `aod`, `lme4`, `openxlsx`, `tidyverse`, `GenomicRanges`, `rtracklayer`. Read the `.out` file when the job finishes. If any mandatory package prints `FALSE`, **stop** and ask the user to rebuild the `bulkrnaseq` image; do not substitute another model or package. (`VGAM`, `glmmTMB`, `MBASED` and `VariantAnnotation` are absent from the image and must never be required.)

---

## Step 3 — Mode

Ask the user (numbered options):

"Which experimental design is this?
1. F1 cross
2. Outbred / human"

- **1. F1 cross** — the animals are F1 offspring of two inbred parental strains. Every SNP that differs between the strains is heterozygous by construction, so genotypes come from a parental-difference VCF (or the mouse helper), not from the animals. Reads are mapped to a third-allele masked reference so neither strain is favoured in mapping.
- **2. Outbred / human** — each individual has its own genotype. Heterozygous SNPs come from a per-individual VCF (external WGS/array data, or the rnavar VCF with a stated caveat). Reads are mapped to the normal reference and STAR's WASP filter removes reads whose mapping depends on the allele they carry.

Store the answer as `{MODE}` = `f1` or `outbred`. Every later step that differs between the modes branches on `{MODE}`.
