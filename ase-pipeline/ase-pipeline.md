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


---

## Step 10 — Reference and genotype preparation scripts

Write every script below into `{RESULTS_DIR}/scripts/` (create `{RESULTS_DIR}/scripts`, `{RESULTS_DIR}/logs`, `{RESULTS_DIR}/tmp` first: `sbatch` silently drops output when the `-o` directory is missing, and `/tmp` is node-local, so job outputs and temporary files always go under `{RESULTS_DIR}`). Substitute the placeholders from Steps 0-9; never leave a `{...}` in a script. Every generated script starts with `set -uo pipefail` and checks the exit code of every step itself (`|| { echo "ERROR: ..." >&2; exit 1; }`), never a blind `set -e`. The prep jobs use the prep tier from Step 9 (`{PREP_RESOURCES}`, for example `-n 4 --mem=8G -t 0:30:00`) and the array jobs use `{ARRAY_RESOURCES}` (the `-n/--mem/-t` triple chosen in Step 9); every job is submitted with `sbatch -p bcc`.

| Script | Mode | Purpose |
|---|---|---|
| `prep_f1_reference.sh` | F1 | third-allele masked genome, het-sites VCF, counting reference files, STAR index |
| `prep_genotypes.sh` | outbred | per-individual heterozygous VCFs, counting reference files, STAR index |
| `extract_mgp_parental_vcf.sh` | F1, optional | mouse helper: builds `{PARENTAL_VCF}` (array job plus a dependent concat job) |
| `align_count_f1.sh`, `align_wasp_count.sh` | F1 / outbred | per-sample array job (Step 11) |

### Shared block C — container wrappers (start of every script, after the `#SBATCH` header and `set -uo pipefail`)

Tools come from cached Singularity biocontainers (`{STAR_SIF}`, `{GATK_SIF}`, `{BCFTOOLS_SIF}`, `{SAMTOOLS_SIF}`, `{PICARD_SIF}` from Step 2); never `module add` anything except singularity. `bcftools`, `bgzip` and `tabix` all come from the bcftools 1.20 container. Keep only the wrappers a script uses. `--bind` is required because `/net/...` paths are not auto-bound; `BIND` lists each directory once (Singularity prints "destination is already in the mount point list" for a repeated one), so add directories with `add_bind`:

```bash
module add singularity/3.10.4 || exit 1
command -v singularity >/dev/null || { echo "ERROR: singularity not on PATH" >&2; exit 1; }
STAR_SIF="{STAR_SIF}"; GATK_SIF="{GATK_SIF}"; BCFTOOLS_SIF="{BCFTOOLS_SIF}"; SAMTOOLS_SIF="{SAMTOOLS_SIF}"; PICARD_SIF="{PICARD_SIF}"
fetch_sif() {   # $1 = local path, $2 = verified URL (Step 2); downloads only when the file is missing
  [ -s "$1" ] && return 0
  mkdir -p "$(dirname "$1")"
  { wget -c -O "$1.part" "$2" && mv "$1.part" "$1"; } || { echo "ERROR: could not download $1" >&2; exit 1; }
}
fetch_sif "$BCFTOOLS_SIF" "https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0"
# ... one fetch_sif line per container this script uses (URLs from the Step 2 table)
BIND=""
add_bind() { case ",$BIND," in *",$1,"*) ;; *) BIND="${BIND:+$BIND,}$1" ;; esac; }
add_bind "{CWD}"; add_bind "{GENOME_DIR}"; add_bind "$(dirname "{FASTA_PATH}")"
star()     { local SIF="$STAR_SIF";     singularity exec --bind "$BIND" "$SIF" STAR "$@"; }
gatk()     { local SIF="$GATK_SIF";     singularity exec --bind "$BIND" "$SIF" gatk "$@"; }
bcftools() { local SIF="$BCFTOOLS_SIF"; singularity exec --bind "$BIND" "$SIF" bcftools "$@"; }
bgzip()    { local SIF="$BCFTOOLS_SIF"; singularity exec --bind "$BIND" "$SIF" bgzip "$@"; }
tabix()    { local SIF="$BCFTOOLS_SIF"; singularity exec --bind "$BIND" "$SIF" tabix "$@"; }
samtools() { local SIF="$SAMTOOLS_SIF"; singularity exec --bind "$BIND" "$SIF" samtools "$@"; }
picard()    { local SIF="$PICARD_SIF";  singularity exec --bind "$BIND" "$SIF" picard "$@"; }
```

Add `add_bind "$(dirname <path>)"` for every other input directory the script reads (parental or genotype VCF, FASTQ directories). If `{FASTA_PATH}` or `{GTF_PATH}` do not exist yet (Step 4 option 1), put the Step 4 download-and-decompress commands (URLs shown to the user first) before block R, with the same `wget -c -O FILE.part URL && mv FILE.part FILE || { echo "ERROR: ..." >&2; exit 1; }` pattern.

### Shared block R — reference files for ASEReadCounter (both modes; `FASTA="{FASTA_PATH}"` is the ORIGINAL, unmasked FASTA)

`gatk ASEReadCounter -R` needs `<fasta>.fai` and `<fasta>.dict` next to the FASTA, and counting must use the original FASTA so its sequence names match the BAM header. If the genome folder is read-only, the job fails with a clear message; copy the FASTA to a writable folder and point `{FASTA_PATH}` there.

```bash
FASTA="{FASTA_PATH}"
[ -s "$FASTA.fai" ] || samtools faidx "$FASTA" || { echo "ERROR: cannot create $FASTA.fai (read-only genome folder?)" >&2; exit 1; }
DICT="${FASTA%.*}.dict"
[ -s "$DICT" ] || gatk CreateSequenceDictionary -R "$FASTA" -O "$DICT" || { echo "ERROR: cannot create $DICT" >&2; exit 1; }
```

### Shared block I — STAR index (`INDEX_FASTA` = `$REF_DIR/masked.fa` in F1 mode, `$FASTA` in outbred mode; `{STAR_INDEX}` from Step 5, part 2)

```bash
[ -s "$INDEX_FASTA.fai" ] || samtools faidx "$INDEX_FASTA" || { echo "ERROR: faidx failed for $INDEX_FASTA" >&2; exit 1; }
SA_INDEX_NBASES=$(awk '{L+=$2} END{n=int(log(L)/log(2)/2-1); if(n>14)n=14; if(n<4)n=4; print n}' "$INDEX_FASTA.fai")
if [ -s "{STAR_INDEX}/SA" ] && [ -s "{STAR_INDEX}/Genome" ] && [ -s "{STAR_INDEX}/sjdbList.out.tab" ]; then
  echo "Reusing STAR index {STAR_INDEX}"
else
  mkdir -p "{STAR_INDEX}" || exit 1
  star --runMode genomeGenerate --genomeDir "{STAR_INDEX}" --genomeFastaFiles "$INDEX_FASTA" \
       --sjdbGTFfile "{GTF_PATH}" --sjdbOverhang {SJDB_OVERHANG} --genomeSAindexNbases "$SA_INDEX_NBASES" \
       --runThreadN "${SLURM_NTASKS:-4}" --outFileNamePrefix "{RESULTS_DIR}/logs/star_index_" \
    || { echo "ERROR: STAR genomeGenerate failed" >&2; exit 1; }
  [ -s "{STAR_INDEX}/SA" ] || { echo "ERROR: STAR index incomplete in {STAR_INDEX}" >&2; exit 1; }
fi
```

STAR prints "Could not move Log.out" for the index build when `--outFileNamePrefix` is set; it is harmless (the log stays in `{RESULTS_DIR}/logs`). The suffix-array parameter is computed in shell from the `.fai` (light; never read the FASTA on the login node), so the wizard never substitutes it.

### `prep_f1_reference.sh` (F1 mode)

Header: `#!/bin/bash`, `#SBATCH -N 1 -p bcc`, `#SBATCH {PREP_RESOURCES}`, `#SBATCH --mail-type=END,FAIL`, `#SBATCH --mail-user={USER_EMAIL}`, `#SBATCH -o {RESULTS_DIR}/logs/prep_f1_reference_%j.out`, `set -uo pipefail`, block C (STAR, GATK, BCFTOOLS, SAMTOOLS; plus `add_bind "$(dirname "{PARENTAL_VCF}")"`), then:

```bash
REF_DIR="{RESULTS_DIR}/reference"; mkdir -p "$REF_DIR" || exit 1
FASTA="{FASTA_PATH}"; PARENTAL="{PARENTAL_VCF}"
[ -s "$FASTA.fai" ] || samtools faidx "$FASTA" || { echo "ERROR: cannot create $FASTA.fai" >&2; exit 1; }

# 1. biallelic SNP sites of the parental VCF (REF = strain A, ALT = strain B), sites only
bcftools view -G -m2 -M2 -v snps -O v -o "$REF_DIR/parental_snps.sites.vcf" "$PARENTAL" \
  || { echo "ERROR: bcftools view failed on $PARENTAL" >&2; exit 1; }
N_SITES=$(grep -vc '^#' "$REF_DIR/parental_snps.sites.vcf")
[ "$N_SITES" -gt 0 ] || { echo "ERROR: no biallelic SNP sites in $PARENTAL (contig names must match the FASTA)" >&2; exit 1; }
echo "Parental biallelic SNP sites: $N_SITES"

# 2. strain A must be the reference strain: the VCF REF must equal the FASTA base at every site
bcftools norm -f "$FASTA" -c e -o /dev/null "$REF_DIR/parental_snps.sites.vcf" \
  || { echo "ERROR: the REF allele of $PARENTAL differs from the FASTA; F1 mode needs strain A = the reference strain, REF = strain A, ALT = strain B (and matching contig names)" >&2; exit 1; }

# 3. third-allele masking: at each site the masked genome carries the first of A, C, G, T that is
#    neither the REF nor the ALT allele (same rule everywhere), so neither strain is favoured in mapping
awk 'BEGIN{OFS="\t"}
  function third(r,a,  i,c,s){ s="ACGT"; for(i=1;i<=4;i++){ c=substr(s,i,1); if(c!=r && c!=a) return c } }
  /^##/ {print; next}
  /^#/  {print; next}
  { print $1,$2,$3,$4,third(toupper($4),toupper($5)),".",".","." }' \
  "$REF_DIR/parental_snps.sites.vcf" > "$REF_DIR/masked_sites.vcf" || { echo "ERROR: masked sites failed" >&2; exit 1; }
awk 'BEGIN{OFS="\t"} !/^#/ {print $1,$2-1,$2,$4 ">" $5}' "$REF_DIR/masked_sites.vcf" > "$REF_DIR/masked_positions.bed"
bgzip -c "$REF_DIR/masked_sites.vcf" > "$REF_DIR/masked_sites.vcf.gz" && tabix -f -p vcf "$REF_DIR/masked_sites.vcf.gz" \
  || { echo "ERROR: bgzip/tabix failed for masked_sites.vcf" >&2; exit 1; }
bcftools consensus -H A -f "$FASTA" "$REF_DIR/masked_sites.vcf.gz" > "$REF_DIR/masked.fa" \
  || { echo "ERROR: bcftools consensus failed" >&2; exit 1; }

# 4. VERIFY the masked genome: it must differ from the FASTA at exactly the masked sites, and at every
#    site the masked base must be the chosen third allele (cmp -l lists differing bytes; the .fai
#    layouts must be identical so that a byte offset maps to one contig position)
samtools faidx "$REF_DIR/masked.fa" || { echo "ERROR: faidx failed for masked.fa" >&2; exit 1; }
cmp -s <(cut -f1-5 "$FASTA.fai") <(cut -f1-5 "$REF_DIR/masked.fa.fai") \
  || { echo "ERROR: masked.fa has a different contig layout from $FASTA" >&2; exit 1; }
VERIFY=$(cmp -l "$FASTA" "$REF_DIR/masked.fa" 2> "$REF_DIR/cmp.err" | awk -v FAI="$FASTA.fai" -v SITES="$REF_DIR/masked_sites.vcf" '
  BEGIN { nc=0
    while ((getline line < FAI) > 0) { split(line, f, "\t"); nc++; name[nc]=f[1]; off[nc]=f[3]+0; lb[nc]=f[4]+0; lw[nc]=f[5]+0 }
    while ((getline line < SITES) > 0) { if (line ~ /^#/) continue; split(line, f, "\t"); k=f[1] ":" f[2]; if (!(k in want)) nwant++; want[k]=f[5]; nrec++ }
    oct["101"]="A"; oct["103"]="C"; oct["107"]="G"; oct["124"]="T"; oct["141"]="A"; oct["143"]="C"; oct["147"]="G"; oct["164"]="T" }
  { b=$1-1; lo=1; hi=nc
    while (lo<hi) { mid=int((lo+hi+1)/2); if (off[mid]<=b) lo=mid; else hi=mid-1 }
    i=b-off[lo]; pos=int(i/lw[lo])*lb[lo] + (i%lw[lo]) + 1; k=name[lo] ":" pos; ndiff++
    if ((k in want) && oct[$3]==want[k]) ok++; else bad++ }
  END { printf "%d %d %d %d %d\n", ndiff+0, ok+0, bad+0, nwant+0, nrec+0 }')
read -r N_DIFF N_OK N_BAD N_WANT N_REC <<< "$VERIFY"
echo "Masked genome check: differing positions=$N_DIFF, matching the third allele=$N_OK, wrong=$N_BAD, distinct sites=$N_WANT, site records=$N_REC"
if grep -q EOF "$REF_DIR/cmp.err" || [ "$N_DIFF" -ne "$N_REC" ] || [ "$N_WANT" -ne "$N_REC" ] || [ "$N_OK" -ne "$N_DIFF" ] || [ "$N_BAD" -ne 0 ]; then
  echo "ERROR: masked.fa does not match the third allele masking of the $N_REC sites; stopping" >&2; exit 1
fi

# 5. het-sites VCF for ASEReadCounter: single sample F1, genotype 0/1 at every parental SNP
#    (bgzipped + tabix; a sites-only VCF gives 0 rows and homozygous sites are skipped)
{ echo "##fileformat=VCFv4.2"
  awk '{printf "##contig=<ID=%s,length=%s>\n", $1, $2}' "$FASTA.fai"
  echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tF1\n'
  awk 'BEGIN{OFS="\t"} !/^#/ {print $1,$2,".",$4,$5,".",".",".","GT","0/1"}' "$REF_DIR/parental_snps.sites.vcf"
} > "$REF_DIR/f1_het_sites.vcf" || { echo "ERROR: f1_het_sites.vcf failed" >&2; exit 1; }
bgzip -c "$REF_DIR/f1_het_sites.vcf" > "$REF_DIR/f1_het_sites.vcf.gz" && tabix -f -p vcf "$REF_DIR/f1_het_sites.vcf.gz" \
  || { echo "ERROR: bgzip/tabix failed for f1_het_sites.vcf" >&2; exit 1; }
```

`-H A` makes `bcftools consensus` apply the ALT (third) allele of the sites-only VCF; without it nothing is applied and the check in step 4 stops the job. Then block R (`.fai` and `.dict` of the original FASTA), then block I with `INDEX_FASTA="$REF_DIR/masked.fa"`, then `echo "prep_f1_reference done"`. The array job later counts against the ORIGINAL FASTA with `f1_het_sites.vcf.gz`. Third-allele masking example: REF `A`, ALT `G` gives `C`; REF `A`, ALT `C` gives `G`. Sites are the parental SNPs only (no indels), and `masked.fa` is only used to build the STAR index.

### `prep_genotypes.sh` (outbred mode)

Header as above (`prep_genotypes_%j.out`), block C (STAR, GATK, BCFTOOLS, SAMTOOLS plus `add_bind` for each genotype VCF directory), then one heterozygous single-sample VCF per individual. The wizard writes the `MAP` lines (individual, path, sample name inside that VCF) from `{GENOTYPE_VCFS}` and shows the two filter settings for editing: rnavar VCFs use `FILTER_EXPR='FMT/DP>=10'`; external genotypes use a genotype-quality filter such as `FILTER_EXPR='FMT/GQ>=20'`; an empty value skips the expression filter (test data). `-f PASS,.` keeps PASS and unfiltered records.

```bash
GENO_DIR="{RESULTS_DIR}/genotypes"; mkdir -p "$GENO_DIR" || exit 1
FILTER_EXPR='FMT/DP>=10'
MIN_SITES_WARN=1000      # human; use 20 for small test data
cat > "$GENO_DIR/genotype_map.tsv" <<'MAP'
{INDIVIDUAL}	{GENOTYPE_VCF}	{SAMPLE_IN_VCF}
MAP
while IFS=$'\t' read -r IND VCF VSAMPLE; do
  [ -n "$IND" ] || continue
  OUT="$GENO_DIR/$IND.het"
  bcftools view -s "$VSAMPLE" -g het -v snps -m2 -M2 -f PASS,. -O z -o "$OUT.stage1.vcf.gz" "$VCF" \
    || { echo "ERROR: bcftools view failed for individual $IND ($VCF, sample $VSAMPLE)" >&2; exit 1; }
  if [ -n "$FILTER_EXPR" ]; then
    bcftools view -i "$FILTER_EXPR" -O z -o "$OUT.vcf.gz" "$OUT.stage1.vcf.gz" \
      || { echo "ERROR: filter '$FILTER_EXPR' failed for individual $IND" >&2; exit 1; }
  else
    mv "$OUT.stage1.vcf.gz" "$OUT.vcf.gz" || exit 1
  fi
  rm -f "$OUT.stage1.vcf.gz"
  tabix -f -p vcf "$OUT.vcf.gz" || { echo "ERROR: tabix failed for $OUT.vcf.gz" >&2; exit 1; }
  # plain-text copy for STAR --varVCFfile (heterozygous-only, single sample: STAR uses only the first sample)
  bcftools view -O v -o "$OUT.vcf" "$OUT.vcf.gz" || { echo "ERROR: cannot write $OUT.vcf" >&2; exit 1; }
  N=$(bcftools view -H "$OUT.vcf.gz" | wc -l)
  echo "Individual $IND: $N heterozygous SNP sites"
  [ "$N" -gt 0 ] || { echo "ERROR: individual $IND has no heterozygous SNPs after filtering" >&2; exit 1; }
  [ "$N" -ge "$MIN_SITES_WARN" ] || echo "WARNING: individual $IND has only $N heterozygous sites (fewer than $MIN_SITES_WARN)"
done < "$GENO_DIR/genotype_map.tsv"
```

The array job reads `{RESULTS_DIR}/genotypes/{individual}.het.vcf` (STAR) and `{individual}.het.vcf.gz` (ASEReadCounter; a plain VCF fails there). Then block R and block I with `INDEX_FASTA="$FASTA"` (the unmasked genome), then `echo "prep_genotypes done"`.

### `extract_mgp_parental_vcf.sh` (F1 mode, mouse helper; optional)

For the two named strains (`{STRAIN_A}` = the reference strain, for example `C57BL_6NJ`, and `{STRAIN_B}`, for example `A_J`; names exactly as in the VCF header) the helper queries the Mouse Genomes Project VCF remotely. Remote access costs about 80 s fixed overhead plus about 80 s per 10 Mb of region, so a whole-genome serial run would take about 6 hours: it is therefore **one array task per chromosome** (each well under 4 h) followed by a final dependent `bcftools concat` job. The wizard must verify the URL with a HEAD request in the session before writing it (this was verified on 2026-09-29: HTTP 200, `ase-pipeline/tests/fixtures/verified_urls.txt`):

```bash
curl -sI "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz" | head -1
```

The release is GRCm39 and its contigs are `1`, `2`, ... (no `chr` prefix), so the FASTA must be GRCm39 with the same names; the script stops on a first-contig mismatch (the same guard as in `nfcore-rnavar-setup`). The wizard writes `{RESULTS_DIR}/mgp/chromosomes.txt` (one FASTA chromosome per line, taken from the `.fai`: `1` ... `19`, `X`) and sets `{N_CHROM}` to its number of lines. Header: `#SBATCH -N 1 -p bcc`, `#SBATCH --array=1-{N_CHROM}`, `#SBATCH -n 2 --mem=8G -t 4:00:00`, mail lines, `#SBATCH -o {RESULTS_DIR}/logs/extract_mgp_%A_%a.out`, `set -uo pipefail`, block C (BCFTOOLS only, plus `add_bind "{RESULTS_DIR}"`), then:

```bash
URL="https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz"
A="{STRAIN_A}"; B="{STRAIN_B}"
MGP_DIR="{RESULTS_DIR}/mgp"; OUT_VCF="{PARENTAL_VCF}"      # {RESULTS_DIR}/mgp/parental_{STRAIN_A}_{STRAIN_B}.vcf.gz
if [ "${1:-}" = "concat" ]; then                            # final dependent job: sbatch --dependency=afterok:<array_jobid> extract_mgp_parental_vcf.sh concat
  : > "$MGP_DIR/concat_list.txt"
  while read -r CHR; do
    [ -s "$MGP_DIR/$CHR.parental.vcf.gz.tbi" ] || { echo "ERROR: missing per-chromosome result for $CHR" >&2; exit 1; }
    echo "$MGP_DIR/$CHR.parental.vcf.gz" >> "$MGP_DIR/concat_list.txt"
  done < "$MGP_DIR/chromosomes.txt"
  bcftools concat -f "$MGP_DIR/concat_list.txt" -O z -o "$OUT_VCF" && tabix -f -p vcf "$OUT_VCF" \
    || { echo "ERROR: bcftools concat failed" >&2; exit 1; }
  echo "Parental VCF: $(bcftools view -H "$OUT_VCF" | wc -l) sites"; exit 0
fi
CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MGP_DIR/chromosomes.txt")
[ -n "$CHR" ] || { echo "ERROR: no chromosome for array task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
# contig guard: first contig of the MGP header versus first FASTA header
VCF_CONTIG=$(bcftools view -h "$URL" | awk -F'[=,>]' '/^##contig/ {print $3; exit}')
FASTA_CONTIG=$(grep -m1 '^>' "{FASTA_PATH}" | cut -d' ' -f1 | sed 's/^>//')
if [ -z "$VCF_CONTIG" ] || [ -z "$FASTA_CONTIG" ] || [ "$VCF_CONTIG" != "$FASTA_CONTIG" ]; then
  echo "ERROR: contig mismatch (MGP: '$VCF_CONTIG', FASTA: '$FASTA_CONTIG'); use a GRCm39 FASTA with names 1, 2, ..." >&2; exit 1
fi
RAW="$MGP_DIR/$CHR.raw.vcf.gz"
bcftools view -r "$CHR" -s "$A,$B" -f PASS -m2 -M2 -v snps -O z -o "$RAW" "$URL" \
  || { echo "ERROR: remote query failed for chromosome $CHR" >&2; exit 1; }
IA=$(bcftools query -l "$RAW" | awk -v s="$A" '$0==s {print NR-1}')
IB=$(bcftools query -l "$RAW" | awk -v s="$B" '$0==s {print NR-1}')
[ -n "$IA" ] && [ -n "$IB" ] || { echo "ERROR: strain '$A' or '$B' not found in the VCF header" >&2; exit 1; }
# keep sites where the two strains differ and both have FMT/FI=1; strain A must carry the reference allele
# (parental VCF convention: REF = strain A = FASTA base), strain B the alternative allele
bcftools view -G -i "FMT/FI[$IA]=1 && FMT/FI[$IB]=1 && GT[$IA]=\"0/0\" && GT[$IB]=\"1/1\"" -O z -o "$MGP_DIR/$CHR.parental.vcf.gz" "$RAW" \
  && tabix -f -p vcf "$MGP_DIR/$CHR.parental.vcf.gz" || { echo "ERROR: filtering failed for chromosome $CHR" >&2; exit 1; }
echo "Chromosome $CHR: $(bcftools view -H "$RAW" | wc -l) PASS SNPs, $(bcftools view -H "$MGP_DIR/$CHR.parental.vcf.gz" | wc -l) kept"
```

Submit: `J=$(sbatch -p bcc --parsable extract_mgp_parental_vcf.sh)`, then `sbatch -p bcc --dependency=afterok:$J extract_mgp_parental_vcf.sh concat`, and `prep_f1_reference.sh` only after the concat job (`--dependency=afterok:<concat_jobid>`). Sites where strain A carries the alternative allele (the reference genome is C57BL/6J, so a substrain can differ from it) are dropped: the masked reference and the REF = strain A convention require the FASTA base to be the strain A allele.

---

## Step 11 — Per-sample array scripts

Both scripts are SLURM array jobs (`#SBATCH --array=1-{ARRAY_N}`, `#SBATCH {ARRAY_RESOURCES}`), take row `SLURM_ARRAY_TASK_ID` of `{SAMPLES_CSV}` (data row 1 = task 1), and write the outputs below. They use block C, but the containers must already be cached by the prep job, so parallel tasks never download the same file: in these scripts define `fetch_sif` as a check only, `fetch_sif() { [ -s "$1" ] || { echo "ERROR: missing container $1 (run the prep job first)" >&2; exit 1; }; }`. `{MIN_MAPQ}` and `{MIN_BASEQ}` are the ASEReadCounter thresholds recorded in Step 8.

| Output | Content |
|---|---|
| `{RESULTS_DIR}/bam/{sample}.bam` (+ `.bai`) | coordinate-sorted, duplicates marked (not removed: ASEReadCounter skips them), read group kept |
| `{RESULTS_DIR}/ase_counts/{sample}.table` | ASEReadCounter table: `contig position variantID refAllele altAllele refCount altCount totalCount lowMAPQDepth lowBaseQDepth rawDepth otherBases improperPairs` |
| `{RESULTS_DIR}/ase_counts/{sample}.wasp_stats.tsv` | outbred only: alignments by `vW` and first `vA` value before filtering |

**Read group is mandatory.** STAR writes no read group by default, and ASEReadCounter's read-group filter then silently drops every read (empty table, exit 0). Both scripts therefore pass `--outSAMattrRGline ID:$SAMPLE SM:$SAMPLE PL:ILLUMINA`; Picard keeps the read group.

**Common start of both scripts** (after the `#SBATCH` header, `set -uo pipefail` and block C with STAR, GATK, SAMTOOLS and PICARD wrappers):

```bash
ROW=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1' "{SAMPLES_CSV}" | tr -d '\r')
[ -n "$ROW" ] || { echo "ERROR: no row $SLURM_ARRAY_TASK_ID in {SAMPLES_CSV}" >&2; exit 1; }
IFS=, read -r SAMPLE FQ1 FQ2 CONDITION CROSS INDIVIDUAL <<< "$ROW"
case "$FQ1" in /*) ;; *) FQ1="{CWD}/$FQ1" ;; esac
[ -z "$FQ2" ] || case "$FQ2" in /*) ;; *) FQ2="{CWD}/$FQ2" ;; esac
add_bind "$(dirname "$FQ1")"; [ -z "$FQ2" ] || add_bind "$(dirname "$FQ2")"
die() { echo "ERROR: sample $SAMPLE: $*" >&2; exit 1; }
[ -s "$FQ1" ] || die "missing $FQ1"
FASTA="{FASTA_PATH}"; R="{RESULTS_DIR}"
STAR_DIR="$R/star/$SAMPLE"; TMPD="$R/tmp/$SAMPLE"; BAM="$R/bam/$SAMPLE.bam"; TABLE="$R/ase_counts/$SAMPLE.table"
mkdir -p "$STAR_DIR" "$TMPD" "$R/bam" "$R/ase_counts" || die "cannot create output directories"
rm -rf "$STAR_DIR/_STARtmp"
```

### `align_count_f1.sh` (F1 mode)

After the common start (`--readFilesIn` takes one file for single-end data and two for paired-end data; the wizard writes `"$FQ1" "$FQ2"` or only `"$FQ1"` accordingly):

```bash
star --runThreadN "${SLURM_NTASKS:-4}" --genomeDir "{STAR_INDEX}" \
     --readFilesIn "$FQ1" "$FQ2" --readFilesCommand zcat \
     --outFileNamePrefix "$STAR_DIR/" --outSAMtype BAM SortedByCoordinate \
     --outSAMattrRGline ID:$SAMPLE SM:$SAMPLE PL:ILLUMINA || die "STAR failed"
RAW="$STAR_DIR/Aligned.sortedByCoord.out.bam"
[ -s "$RAW" ] || die "STAR wrote no BAM"
picard MarkDuplicates I="$RAW" O="$BAM" M="$R/bam/$SAMPLE.markdup_metrics.txt" \
       TMP_DIR="$TMPD" REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT || die "MarkDuplicates failed"
samtools index "$BAM" || die "samtools index failed"
gatk ASEReadCounter -R "$FASTA" -I "$BAM" -V "$R/reference/f1_het_sites.vcf.gz" \
     --min-mapping-quality {MIN_MAPQ} --min-base-quality {MIN_BASEQ} \
     --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE \
     --tmp-dir "$TMPD" -O "$TABLE" || die "ASEReadCounter failed"
```

Picard 3.1.1 accepts the legacy `KEY=VALUE` form used above (it prints a syntax-change notice; verified in the 3.1.1 container). `STAR --outSAMtype BAM SortedByCoordinate` already sorts, so no separate `samtools sort` step is needed. The BAM header must match the ORIGINAL FASTA (`masked.fa` only changes bases, never names).

### `align_wasp_count.sh` (outbred mode)

After the common start (also `[ -n "$INDIVIDUAL" ] || die "empty individual"`), with `HET_PLAIN="$R/genotypes/$INDIVIDUAL.het.vcf"` and `HET_GZ="$R/genotypes/$INDIVIDUAL.het.vcf.gz"` (both from `prep_genotypes.sh`; the STAR file must be the heterozygous-only plain VCF):

```bash
[ -s "$HET_PLAIN" ] && [ -s "$HET_GZ" ] || die "missing genotypes for individual $INDIVIDUAL (run prep_genotypes.sh)"
star --runThreadN "${SLURM_NTASKS:-4}" --genomeDir "{STAR_INDEX}" \
     --readFilesIn "$FQ1" "$FQ2" --readFilesCommand zcat \
     --outFileNamePrefix "$STAR_DIR/" --varVCFfile "$HET_PLAIN" --waspOutputMode SAMtag \
     --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM vA vG vW \
     --outSAMattrRGline ID:$SAMPLE SM:$SAMPLE PL:ILLUMINA || die "STAR failed"
RAW="$STAR_DIR/Aligned.sortedByCoord.out.bam"
[ -s "$RAW" ] || die "STAR wrote no BAM"

# WASP statistics from the STAR BAM BEFORE filtering: alignments by vW value (none = no tag; 1 pass;
# 2 multi-mapping; 3 variant base N; 4 remap failed; 5 remap multi-maps; 6 remap to a different locus;
# 7 too many variants) and by the first vA value (1 ref, 2 alt, 3 no match; - = no tag)
STATS="$R/ase_counts/$SAMPLE.wasp_stats.tsv"
printf 'vW\tvA\tn\n' > "$STATS" || die "cannot write $STATS"
samtools view "$RAW" | awk 'BEGIN{OFS="\t"} { w="none"; a="-"
    for (i=12; i<=NF; i++) { if ($i ~ /^vW:i:/) w=substr($i,6); else if ($i ~ /^vA:B:c,/) { split($i,p,","); a=p[2] } }
    n[w OFS a]++ } END { for (k in n) print k, n[k] }' | sort -k1,1 -k2,2 >> "$STATS" || die "WASP statistics failed"

# keep alignments with vW:i:1 and alignments without a vW tag; drop vW 2-7
FILT="$STAR_DIR/wasp_filtered.bam"
samtools view -h "$RAW" | awk '/^@/ {print; next} { keep=1
    for (i=12; i<=NF; i++) if ($i ~ /^vW:i:/ && $i != "vW:i:1") { keep=0; break }
    if (keep) print }' | samtools view -b -o "$FILT" - || die "vW filtering failed"
picard MarkDuplicates I="$FILT" O="$BAM" M="$R/bam/$SAMPLE.markdup_metrics.txt" \
       TMP_DIR="$TMPD" REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT || die "MarkDuplicates failed"
samtools index "$BAM" || die "samtools index failed"
gatk ASEReadCounter -R "$FASTA" -I "$BAM" -V "$HET_GZ" \
     --min-mapping-quality {MIN_MAPQ} --min-base-quality {MIN_BASEQ} \
     --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE \
     --tmp-dir "$TMPD" -O "$TABLE" || die "ASEReadCounter failed"
```

Never use `--outSAMattributes` values other than the explicit list above: WASP needs `vA vG vW` spelled out. `vG` is 0-based. WASP does not promise less bias; the stats file lets Rmd 01 show how many alignments were removed and from which allele.

### Empty-table guard (last block of both scripts)

The tool exits 0 even when it counted nothing, so the script itself fails, and an empty table can never reach Rmd 01:

```bash
# ase_counts/$SAMPLE.table needs non-empty data rows (header only = ASEReadCounter counted nothing)
N_ROWS=0
[ -s "$TABLE" ] && N_ROWS=$(awk 'NR>1' "$TABLE" | wc -l)
if [ "$N_ROWS" -lt 1 ]; then
  echo "ERROR: sample $SAMPLE: ase_counts/$SAMPLE.table is missing or has no data rows (check the read group, the het-sites VCF and the contig names)" >&2
  exit 1
fi
echo "Sample $SAMPLE: $N_ROWS sites counted"
```

### Submission order

Write scripts to `{RESULTS_DIR}/scripts/`, run `mkdir -p {RESULTS_DIR}/logs {RESULTS_DIR}/tmp`, then submit with `--parsable` and `--dependency=afterok` so the array never starts before its inputs exist (when the mouse helper is used, the prep job depends on its concat job):

```bash
P=$(sbatch -p bcc --parsable {RESULTS_DIR}/scripts/prep_f1_reference.sh)      # outbred: prep_genotypes.sh
sbatch -p bcc --dependency=afterok:$P {RESULTS_DIR}/scripts/align_count_f1.sh  # outbred: align_wasp_count.sh
```

Wait for the jobs with a bounded loop and read the job logs; on a failure show the error line and stop. Rmd 01 (next stage) reads `{RESULTS_DIR}/ase_counts/*.table`, `*.wasp_stats.tsv` and `{SAMPLES_CSV}`.
