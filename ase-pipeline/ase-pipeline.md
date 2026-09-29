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
| `RHO_MIN` | 0.01 | F1 gene-level tests only: floor for the per-sample overdispersion (`rho_gene = max(rho_corrected, RHO_MIN)`), so a noisy or zero estimate can never make the gene tests anti-conservative |

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
| `extract_mgp_parental_vcf.sh` | F1, optional | mouse helper: builds `{PARENTAL_VCF}` (array job plus a separate dependent concat job, `concat_mgp_parental_vcf.sh`) |
| `concat_mgp_parental_vcf.sh` | F1, optional | dependent single job that concatenates the per-chromosome results (no array header) |
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

# 4. VERIFY the masked genome: it must differ from the original FASTA at exactly the masked sites, and at
#    every site the masked base must be the chosen third allele. Both FASTAs are first normalised (bare-name
#    headers, 60 bases per line) so that any line width or descriptive header (Ensembl, NCBI, UCSC) works,
#    then cmp -l lists the differing bytes and the normalised .fai maps a byte offset to contig:position.
#    FASTAs that are soft-masked (lowercase blocks) work too: the octal map below accepts both cases.
samtools faidx "$REF_DIR/masked.fa" || { echo "ERROR: faidx failed for masked.fa" >&2; exit 1; }
NORM_A="$REF_DIR/verify_original.fa"; NORM_B="$REF_DIR/verify_masked.fa"
samtools faidx -n 60 -o "$NORM_A" "$FASTA" $(cut -f1 "$FASTA.fai") \
  || { echo "ERROR: cannot normalise $FASTA" >&2; exit 1; }
samtools faidx -n 60 -o "$NORM_B" "$REF_DIR/masked.fa" $(cut -f1 "$REF_DIR/masked.fa.fai") \
  || { echo "ERROR: cannot normalise masked.fa" >&2; exit 1; }
samtools faidx "$NORM_A" && samtools faidx "$NORM_B" || { echo "ERROR: faidx failed on the normalised FASTAs" >&2; exit 1; }
cmp -s <(cut -f1-5 "$NORM_A.fai") <(cut -f1-5 "$NORM_B.fai") \
  || { echo "ERROR: masked.fa has different contigs or lengths from $FASTA" >&2; exit 1; }
VERIFY=$(cmp -l "$NORM_A" "$NORM_B" 2> "$REF_DIR/cmp.err" | awk -v FAI="$NORM_A.fai" -v SITES="$REF_DIR/masked_sites.vcf" '
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
rm -f "$NORM_A" "$NORM_A.fai" "$NORM_B" "$NORM_B.fai"

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
MGP_DIR="{RESULTS_DIR}/mgp"; OUT_VCF="{PARENTAL_VCF}"      # {RESULTS_DIR}/mgp/parental_{STRAIN_A}_{STRAIN_B}.vcf.gz (written by the concat job)
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

Submit the array job first, then the concat job as its own script (`concat_mgp_parental_vcf.sh` has no `--array` header, so exactly one task runs and nothing races on `{PARENTAL_VCF}`): `J=$(sbatch -p bcc --parsable extract_mgp_parental_vcf.sh)`, `C=$(sbatch -p bcc --parsable --dependency=afterok:$J concat_mgp_parental_vcf.sh)`, and `prep_f1_reference.sh` only after the concat job (`--dependency=afterok:$C`).

**`concat_mgp_parental_vcf.sh`** header: `#SBATCH -N 1 -p bcc`, `#SBATCH -n 1 --mem=4G -t 1:00:00`, mail lines, `#SBATCH -o {RESULTS_DIR}/logs/concat_mgp_%j.out`, `set -uo pipefail`, block C (BCFTOOLS only, plus `add_bind "{RESULTS_DIR}"`), then:

```bash
MGP_DIR="{RESULTS_DIR}/mgp"; OUT_VCF="{PARENTAL_VCF}"
: > "$MGP_DIR/concat_list.txt"
while read -r CHR; do
  [ -s "$MGP_DIR/$CHR.parental.vcf.gz.tbi" ] || { echo "ERROR: missing per-chromosome result for $CHR" >&2; exit 1; }
  echo "$MGP_DIR/$CHR.parental.vcf.gz" >> "$MGP_DIR/concat_list.txt"
done < "$MGP_DIR/chromosomes.txt"
bcftools concat -f "$MGP_DIR/concat_list.txt" -O z -o "$OUT_VCF" && tabix -f -p vcf "$OUT_VCF" \
  || { echo "ERROR: bcftools concat failed" >&2; exit 1; }
echo "Parental VCF $OUT_VCF: $(bcftools view -H "$OUT_VCF" | wc -l) sites kept"
```

**Wizard instructions.** Before submission, tell a C57BL_6NJ x A_J user: "Sites where C57BL_6NJ differs from the GRCm39 (C57BL/6J) assembly are dropped, because the masked reference needs REF = strain A = the assembly base." After the helper finishes, report the number of kept sites (the last line of the concat log) to the user.

Sites where strain A carries the alternative allele (the reference genome is C57BL/6J, so a substrain can differ from it) are dropped: the masked reference and the REF = strain A convention require the FASTA base to be the strain A allele.

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
FASTA="{FASTA_PATH}"; R="{RESULTS_DIR}"
STAR_DIR="$R/star/$SAMPLE"; TMPD="$R/tmp/$SAMPLE"; BAM="$R/bam/$SAMPLE.bam"; TABLE="$R/ase_counts/$SAMPLE.table"; STATS="$R/ase_counts/$SAMPLE.wasp_stats.tsv"
mkdir -p "$STAR_DIR" "$TMPD" "$R/bam" "$R/ase_counts" || die "cannot create output directories"
rm -rf "$STAR_DIR/_STARtmp"
rm -f "$TABLE" "$STATS"      # never keep a table or stats file from an earlier run: a failed sample must have none
[ -s "$FQ1" ] || die "missing $FQ1"
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
  rm -f "$TABLE"          # a header-only table must not stay on disk for Rmd 01 to glob
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

Wait for the jobs with a bounded loop and read the job logs; on a failure show the error line and stop. Rmd 01 (Step 12) reads `{RESULTS_DIR}/ase_counts/{sample}.table` (and `{sample}.wasp_stats.tsv` in outbred mode) for every sample in `{SAMPLES_CSV}`.

---

## Step 12 — Rmd 01: import, site filters and reference-bias QC

Write `{CWD}/{TODAY_YYMMDD}_{WD_NAME}_01_import_qc.Rmd` (ask once for `{AUTHOR}` and `{PROJECT_TITLE}`; the same two values go into Rmd 02). The Rmd is self-contained (it sources no helper files), uses `knitr::opts_chunk$set(cache = FALSE)` and `options(scipen = 9)`, loads the Bioconductor packages before `tidyverse`, and calls every dplyr verb with the `dplyr::` prefix. Substitute the placeholders from Steps 0-9 (`{MODE}`, `{STRAIN_A}`, `{STRAIN_B}`, `{SAMPLES_CSV}`, `{RESULTS_DIR}`, `{GTF_PATH}`, and the four constants of Step 8 as bare numbers); for outbred mode leave `{STRAIN_A}` and `{STRAIN_B}` as the words `REF` and `ALT`. The constants block below is the only place the four thresholds are defined; every table and figure uses these objects.

**It reads `ase_counts/{sample}.table` for every sample listed in `samples.csv`, never by globbing `*.table`.** A stray or stale table from an earlier run is therefore never accepted, and the Rmd stops with an error that names every sample whose table is missing, header-only or without data rows (the Step 11 scripts delete failed tables; this is the second line of defence). In outbred mode the same rule applies to `ase_counts/{sample}.wasp_stats.tsv`.

````rmd
---
title: "{PROJECT_TITLE} - ASE import and QC"
author: "{AUTHOR}"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(cache = FALSE, echo = TRUE, message = FALSE, warning = FALSE, fig.width = 9, fig.height = 6)
options(scipen = 9)
library(GenomicRanges); library(rtracklayer); library(openxlsx)   # Bioconductor first
library(tidyverse)
```

## Constants

```{r constants}
MODE        <- "{MODE}"                 # "f1" or "outbred"
STRAIN_A    <- "{STRAIN_A}"             # F1: REF count = strain A; outbred: REF
STRAIN_B    <- "{STRAIN_B}"             # F1: ALT count = strain B; outbred: ALT
SAMPLES_CSV <- "{SAMPLES_CSV}"
RESULTS_DIR <- "{RESULTS_DIR}"
COUNTS_DIR  <- file.path(RESULTS_DIR, "ase_counts")
GTF_PATH    <- "{GTF_PATH}"
MIN_DEPTH   <- {MIN_DEPTH}              # minimum total reads at a site
FDR_SIG     <- {FDR_SIG}                # BH threshold (used in Rmd 02)
ABS_DEV_SIG <- {ABS_DEV_SIG}            # minimum |ALT fraction - 0.5| (used in Rmd 02)
BIAS_TOL    <- {BIAS_TOL}               # reference-bias tolerance
stopifnot(MODE %in% c("f1", "outbred"))
ref_label <- STRAIN_A; alt_label <- STRAIN_B
if (!file.exists(GTF_PATH)) stop("GTF not found: ", GTF_PATH, call. = FALSE)
```

## Samples and count tables

```{r read}
samples <- read.csv(SAMPLES_CSV, stringsAsFactors = FALSE, colClasses = "character")
need_cols <- c("sample", "fastq_1", "fastq_2", "condition", "cross_direction", "individual")
if (!all(need_cols %in% names(samples)))
  stop("samples.csv lacks columns: ", paste(setdiff(need_cols, names(samples)), collapse = ", "), call. = FALSE)
if (anyDuplicated(samples$sample)) stop("duplicated sample names in samples.csv", call. = FALSE)

EXPECTED <- c("contig", "position", "variantID", "refAllele", "altAllele", "refCount", "altCount",
              "totalCount", "lowMAPQDepth", "lowBaseQDepth", "rawDepth", "otherBases", "improperPairs")
read_counts <- function(s) {
  f <- file.path(COUNTS_DIR, paste0(s, ".table"))
  if (!file.exists(f)) return(list(data = NULL, problem = "table missing"))
  d <- tryCatch(read.delim(f, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(d)) return(list(data = NULL, problem = "table is empty or has a header but no data rows"))
  if (!all(EXPECTED %in% names(d))) return(list(data = NULL, problem = "table lacks the ASEReadCounter columns"))
  if (nrow(d) == 0) return(list(data = NULL, problem = "table has a header but no data rows"))
  d$sample <- s
  list(data = d, problem = NA_character_)
}
read_wasp <- function(s) {
  f <- file.path(COUNTS_DIR, paste0(s, ".wasp_stats.tsv"))
  if (!file.exists(f)) return(list(data = NULL, problem = "wasp_stats.tsv missing"))
  d <- tryCatch(read.delim(f, stringsAsFactors = FALSE, colClasses = c("character", "character", "numeric")),
                error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0 || !identical(names(d), c("vW", "vA", "n")))
    return(list(data = NULL, problem = "wasp_stats.tsv empty or malformed"))
  d$sample <- s
  list(data = d, problem = NA_character_)
}
res <- lapply(samples$sample, read_counts); names(res) <- samples$sample
problems <- vapply(res, function(r) r$problem, character(1))
if (MODE == "outbred") {
  wres <- lapply(samples$sample, read_wasp); names(wres) <- samples$sample
  wprob <- vapply(wres, function(r) r$problem, character(1))
  problems <- ifelse(is.na(problems), wprob, ifelse(is.na(wprob), problems, paste(problems, wprob, sep = "; ")))
}
problems <- problems[!is.na(problems)]
if (length(problems) > 0)
  stop(length(problems), " sample(s) have no usable ASE count table in ", COUNTS_DIR, ":\n",
       paste0("  - ", names(problems), ": ", problems, collapse = "\n"),
       "\nRe-run the per-sample array script for these samples; tables of samples not listed in samples.csv are ignored.",
       call. = FALSE)
extra <- setdiff(sub("\\.table$", "", list.files(COUNTS_DIR, pattern = "\\.table$")), samples$sample)
if (length(extra) > 0) message("Ignoring tables of samples not in samples.csv: ", paste(extra, collapse = ", "))
sites_raw <- dplyr::bind_rows(lapply(res, function(r) r$data))
knitr::kable(dplyr::count(sites_raw, sample, name = "sites_in_table"), caption = "Sites per sample in the ASEReadCounter tables")
```

## Site filters

Sites are kept when total depth is at least `MIN_DEPTH`, when low-MAPQ + low-base-quality + other-base depth is at most 10 percent of the raw depth, and when no improper pairs overlap the site. Each filter is applied to the sites left by the previous one; the table shows how many sites each removed.

```{r filters}
s0 <- sites_raw
s1 <- dplyr::filter(s0, totalCount >= MIN_DEPTH)
s2 <- dplyr::filter(s1, (lowMAPQDepth + lowBaseQDepth + otherBases) <= 0.1 * rawDepth)
s3 <- dplyr::filter(s2, improperPairs == 0)
cnt <- function(d) as.integer(table(factor(d$sample, levels = samples$sample)))
filter_log <- data.frame(sample = samples$sample, sites_in = cnt(s0),
  removed_low_depth = cnt(s0) - cnt(s1), removed_low_quality_or_other_bases = cnt(s1) - cnt(s2),
  removed_improper_pairs = cnt(s2) - cnt(s3), sites_kept = cnt(s3))
knitr::kable(filter_log, caption = paste0("Sites removed by each filter (MIN_DEPTH = ", MIN_DEPTH, ")"))
if (any(filter_log$sites_kept == 0))
  stop("no sites left after filtering for: ", paste(filter_log$sample[filter_log$sites_kept == 0], collapse = ", "), call. = FALSE)
n_mismatch <- sum(s3$refCount + s3$altCount != s3$totalCount)
if (n_mismatch > 0) message(n_mismatch, " kept sites have refCount + altCount different from totalCount; ref + alt is used")
```

## Allele mapping

F1 mode: the REF count is strain A and the ALT count is strain B (the het-sites VCF was built with REF = strain A). Outbred mode: REF and ALT. All later tables use `ref_n`, `alt_n`, `total = ref_n + alt_n` and `alt_frac`.

```{r alleles}
sites <- s3 %>%
  dplyr::transmute(sample, contig, position,
                   SNP = paste0(contig, ":", position, ":", refAllele, ">", altAllele),
                   refAllele, altAllele, ref_n = refCount, alt_n = altCount,
                   total = refCount + altCount, ref_frac = refCount / (refCount + altCount),
                   alt_frac = altCount / (refCount + altCount)) %>%
  dplyr::left_join(samples[, c("sample", "condition", "cross_direction", "individual")], by = "sample")
cat("REF count =", ref_label, "; ALT count =", alt_label, "\n")
```

## Coverage and site counts

```{r coverage}
coverage <- sites %>% dplyr::group_by(sample) %>%
  dplyr::summarise(sites = dplyr::n(), total_reads = sum(total), median_depth = median(total),
                   min_depth = min(total), max_depth = max(total), .groups = "drop")
knitr::kable(coverage, caption = "Filtered sites and coverage per sample")
ggplot(sites, aes(total)) + geom_histogram(bins = 40) + scale_x_log10() + facet_wrap(~sample) +
  labs(x = "total reads at site (log10)", y = "sites") + theme_bw()
```

## Reference-bias diagnostic

For every sample the mean REF fraction is REF reads divided by REF + ALT reads summed over the filtered sites. A sample is flagged when its deviation from 0.5 exceeds `max(BIAS_TOL, 3 * SE)` with `SE = sqrt(0.25 / total reads)`. **This is a screen, not a test.** On data with planted or real allelic imbalance the all-sites mean is not expected to be 0.5, so a flag says "look at the ratios with caution", not "the sample is wrong"; a benign spread of about +0.02 to -0.03 was seen on synthetic data with no imbalance at about 2000 reads per sample.

```{r bias}
bias <- sites %>% dplyr::group_by(sample) %>%
  dplyr::summarise(n_sites = dplyr::n(), total_reads = sum(total), ref_reads = sum(ref_n), .groups = "drop") %>%
  dplyr::mutate(mean_ref_frac = ref_reads / total_reads, deviation = mean_ref_frac - 0.5,
                SE = sqrt(0.25 / total_reads), threshold = pmax(BIAS_TOL, 3 * SE),
                flagged = abs(deviation) > threshold)
knitr::kable(bias, digits = 4, caption = paste0("Reference-bias diagnostic (REF = ", ref_label, ")"))
if (any(bias$flagged)) cat("FLAGGED samples:", paste(bias$sample[bias$flagged], collapse = ", "),
                           "- show this flag next to every ratio table and figure.\n") else cat("No sample is flagged.\n")
ggplot(sites, aes(ref_frac)) + geom_histogram(bins = 25) + geom_vline(xintercept = 0.5, linetype = 2) +
  geom_vline(data = bias, aes(xintercept = mean_ref_frac, colour = flagged)) +
  scale_colour_manual(values = c(`FALSE` = "steelblue", `TRUE` = "firebrick")) +
  facet_wrap(~sample) + labs(x = paste0("REF (", ref_label, ") fraction per site"), y = "sites") + theme_bw()
```

## WASP removal (outbred only)

STAR's `vW` tag: 1 passed; 2 multi-mapping read; 3 variant base N; 4 remapping failed; 5 remapped read multi-maps; 6 remapped read maps to a different locus; 7 too many variants; `none` = the alignment overlaps no variant. `vA` is the first allele class of the alignment: 1 REF, 2 ALT, 3 no match, `-` = no tag. Alignments with `vW` 2-7 are removed, so the tables show how many were removed and from which allele class. WASP does not promise less bias; the REF fraction before and after is reported as measured.

```{r wasp, eval = (MODE == "outbred")}
wasp <- dplyr::bind_rows(lapply(wres, function(r) r$data))
knitr::kable(tidyr::pivot_wider(dplyr::count(wasp, sample, vW, wt = n, name = "alignments"),
                                names_from = vW, values_from = alignments, values_fill = 0),
             caption = "Alignments by vW code")
knitr::kable(tidyr::pivot_wider(dplyr::count(dplyr::filter(wasp, vW != "none"), sample, vA, wt = n, name = "alignments"),
                                names_from = vA, values_from = alignments, values_fill = 0),
             caption = "Tagged alignments by allele class (vA: 1 REF, 2 ALT, 3 no match)")
wasp_ref <- wasp %>% dplyr::filter(vA %in% c("1", "2")) %>% dplyr::group_by(sample) %>%
  dplyr::summarise(tagged_before = sum(n), ref_frac_before = sum(n[vA == "1"]) / sum(n),
                   tagged_after_vW1 = sum(n[vW == "1"]), ref_frac_after_vW1 = sum(n[vW == "1" & vA == "1"]) / sum(n[vW == "1"]),
                   removed_ref = sum(n[vW != "1" & vA == "1"]), removed_alt = sum(n[vW != "1" & vA == "2"]), .groups = "drop")
knitr::kable(wasp_ref, digits = 4, caption = "REF fraction among tagged alignments before and after keeping vW = 1")
```

## Checkpoint

```{r checkpoint}
constants <- list(MODE = MODE, STRAIN_A = STRAIN_A, STRAIN_B = STRAIN_B, MIN_DEPTH = MIN_DEPTH,
                  FDR_SIG = FDR_SIG, ABS_DEV_SIG = ABS_DEV_SIG, BIAS_TOL = BIAS_TOL, GTF_PATH = GTF_PATH)
saveRDS(list(sites = sites, samples = samples, constants = constants, bias = bias, filter_log = filter_log,
             wasp = if (MODE == "outbred") wasp else NULL),
        file.path(RESULTS_DIR, "ase_checkpoint.rds"))
cat("wrote", file.path(RESULTS_DIR, "ase_checkpoint.rds"), "\n")
sessionInfo()
```
````

Render it as an `sbatch -p bcc` job from the container (`run_01_import_qc.sh`, `-n 2 --mem=8G -t 0:30:00`, output under `{RESULTS_DIR}/logs`):

```bash
module add singularity/3.10.4 || { echo "ERROR: cannot load singularity" >&2; exit 1; }
singularity exec --bind {CWD} {R_SIF} Rscript -e 'rmarkdown::render("{CWD}/{TODAY_YYMMDD}_{WD_NAME}_01_import_qc.Rmd", output_dir = "{RESULTS_DIR}")' || { echo "ERROR: Rmd 01 failed" >&2; exit 1; }
```

If the job stops with "sample(s) have no usable ASE count table", show the listed samples to the user; do not delete anything, and do not work around it by globbing.

---

## Step 13 — Rmd 02: per-sample allelic imbalance

Write `{CWD}/{TODAY_YYMMDD}_{WD_NAME}_02_imbalance.Rmd` with the same header conventions as Step 12. It loads `ase_checkpoint.rds`, defines its own constants block (same names, same values; the Rmd stops if they differ from the checkpoint, so Rmd 01 and Rmd 02 can never disagree), pastes the statistics functions of Step 14 (the code between the two marker lines, without the marker lines) into the chunk marked below, and writes `{TODAY_YYMMDD}_{WD_NAME}_ASE_imbalance.xlsx` (sheets `SNP`, `Gene`, `Summary`) and `ase_imbalance_checkpoint.rds` to `{RESULTS_DIR}`.

- Per sample, `rho_h0 <- bb_estimate_rho(alt, total)` is estimated under H0 (p = 0.5) from all filtered sites (conservative when true imbalance exists, see Step 14). **Outbred** uses `rho_h0`. **F1** uses `bb_estimate_rho_gene(alt, total, gene)`, estimated with a free mean per gene and bias-corrected (`rho_corrected`), because real imbalance inflates `rho_h0`; the gene-level tests use `rho_gene = max(rho_corrected, RHO_MIN)`; if it is `NA` (fewer than 5 genes with 2 or more SNPs) the Rmd prints a WARNING, falls back to `rho_h0` and records that in the Summary (`rho_source`). The Summary shows all the estimates, and any estimator boundary warning per sample. With fewer than 20 sites `rho_h0` is `NA` and the sample is not tested.
- Per SNP, `bb_pvalue`, then Benjamini-Hochberg within each sample. The single column `sig` is `padj < FDR_SIG` and `|ALT fraction - 0.5| >= ABS_DEV_SIG` (F1: ALT is strain B). The Summary table and every plot use this column and nothing else; the Rmd checks that the counts drawn in the figures equal the Summary counts.
- Gene level, both modes: SNP positions are overlapped with the GTF exons by `GenomicRanges::findOverlaps` (`gene_id` from the GTF). **F1:** counts are not summed before testing (summing and applying the per-SNP `rho` to the total inflates the variance by about `1 + (n - 1) * rho` and destroys power); each gene's SNPs share one strain-B fraction and `bb_gene_lrt` tests it against 0.5 with `rho_gene`; the table reports `phat` (strain-B fraction), p, BH within the sample, and `sig` from `FDR_SIG` and `ABS_DEV_SIG` on `|phat - 0.5|`. **Outbred:** SNPs cannot be pooled without phasing, so the gene p-value is the `acat` combination of the SNP p-values, labelled "unphased, no direction" (no direction column); the gene is `sig` when its BH-adjusted `acat` p-value is below `FDR_SIG` and at least one of its SNPs deviates by `ABS_DEV_SIG` or more.
- The reference-bias flag from Rmd 01 is printed with the tables and written into the Summary sheet; if a sample is flagged, say so next to its ratios.

````rmd
---
title: "{PROJECT_TITLE} - ASE per-sample imbalance"
author: "{AUTHOR}"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(cache = FALSE, echo = TRUE, message = FALSE, warning = FALSE, fig.width = 10, fig.height = 7)
options(scipen = 9)
library(GenomicRanges); library(rtracklayer); library(openxlsx)   # Bioconductor first
library(tidyverse)
```

## Constants and checkpoint

```{r constants}
MODE        <- "{MODE}"
STRAIN_A    <- "{STRAIN_A}"
STRAIN_B    <- "{STRAIN_B}"
RESULTS_DIR <- "{RESULTS_DIR}"
GTF_PATH    <- "{GTF_PATH}"
DATE_TAG    <- "{TODAY_YYMMDD}_{WD_NAME}"
MIN_DEPTH   <- {MIN_DEPTH}
FDR_SIG     <- {FDR_SIG}
ABS_DEV_SIG <- {ABS_DEV_SIG}
RHO_MIN     <- {RHO_MIN}             # F1 gene level: floor for the bias-corrected overdispersion
BIAS_TOL    <- {BIAS_TOL}
ref_label <- STRAIN_A; alt_label <- STRAIN_B
ck <- readRDS(file.path(RESULTS_DIR, "ase_checkpoint.rds"))
same <- c(MODE = identical(ck$constants$MODE, MODE),
          vapply(c("MIN_DEPTH", "FDR_SIG", "ABS_DEV_SIG", "BIAS_TOL"),
                 function(k) isTRUE(all.equal(ck$constants[[k]], get(k))), logical(1)))
if (!all(same)) stop("constants differ from the Rmd 01 checkpoint (", paste(names(same)[!same], collapse = ", "),
                     "); re-render Rmd 01 with the same values", call. = FALSE)
sites <- ck$sites; bias <- ck$bias
knitr::kable(bias, digits = 4, caption = paste0("Reference-bias diagnostic from Rmd 01 (REF = ", ref_label, "); a screen, not a test"))
bias_note <- ifelse(bias$flagged, "FLAGGED: reference bias above tolerance, read ratios with caution", "not flagged")
names(bias_note) <- bias$sample
if (any(bias$flagged)) cat("FLAGGED samples:", paste(bias$sample[bias$flagged], collapse = ", "), "\n")
```

## Statistics functions

```{r stats}
# <<< paste here the code of Step 14 between the two marker lines (the marker lines themselves are not pasted) >>>
```

## SNP to gene map

```{r genemap}
gtf <- rtracklayer::import(GTF_PATH)
if (is.null(gtf$gene_id)) stop("the GTF has no gene_id attribute", call. = FALSE)
ex <- gtf[gtf$type == "exon" & !is.na(gtf$gene_id)]
pos <- unique(sites[, c("contig", "position")])
gr <- GenomicRanges::GRanges(pos$contig, IRanges::IRanges(pos$position, width = 1))
gr_chr <- unique(as.character(GenomicRanges::seqnames(gr))); ex_chr <- unique(as.character(GenomicRanges::seqnames(ex)))
if (length(intersect(gr_chr, ex_chr)) == 0)
  stop("no contig names in common between the count tables (", paste(head(gr_chr, 3), collapse = ", "),
       ") and the GTF (", paste(head(ex_chr, 3), collapse = ", "), ")", call. = FALSE)
hits <- GenomicRanges::findOverlaps(gr, ex)
snp_gene <- data.frame(contig = pos$contig[S4Vectors::queryHits(hits)], position = pos$position[S4Vectors::queryHits(hits)],
                       gene_id = ex$gene_id[S4Vectors::subjectHits(hits)], stringsAsFactors = FALSE) %>% dplyr::distinct()
cat(nrow(pos), "SNP positions,", length(unique(snp_gene$gene_id)), "genes with at least one SNP\n")
```

## Per-SNP beta-binomial test

The overdispersion `rho` is estimated per sample. **Outbred:** under H0 (p = 0.5) from all filtered sites (`bb_estimate_rho`); when real imbalance exists this inflates `rho`, which lowers power and never inflates false positives (conservative). **F1:** the H0-based value (`rho_h0`) is inflated by every truly imbalanced gene, so `rho` is instead estimated with a free mean per gene (`bb_estimate_rho_gene`, genes with at least 2 SNPs). That estimate (`rho_free`) is biased low, because a free mean per gene absorbs part of the variance (by about (k - 1) / k for genes with k SNPs), so it is bias-corrected (`rho_corrected`, see Step 14). The SNP-level tests use `rho_corrected`; the gene-level tests use `rho_gene = max(rho_corrected, RHO_MIN)`. If the free-mean estimate is not available (fewer than 5 usable genes) the H0-based value is used and a WARNING is printed and written to the Summary (`rho_source`). Any estimator warning (for example an estimate at the upper boundary) is printed per sample and written to `rho_warning`.

```{r snp}
bb_p_safe <- function(x, n, rho) if (is.na(rho)) NA_real_ else bb_pvalue(x, n, rho = rho)
collect_warnings <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") })
  list(value = v, warn = paste(unique(w), collapse = "; "))
}
snp <- dplyr::bind_rows(lapply(split(sites, sites$sample), function(d) {
  r0 <- collect_warnings(bb_estimate_rho(d$alt_n, d$total)); rho_h0 <- r0$value; warn <- r0$warn
  rho_free <- NA_real_; rho_corrected <- NA_real_; rho <- rho_h0; src <- "H0-based"
  if (MODE == "f1") {
    g <- snp_gene$gene_id[match(paste(d$contig, d$position), paste(snp_gene$contig, snp_gene$position))]
    r1 <- collect_warnings(bb_estimate_rho_gene(d$alt_n, d$total, g))
    warn <- paste(c(warn, r1$warn)[nzchar(c(warn, r1$warn))], collapse = "; ")
    if (!is.na(r1$value[["corrected"]])) {
      rho_free <- r1$value[["free"]]; rho_corrected <- r1$value[["corrected"]]; rho <- rho_corrected
      src <- "free-mean per gene, bias-corrected"
    } else {
      cat("WARNING: sample", d$sample[1], "has fewer than 5 genes with 2 or more SNPs; using the H0-based rho\n"); src <- "H0-based (fallback)" }
  }
  if (nzchar(warn)) cat("WARNING (rho estimation), sample", d$sample[1], ":", warn, "\n")
  d$rho_h0 <- rho_h0; d$rho_free <- rho_free; d$rho_corrected <- rho_corrected; d$rho <- rho
  d$rho_gene <- if (MODE == "f1") max(rho, RHO_MIN) else NA_real_   # floor for the gene-level tests (F1)
  d$rho_source <- src; d$rho_warning <- warn
  d$p <- mapply(bb_p_safe, d$alt_n, d$total, MoreArgs = list(rho = rho))
  d$padj <- p.adjust(d$p, method = "BH")
  d
})) %>%
  dplyr::mutate(dev = abs(alt_frac - 0.5),
                sig = !is.na(padj) & padj < FDR_SIG & dev >= ABS_DEV_SIG,
                direction = dplyr::case_when(!sig ~ "none", alt_frac > 0.5 ~ paste(alt_label, "higher"),
                                             TRUE ~ paste(ref_label, "higher")))
knitr::kable(dplyr::distinct(snp, sample, rho_h0, rho_free, rho_corrected, rho, rho_gene, rho_source), digits = 4,
             caption = "Overdispersion per sample: rho_h0 (p = 0.5 at every site, inflated by real imbalance), rho_free (free mean per gene, biased low), rho_corrected (bias-corrected, used for SNP tests), rho_gene = max(rho_corrected, RHO_MIN) (used for gene tests); F1 only for the last three. NA = fewer than 20 sites, sample not tested")
```

## Gene level

**F1: counts are not summed before testing.** Summing the two strains' reads over a gene's SNPs and applying the per-SNP overdispersion to the total would inflate the variance by `1 + (n - 1) * rho` on the summed depth `n` (about 45 at n = 450), which destroys power. Instead each gene's SNPs share one strain-B fraction `p`, and `bb_gene_lrt` runs a likelihood-ratio test of `p = 0.5` against `p` free on the SNP-level counts with the sample's `rho`; `phat` is the strain-B fraction estimate (`alt_frac` in the table). The summed reads are shown for information only. **Outbred:** `acat` over the SNP p-values.

```{r gene}
snp_by_gene <- dplyr::inner_join(snp, snp_gene, by = c("contig", "position"))

if (MODE == "f1") {
  gene <- dplyr::bind_rows(lapply(split(snp_by_gene, snp_by_gene$sample), function(d) {
    dplyr::bind_rows(lapply(split(d, d$gene_id), function(g) {
      r <- if (is.na(g$rho_gene[1])) list(p = NA_real_, phat = NA_real_) else bb_gene_lrt(g$alt_n, g$total, g$rho_gene[1])
      data.frame(sample = g$sample[1], gene_id = g$gene_id[1], n_snps = nrow(g), ref_n = sum(g$ref_n), alt_n = sum(g$alt_n),
                 total = sum(g$total), rho = g$rho_gene[1], alt_frac = r$phat, p = r$p, stringsAsFactors = FALSE)
    }))
  })) %>%
    dplyr::mutate(dev = abs(alt_frac - 0.5)) %>%
    dplyr::group_by(sample) %>% dplyr::mutate(padj = p.adjust(p, method = "BH")) %>% dplyr::ungroup() %>%
    dplyr::mutate(sig = !is.na(padj) & padj < FDR_SIG & dev >= ABS_DEV_SIG,
                  direction = dplyr::case_when(!sig ~ "none", alt_frac > 0.5 ~ paste(alt_label, "higher"),
                                               TRUE ~ paste(ref_label, "higher")))
} else {
  gene <- snp_by_gene %>% dplyr::filter(!is.na(p)) %>% dplyr::group_by(sample, gene_id) %>%
    dplyr::summarise(n_snps = dplyr::n(), acat_p = acat(p), max_dev = max(dev), .groups = "drop") %>%
    dplyr::group_by(sample) %>% dplyr::mutate(padj = p.adjust(acat_p, method = "BH")) %>% dplyr::ungroup() %>%
    dplyr::mutate(sig = !is.na(padj) & padj < FDR_SIG & max_dev >= ABS_DEV_SIG,
                  note = "unphased, no direction")
}
gene <- dplyr::arrange(gene, sample, gene_id)
snp_annot <- snp_gene %>% dplyr::group_by(contig, position) %>%
  dplyr::summarise(Gene = paste(unique(gene_id), collapse = ";"), .groups = "drop")
snp <- dplyr::left_join(snp, snp_annot, by = c("contig", "position"))
knitr::kable(head(dplyr::filter(gene, sig), 30), digits = 4, caption = "Significant genes (first 30)")
```

## Summary

```{r summary}
snp_plot <- dplyr::filter(snp, !is.na(p))      # data behind the SNP figure
gene_plot <- dplyr::filter(gene, if (MODE == "f1") !is.na(p) else !is.na(acat_p))   # data behind the gene figure
summary_tbl <- snp %>% dplyr::group_by(sample) %>%
  dplyr::summarise(rho_h0 = dplyr::first(rho_h0), rho_free = dplyr::first(rho_free), rho_corrected = dplyr::first(rho_corrected), rho_used = dplyr::first(rho), rho_gene = dplyr::first(rho_gene), rho_source = dplyr::first(rho_source), rho_warning = dplyr::first(rho_warning), sites_tested = sum(!is.na(p)), sig_sites = sum(sig),
                   sig_sites_alt_higher = sum(sig & alt_frac > 0.5), sig_sites_ref_higher = sum(sig & alt_frac < 0.5),
                   .groups = "drop")
gsum <- if (MODE == "f1") {
  gene %>% dplyr::group_by(sample) %>%
    dplyr::summarise(genes_tested = sum(!is.na(p)), sig_genes = sum(sig),
                     sig_genes_alt_higher = sum(sig & alt_frac > 0.5), sig_genes_ref_higher = sum(sig & alt_frac < 0.5), .groups = "drop")
} else {
  gene %>% dplyr::group_by(sample) %>%
    dplyr::summarise(genes_tested = sum(!is.na(acat_p)), sig_genes = sum(sig), .groups = "drop")
}
summary_tbl <- summary_tbl %>% dplyr::left_join(gsum, by = "sample") %>%
  dplyr::left_join(dplyr::select(bias, sample, mean_ref_frac, ref_bias_flagged = flagged), by = "sample") %>%
  dplyr::mutate(ref_bias_note = bias_note[sample],
                gene_level = if (MODE == "f1") "per-gene LRT on SNP-level counts (shared strain fraction)" else "unphased, no direction (acat of SNP p-values)")
knitr::kable(summary_tbl, digits = 4, caption = "Summary per sample (FDR_SIG, ABS_DEV_SIG as in the constants block)")
# the figures are drawn from snp_plot / gene_plot: their sig counts must equal the Summary counts
fig_sites <- tapply(snp_plot$sig, factor(snp_plot$sample, levels = summary_tbl$sample), sum)
fig_genes <- tapply(gene_plot$sig, factor(gene_plot$sample, levels = summary_tbl$sample), sum)
fig_sites[is.na(fig_sites)] <- 0; fig_genes[is.na(fig_genes)] <- 0
stopifnot(identical(as.integer(fig_sites), as.integer(summary_tbl$sig_sites)),
          identical(as.integer(fig_genes), as.integer(summary_tbl$sig_genes)))
cat("Figure sig counts equal Summary counts: TRUE\n")
```

## Figures

```{r figures}
cols <- c(`FALSE` = "grey60", `TRUE` = "firebrick")
lab_sites <- snp_plot %>% dplyr::group_by(sample) %>%
  dplyr::summarise(lab = paste0(dplyr::first(sample), ": ", sum(sig), " significant of ", dplyr::n(), " sites"), .groups = "drop")
p1 <- dplyr::left_join(snp_plot, lab_sites, by = "sample") %>%
  ggplot(aes(total, alt_frac, colour = sig)) + geom_hline(yintercept = 0.5, linetype = 2) + geom_point(alpha = 0.7) +
  scale_x_log10() + scale_colour_manual(values = cols) + facet_wrap(~lab) +
  labs(x = "total reads at site (log10)", y = paste0(alt_label, " fraction"), colour = "sig",
       caption = paste0("sig: BH < ", FDR_SIG, " and |fraction - 0.5| >= ", ABS_DEV_SIG,
                        if (any(bias$flagged)) "; some samples flagged for reference bias" else "")) + theme_bw()
print(p1)
ggsave(file.path(RESULTS_DIR, paste0(DATE_TAG, "_ASE_sites_vs_depth.pdf")), p1, width = 10, height = 7)
if (MODE == "f1") {
  lab_genes <- gene_plot %>% dplyr::group_by(sample) %>%
    dplyr::summarise(lab = paste0(dplyr::first(sample), ": ", sum(sig), " significant of ", dplyr::n(), " genes"), .groups = "drop")
  p2 <- dplyr::left_join(gene_plot, lab_genes, by = "sample") %>%
    ggplot(aes(gene_id, alt_frac, colour = sig)) + geom_hline(yintercept = 0.5, linetype = 2) + geom_point(size = 2) +
    scale_colour_manual(values = cols) + facet_wrap(~lab) + coord_flip() +
    labs(x = NULL, y = paste0(alt_label, " fraction (LRT estimate)"), colour = "sig") + theme_bw()
} else {
  lab_genes <- gene_plot %>% dplyr::group_by(sample) %>%
    dplyr::summarise(lab = paste0(dplyr::first(sample), ": ", sum(sig), " significant of ", dplyr::n(), " genes"), .groups = "drop")
  p2 <- dplyr::left_join(gene_plot, lab_genes, by = "sample") %>%
    ggplot(aes(gene_id, -log10(padj), colour = sig)) + geom_point(size = 2) + scale_colour_manual(values = cols) +
    facet_wrap(~lab) + coord_flip() +
    labs(x = NULL, y = "-log10 adjusted acat p (unphased, no direction)", colour = "sig") + theme_bw()
}
print(p2)
ggsave(file.path(RESULTS_DIR, paste0(DATE_TAG, "_ASE_genes.pdf")), p2, width = 10, height = 7)
```

## Export

```{r export}
rename_out <- function(d) {
  names(d) <- sub("^ref_n$", paste0("reads_", ref_label), names(d))
  names(d) <- sub("^alt_n$", paste0("reads_", alt_label), names(d))
  names(d) <- sub("alt_higher", paste0(alt_label, "_higher"), names(d))
  names(d) <- sub("ref_higher", paste0(ref_label, "_higher"), names(d))
  names(d) <- sub("^alt_frac$", paste0("frac_", alt_label), names(d))
  names(d) <- sub("^mean_ref_frac$", paste0("mean_frac_", ref_label), names(d))
  d
}
xlsx_file <- file.path(RESULTS_DIR, paste0(DATE_TAG, "_ASE_imbalance.xlsx"))
wb <- openxlsx::createWorkbook()
for (nm in c("SNP", "Gene", "Summary")) openxlsx::addWorksheet(wb, nm)
openxlsx::writeData(wb, "SNP", rename_out(dplyr::select(snp, SNP, Gene, sample, condition, contig, position, ref_n, alt_n, total, alt_frac, rho, p, padj, dev, sig, direction)))
openxlsx::writeData(wb, "Gene", rename_out(gene))
openxlsx::writeData(wb, "Summary", rename_out(summary_tbl))
openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)
saveRDS(list(snp = snp, gene = gene, summary = summary_tbl, constants = ck$constants),
        file.path(RESULTS_DIR, "ase_imbalance_checkpoint.rds"))
cat("wrote", xlsx_file, "\n")
sessionInfo()
```
````

Render it with `run_02_imbalance.sh` (same pattern as `run_01_import_qc.sh`, Rmd 02 file name, `-n 2 --mem=16G -t 1:00:00`, submitted with `sbatch -p bcc --dependency=afterok:<Rmd 01 job>`).

---

## Step 14 — Statistics functions

These base-R functions are pasted verbatim into Rmd 02 (`{TODAY_YYMMDD}_{WD_NAME}_02_imbalance.Rmd`); the marker lines delimit the code that `ase-pipeline/tests/r/test_ase_stats.R` extracts from this file and unit-tests, so the tested code is the shipped code. They need no packages beyond base R.

```r
# --- ase-stats-begin
bb_negloglik <- function(par, x, n, p0 = 0.5) {
  rho <- plogis(par)
  if (rho < 1e-8) return(-sum(dbinom(x, n, p0, log = TRUE)))
  a <- p0 * (1 - rho) / rho; b <- (1 - p0) * (1 - rho) / rho
  -sum(lchoose(n, x) + lbeta(x + a, n - x + b) - lbeta(a, b))
}
bb_estimate_rho <- function(x, n, p0 = 0.5) {
  keep <- n > 0; x <- x[keep]; n <- n[keep]
  if (length(x) < 20) return(NA_real_)
  fit <- optimize(bb_negloglik, interval = c(-12, 2), x = x, n = n, p0 = p0)
  if (fit$minimum < -12 + 0.05) return(0)   # lower boundary: no overdispersion, pure binomial
  if (fit$minimum > 2 - 0.05) warning("rho estimate at the upper boundary (about 0.88): extreme overdispersion, check the data")
  rho <- plogis(fit$minimum)
  if (rho < 1e-4) 0 else rho
}
bb_pvalue <- function(x, n, rho, p0 = 0.5) {
  if (is.na(n) || is.na(x) || n <= 0 || x < 0 || x > n || x != round(x)) return(NA_real_)
  k <- 0:n
  lp <- if (rho < 1e-8) dbinom(k, n, p0, log = TRUE) else {
    a <- p0 * (1 - rho) / rho; b <- (1 - p0) * (1 - rho) / rho
    lchoose(n, k) + lbeta(k + a, n - k + b) - lbeta(a, b)
  }
  min(1, sum(exp(lp[lp <= lp[x + 1] + 1e-9])))
}
acat <- function(p) {   # Cauchy combination with equal weights; p is capped at 0.99 because tan((0.5 - 1) * pi) is about -1.6e16 and would swamp any signal
  p <- p[!is.na(p)]; if (length(p) == 0) return(NA_real_)
  p <- pmin(pmax(p, 1e-15), 0.99)
  0.5 - atan(mean(tan((0.5 - p) * pi))) / pi
}
bb_ll <- function(x, n, p, rho) {   # summed beta-binomial log-likelihood over vectors; rho < 1e-8 -> binomial
  if (rho < 1e-8) return(sum(dbinom(x, n, p, log = TRUE)))
  a <- p * (1 - rho) / rho; b <- (1 - p) * (1 - rho) / rho
  sum(lchoose(n, x) + lbeta(x + a, n - x + b) - lbeta(a, b))
}
bb_gene_lrt <- function(x, n, rho) {   # H0: p = 0.5 vs H1: p free; the gene's SNPs share ONE p (counts are NOT summed)
  keep <- n > 0; x <- x[keep]; n <- n[keep]
  if (length(x) == 0) return(list(p = NA_real_, phat = NA_real_))
  ll0 <- bb_ll(x, n, 0.5, rho)
  fit <- optimize(function(p) -bb_ll(x, n, p, rho), interval = c(1e-6, 1 - 1e-6))
  stat <- max(0, 2 * (-fit$objective - ll0))
  list(p = pchisq(stat, 1, lower.tail = FALSE), phat = fit$minimum)
}
bb_estimate_rho_gene <- function(x, n, gene) {   # c(free, corrected): rho with a FREE mean per gene (no inflation from true imbalance)
  keep <- n > 0 & !is.na(gene); x <- x[keep]; n <- n[keep]; gene <- gene[keep]
  idx <- split(seq_along(x), gene)
  idx <- idx[vapply(idx, length, integer(1)) >= 2]
  if (length(idx) < 5) return(c(free = NA_real_, corrected = NA_real_))
  nll <- function(par) {
    rho <- plogis(par)
    sum(vapply(idx, function(i) optimize(function(p) -bb_ll(x[i], n[i], p, rho), interval = c(1e-6, 1 - 1e-6))$objective, numeric(1)))
  }
  fit <- optimize(nll, interval = c(-12, 2))
  if (fit$minimum < -12 + 0.05) return(c(free = 0, corrected = 0))
  if (fit$minimum > 2 - 0.05) warning("rho estimate at the upper boundary (about 0.88): extreme overdispersion, check the data")
  free <- plogis(fit$minimum)
  if (free < 1e-4) return(c(free = 0, corrected = 0))
  # profile ML with one free mean per gene shrinks the whole variance factor 1 + (n - 1) * rho by (k - 1) / k, so undo that and solve for rho
  k <- vapply(idx, length, integer(1)); cf <- sum(k) / sum(k - 1); m <- mean(unlist(lapply(idx, function(i) n[i])) - 1)
  c(free = free, corrected = min(cf * free + (cf - 1) / m, 0.99))
}
# --- ase-stats-end
```

- `bb_pvalue` is the two-sided exact beta-binomial test against `p0`: the p-value is the total probability of all outcomes that are no more likely than the observed count. With `rho = 0` it is the exact binomial test.
- **`rho` (overdispersion) is estimated per sample under H0 (p = 0.5) from all filtered sites.** When real imbalance exists, the extra variance is absorbed into `rho`, which inflates it. An inflated `rho` lowers power but never inflates false positives, so the estimate is conservative. It is estimated by maximum likelihood on the logistic scale (`optimize`), and a value below 1e-6 is reported as 0 (pure binomial). Fewer than 20 sites gives `NA` and the sample is not tested.
- `acat` is the Cauchy combination test with equal weights; it is valid for correlated p-values, which is why it is used for the SNPs of one gene. It carries no direction. P-values are capped at 0.99 before the transform: `bb_pvalue` returns exactly 1 at the modal count, and `tan((0.5 - 1) * pi)` is about -1.6e16, which would swamp any signal and make the combined p-value 1. Values are also clipped at 1e-15 from below.
- The unit tests check the null size and uniformity of the p-values, power, recovery of `rho` (including `rho` near 0), the edge cases (`n = 0`, `x = 0`, `x = n`) and `acat`. Run them from the repository root with `singularity exec --bind /net/bmc-lab3 <bulkrnaseq sif> Rscript ase-pipeline/tests/r/test_ase_stats.R ase-pipeline/ase-pipeline.md` inside an `sbatch -p bcc` job.
- **F1 gene level uses `bb_estimate_rho_gene` and `bb_gene_lrt`, not summed counts.** `bb_estimate_rho` (H0, p = 0.5 everywhere) treats the between-SNP spread caused by true imbalance as overdispersion, so with many imbalanced genes it is inflated (in the unit test about four times the true value), and applying it to counts summed over a gene multiplies the variance by about `1 + (n - 1) * rho`. `bb_estimate_rho_gene` maximises the likelihood with a free mean for every gene (genes with at least 2 SNPs; `NA` below 5 such genes) and returns `c(free, corrected)`. The free estimate is biased low: a free mean per gene shrinks the whole variance factor `1 + (n - 1) * rho` by `(k - 1) / k` (Neyman-Scott), which for 4-SNP genes gives 0.0115 for a true 0.02. The corrected value undoes that shrinkage, `rho_corrected = cf * rho_free + (cf - 1) / m` with `cf = sum(k) / sum(k - 1)` and `m` the mean of `n - 1` (0.021 in the same test). Rmd 02 then uses `rho_gene = max(rho_corrected, RHO_MIN)` for the gene tests (`RHO_MIN`, default 0.01, from Step 8), so a noisy or zero estimate can never make them anti-conservative. `bb_gene_lrt` tests `p = 0.5` against a shared free `p` per gene on the SNP-level counts and returns `phat`. Outbred mode keeps `bb_estimate_rho`, `bb_pvalue` and `acat`.
- Boundary handling: both `rho` estimators return exactly 0 when the optimum sits at the lower boundary (or below 1e-4), and emit a warning when it sits at the upper boundary; Rmd 02 prints such warnings per sample and writes them to the Summary (`rho_warning`). `bb_pvalue` returns `NA` for `x` below 0, above `n`, non-integer or missing.
