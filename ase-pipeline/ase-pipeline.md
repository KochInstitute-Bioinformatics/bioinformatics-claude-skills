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
- `{TODAY}` = today's date formatted as `YYYY-MM-DD` (for example `2026-09-29`; the same format as the `results/YYYY-MM-DD_*` data-safety convention), used for the results directory and as the prefix of the Rmd and report file names
- `{RESULTS_DIR}` = `{CWD}/results/{TODAY}_{WD_NAME}`, for example `results/2026-09-29_proj` (create it with `mkdir -p` when the first output is written; an existing results directory is reused and never deleted, see Step 15)

**Paths with spaces or commas are not supported.** The generated scripts put paths into `#SBATCH -o` lines, comma-separated `--bind` lists and unquoted shell words, so a space or a comma breaks them. If `{CWD}` contains a space or a comma, stop here and tell the user: "ERROR: the path '{CWD}' contains a space or a comma, which the generated SLURM and Singularity commands cannot handle. Move or link the project to a path without spaces or commas and start again." Apply the same check (and the same stop) to every path the user gives later: `{GENOME_DIR}`, `{FASTA_PATH}`, `{GTF_PATH}`, `{PARENTAL_VCF}`, each genotype VCF and each FASTQ path in `{SAMPLES_CSV}` (Step 5).

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
  Rscript -e 'for (p in c("openxlsx","tidyverse","GenomicRanges","rtracklayer")) cat(p, requireNamespace(p, quietly = TRUE), "\n")'
```

Mandatory packages: `openxlsx`, `tidyverse`, `GenomicRanges`, `rtracklayer` (the beta-binomial models are base R; `aod` is used only by the unit tests and `lme4` is not used). Submit it with `sbatch -p bcc check_r_packages.sh` (every job is submitted with `-p bcc`), wait for the job to finish, then read the result from the job output file `check_r_packages_<jobid>.out`. If any mandatory package prints `FALSE`, **stop** and ask the user to rebuild the `bulkrnaseq` image; do not substitute another model or package. (`VGAM`, `glmmTMB`, `MBASED` and `VariantAnnotation` are absent from the image and must never be required.)

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
- **F1:** ask once for the two strain names, `{STRAIN_A}` (the reference strain) and `{STRAIN_B}` (names from the VCF header, or the fallback of Step 6), and, per sample, `cross_direction`: `cross_direction` is exactly `AxB` (mother `{STRAIN_A}`, father `{STRAIN_B}`) or `BxA` (mother `{STRAIN_B}`). The letters are literally A and B, not strain names, because the parent-of-origin sign of Rmd 03 depends on it. Ask per group of samples (numbered: 1. `AxB` · 2. `BxA`), never as free text. Leave `individual` equal to `sample`: one sample per animal (each F1 animal is its own individual), because the reciprocal and differential tests treat every sample as an independent animal. If an animal has technical replicates (several lanes or libraries), merge the FASTQ files of the technical replicates of an animal before the wizard (for example `cat` the `_R1` files and the `_R2` files in the same order) and give the merged files as one sample; several tissues or time points of one animal are not supported in F1 mode.
- **Outbred:** `individual` per row (which person/animal the sample comes from; several samples may share one); leave `cross_direction` as `NA`.

Show the full table, ask "Does this look correct?", and write the file only after confirmation (if it exists, ask: overwrite or choose another filename).

**Validation rules:** no dashes or spaces in `sample`, `condition`, `cross_direction` or `individual`; no space or comma in any `fastq_1` / `fastq_2` path (stop with the Step 0 message naming the path); sample names unique; each condition should have replicates (defined below; warn, do not stop, when it has only one sample); in F1 mode `cross_direction` is exactly `AxB` or `BxA` for every sample (stop and ask again otherwise); in F1 mode every `individual` value must be unique: run `awk -F, 'NR > 1 {c[$6]++} END {for (i in c) if (c[i] > 1) print i}' {SAMPLES_CSV}` and, if it prints anything, stop and tell the user: "F1 mode needs one sample per animal; these individual values occur more than once: <the printed values>. Merge the FASTQ files of the technical replicates of each animal into one sample and start again." (outbred mode may repeat an `individual`); reciprocal analysis needs at least 2 samples of each cross direction (checked again in Step 7).

**Replicates (definition used here and in Step 7):** replicates = at least 2 samples in the condition. A condition with a single sample triggers the warning above.

### Step 5, part 2 — Read length and STAR index

**Read length (always detect it).** STAR's `sjdbOverhang` is fixed when the index is built. Define `{FASTQ_FILE}` as the first `fastq_1` of each of up to 5 distinct samples in `{SAMPLES_CSV}`. The `zcat | head` read below is bounded (1000 reads) and intentionally light, so it is acceptable on the login node:

```bash
zcat {FASTQ_FILE} | awk 'NR%4==2 {print length($0)}' | head -n 1000 | sort -n | uniq -c | sort -rn | head -3
```

The mode of the lengths is `{READ_LENGTH}`. If lengths differ between samples, report the distribution and warn that the index is tuned to the most common length. Set `{SJDB_OVERHANG}` = `{READ_LENGTH}` − 1 and tell the user: "Detected read length {READ_LENGTH} bp -> sjdbOverhang {SJDB_OVERHANG}."

**STAR index: one index per read length and per mode (F1: also per mask).** Set `{STAR_INDEX}` to:

- F1: the prefix `{GENOME_DIR}/index/star_ase_masked_{STRAIN_A}_{STRAIN_B}_sjdb{SJDB_OVERHANG}` (built from the third-allele masked genome, so it never collides with an unmasked index of the same read length). The masked genome depends on the parental SNP set, which exists only once the prep job (or the mouse helper) has run, so the prep job completes the name at run time: it computes `MASK_KEY`, the first 8 characters of the md5 of the masked site records (`CHROM POS REF ALT` of `masked_sites.vcf`, headers excluded so the same sites always give the same key), and uses `{STAR_INDEX}_<MASK_KEY>/`, for example `star_ase_masked_C57BL_6NJ_A_J_sjdb99_3f9a1c2e/`. A second strain pair or a corrected parental VCF therefore gets its own index and never reuses a stale mask. The prep job writes the full index path to `{RESULTS_DIR}/reference/star_index_path.txt`, and both array jobs (Step 11) read the index path from there;
- outbred: `{GENOME_DIR}/index/star_ase_sjdb{SJDB_OVERHANG}/` (built from the normal reference).

**Built once, reused only when the key matches.** The prep job reuses an existing index only if `SA`, `Genome` and `sjdbList.out.tab` are non-empty there and the file `index_key.txt` inside it, written as the last step of a successful build, holds the same key (F1: `MASK_KEY`; outbred: `unmasked`). Otherwise it builds the index at that path (STAR 2.7.10b from `{STAR_SIF}`, on a compute node), so an interrupted build is never reused either. The prep job computes the genome length and `genomeSAindexNbases` = min(14, max(4, floor(log2(L)/2 - 1))) with L = genome length (for example 300,000 bp gives 8; 40,001 gives 6; 3.1e9 gives 14; 2 kb gives 4) itself, in shell at run time, because the FASTA may not exist yet when the wizard writes the script; the wizard never substitutes that value.

---

## Step 6 — Genotype source

The allelic analysis needs heterozygous SNP positions.

**F1 mode.** Ask (numbered): "1. I have a parental-difference VCF · 2. Build it from the Mouse Genomes Project (mouse helper)".

**Reference strain first.** Ask: "Is one of the two parental strains the strain of the reference assembly (or a substrain of it, for example C57BL/6NJ on GRCm39, which is C57BL/6J)?" If not (for example A/J x CAST/EiJ), say: "F1 mode needs one parent to be the reference strain: the masked reference and the REF = strain A convention assume the assembly base is strain A's allele. A cross of two non-reference strains is not supported by this skill." and stop the F1 setup. Otherwise that parent is `{STRAIN_A}`.

1. `{PARENTAL_VCF}`: a biallelic SNP VCF where the REF allele is `{STRAIN_A}` (the reference strain) and the ALT allele is `{STRAIN_B}`. Either form works: a plain `.vcf` is accepted (the prep job only reads it with `bcftools view`, which needs no index), and a bgzipped `.vcf.gz` works too; no bgzip or tabix step is needed, and the input is never modified. The prep job checks the contig names against the FASTA on a compute node (step 1 of `prep_f1_reference.sh`) and stops before any alignment if they differ; nothing is checked on the login node.
   - **VCF with sample columns** (for example the two strains' own genotypes, or a Mouse Genomes Project extract made outside the helper): both `{STRAIN_A}` and `{STRAIN_B}` must be sample names in the header (the prep job stops otherwise). The prep job keeps **only** the sites where strain A is homozygous REF and strain B homozygous ALT (`0/0` and `1/1`, phased or not). It drops, and counts in its log, the sites where both strains have the same homozygous genotype (for example both `1/1`: both differ from the assembly but not from each other) and the sites where either strain is heterozygous or missing. Sites where strain A is `1/1` and strain B `0/0` (the strains differ, but strain A carries the non-reference allele) break the REF = strain A convention: by default (`A_ALT_SITES=stop`) the job stops and prints their number. If strain A is a substrain of the assembly strain and the user confirms that a few such sites are expected (for example C57BL_6NJ, which differs from the C57BL/6J assembly at some sites), the wizard sets `A_ALT_SITES=drop`, which drops and counts them like the mouse helper does. The job also stops if fewer than `MIN_PARENTAL_SITES` sites remain (default 1000; 20 for small test data).
   - **Sites-only VCF** (no sample columns): every biallelic SNP is taken on trust as REF = strain A, ALT = strain B. Tell the user: "A sites-only VCF cannot be checked for which strain carries which allele; you guarantee that REF is the {STRAIN_A} allele and ALT the {STRAIN_B} allele at every site, and that the two strains differ there." The prep job still checks that REF equals the FASTA base and that at least `MIN_PARENTAL_SITES` sites remain. A sites-only VCF is fine for counting, because the prep job writes its own genotyped het-sites VCF.
   **Strain names.** If the VCF header has two sample columns, offer those names for `{STRAIN_A}` and `{STRAIN_B}` (the reference strain is `{STRAIN_A}`). If it has no names (sites-only VCF), ask the user; when no names are given, fall back to the Mouse Genomes Project spellings of the B6 x AJ example, `{STRAIN_A}` = `C57BL_6NJ` and `{STRAIN_B}` = `A_J`, and say so. For a sites-only VCF the names only label tables and figures; REF is always strain A and ALT always strain B.
2. The mouse helper queries the Mouse Genomes Project (release REL-2112-v8, GRCm39, contigs `1`, `2`, ... without `chr`) for the two strains as a per-chromosome array job and writes `{PARENTAL_VCF}` with the same convention (it already keeps only strain A `0/0` with strain B `1/1` and writes a sites-only VCF). It is a compute job; do not run it on the login node. It needs neither the FASTA nor its `.fai`: the chromosome list is the fixed GRCm39 list of Step 10, and the contig names are compared with the FASTA in the prep job, which runs after it.

Either way the prep job also writes `f1_het_sites.vcf.gz`, a bgzipped and tabix-indexed het-sites VCF (single sample `F1`, genotype `0/1` at every parental SNP), beside the parental-difference VCF. ASEReadCounter needs this genotype column and an indexed VCF; a sites-only VCF gives zero counts.

**Outbred mode.** One VCF per `individual`, `{GENOTYPE_VCFS}` (a mapping individual -> path). Ask (numbered): "1. Matched WGS or array genotypes (preferred) · 2. The rnavar `variant_calling/` filtered VCFs".

- Option 2 requires this caveat, stated verbatim to the user: "RNA-derived genotypes are circular: heterozygous sites with strong imbalance may be called homozygous, so genotypes called from the same RNA-seq bias results toward balance and miss lowly expressed sites. External genotypes are preferred." Continue only after the user confirms.
- The prep job turns each individual's VCF into a single-sample VCF restricted to biallelic heterozygous SNPs (`bcftools view -s SAMPLE -g het -v snps -m2 -M2`, plus quality filters). The per-individual VCF given to STAR must be heterozygous-only: STAR does not ignore homozygous genotypes and uses only the first sample column.
- **Contig names must match the FASTA.** Human WGS or array VCFs often use `chr1` while the Ensembl FASTA uses `1`; STAR would then find no variant on any contig (WASP silently does nothing) and every array task would waste a full alignment. The prep job therefore compares, on a compute node and before any alignment, the contigs of every individual's heterozygous sites with the FASTA's `.fai`, and stops with a non-zero exit and a message that names the individual, the number of sites on contigs missing from the FASTA and examples of both naming styles, if **any** site lies on a contig the FASTA does not have. The skill does not rename contigs itself: if the names differ only by the `chr` prefix, tell the user to supply renamed VCFs (for example made with `bcftools annotate --rename-chrs` and a two-column map `chr1 1`, ..., `chrM MT`, in a compute job) or to use a FASTA with the same naming; if the VCF has extra contigs (alt, decoy, unplaced), restrict it to the FASTA's contigs first.

Input VCFs are read-only; the derived files go to `{RESULTS_DIR}`.

---

## Step 7 — Analysis menu

Ask (numbered, multi-select, for example "1,3"):

1. **Per-sample allelic imbalance** (always on, cannot be deselected): reference-bias diagnostic plus per-gene and per-site allelic ratios for every sample (Rmd 01 and Rmd 02).
2. **Reciprocal F1 analysis** (strain effect versus parent-of-origin effect, Rmd 03, Step 16): offered only when `{MODE}` = `f1` and at least 2 samples have `cross_direction` `AxB` and at least 2 have `BxA`. When there are several conditions, at least one condition must hold samples of both directions; otherwise condition and cross direction are confounded and Rmd 03 stops, so the option is hidden.
3. **Differential ASE between conditions** (Rmd 04, Step 17): offered only when at least two conditions each have replicates and the design allows a paired or adjusted comparison. In F1 mode, condition and cross direction must not be confounded: when both directions are present, it must not be true that every sample of one condition has one direction and every sample of the other condition has the other. In outbred mode, at least one pair of conditions with at least 2 individuals sampled in both is needed (the test pairs each individual's samples); which contrasts can be tested is checked after the reference condition is chosen (below).
4. **phASER haplotype phasing**: offered only when `{MODE}` = `outbred`; available in a later stage.

Show only the options whose preconditions hold, and say why any other is hidden, for example "Differential ASE is hidden: condition and cross direction are confounded (all ctrl samples are AxB and all treat samples BxA), so a condition effect cannot be told apart from a parent-of-origin effect." Store the selection as `{ANALYSES}`, which holds the words, separated by spaces, never the menu numbers: menu number 1 is the word `per-sample`, 2 is `reciprocal`, 3 is `differential` (`per-sample` is always present); for example the answer "1,2,3" is stored as `per-sample reciprocal differential`. Menu number 4 (phASER) runs nothing yet and is not stored in `{ANALYSES}`; `submit_chain.sh` (Step 15) stops when `{ANALYSES}` holds any other word. If the user chooses phASER, say "available in a later stage" and continue.

If differential ASE is selected, ask: "Which condition is the reference (baseline)?" (numbered list of the conditions in `{SAMPLES_CSV}`, in the order they first appear; the default, marked in the list, is the condition named control, ctrl, untreated, wt or vehicle (any case) if there is one, else the first condition alphabetically). Store it as `{REF_CONDITION}`. Guard: before Rmd 04 is written, check that the value is a condition of `{SAMPLES_CSV}` (`awk -F, -v r="{REF_CONDITION}" 'NR > 1 && $4 == r {f = 1} END {exit !f}' {SAMPLES_CSV}`) and ask again if it is not. Every other condition is compared with it, one contrast at a time. For each contrast in outbred mode, list the individuals that are sampled in only one of the two conditions: they are left out of that contrast.

The design rules of option 3 are then checked for every contrast against `{REF_CONDITION}`, the way Rmd 04 checks them (Step 17), and the outcome is told to the user before anything is submitted:

- F1: a contrast in which condition and cross direction are confounded stops the whole of Rmd 04, not only that contrast. Name the contrast; choose another reference condition first, and deselect differential ASE only if no reference gives a design without confounding.
- A contrast with fewer than 2 samples in one of its two conditions (F1), or fewer than 2 individuals sampled in both of its conditions (outbred), is not tested: Rmd 04 lists it under "Contrasts not tested" and writes no row for it. If no contrast can be tested, Rmd 04 stops: choose another reference condition first (for example one with replicates, or with individuals shared with another condition), and deselect differential ASE only if none works.
- A contrast that passes these checks can still end with nothing tested (for example no gene with 2 covered samples per condition); Rmd 04 then prints "nothing tested in contrast ..." and writes its rows with `tested` = 0.

---

## Step 8 — Constants

Show these defaults and let the user edit any of them:

| Constant | Default | Meaning |
|---|---|---|
| `MIN_DEPTH` | 10 | minimum total reads at a site or gene for testing |
| `FDR_SIG` | 0.05 | FDR threshold for significance |
| `ABS_DEV_SIG` | 0.1 | minimum absolute deviation of the reference fraction from 0.5 to call a gene imbalanced |
| `BIAS_TOL` | 0.03 | reference-bias tolerance: a sample is flagged when its reference-fraction deviation from 0.5 exceeds `max(BIAS_TOL, 3 x SE)`, SE = sqrt(0.25 / total reads) |
| `RHO_MIN` | 0.01 | both modes: floor for the per-sample overdispersion used by the tests (F1: `rho = max(rho_corrected, RHO_MIN)` at SNP and gene level; outbred: `rho = max(rho_trim, RHO_MIN)`), so a noisy or zero estimate can never make the tests anti-conservative |
| `THIN_BP` | 500 | F1 gene-level test of Rmd 02 and the gene-level tests of Rmd 03 and 04: SNPs closer than this many exonic bases can be counted in the same read pair (ASEReadCounter counts a fragment at every SNP it covers), so only the deepest SNP per window is used; set it to about the longest fragment (insert) length of the library |

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

Tools come from cached Singularity biocontainers (`{STAR_SIF}`, `{GATK_SIF}`, `{BCFTOOLS_SIF}`, `{SAMTOOLS_SIF}`, `{PICARD_SIF}` from Step 2); never `module add` anything except singularity. `bcftools`, `bgzip` and `tabix` all come from the bcftools 1.20 container. Keep only the wrappers a script uses; the wrapper lines below are one per tool, so choose them by the function name at the start of the line (never by position), with exactly these sets per script:

| Script | Wrapper functions kept (block C) | `fetch_sif` |
|---|---|---|
| `prep_f1_reference.sh`, `prep_genotypes.sh` | `star`, `gatk`, `bcftools`, `bgzip`, `tabix`, `samtools`, `picard` (all) | downloading version, called for all five containers ("Prep container fetch" below) |
| `align_count_f1.sh`, `align_wasp_count.sh` | `star`, `gatk`, `samtools`, `picard` | check-only version, called for each of the four (Step 11) |
| `extract_mgp_parental_vcf.sh`, `concat_mgp_parental_vcf.sh` | `bcftools`, `tabix` | downloading version, called for `$BCFTOOLS_SIF` only (the helper line in "Per-script fetch lines" below); the `add_bind` line is reduced to `add_bind "{CWD}"`, because the helper reads no genome file and the genome folder may not exist yet when it runs (Singularity stops on a missing bind source) |

`--bind` is required because `/net/...` paths are not auto-bound; `BIND` lists each directory once (Singularity prints "destination is already in the mount point list" for a repeated one), so add directories with `add_bind`:

```bash
module add singularity/3.10.4 || exit 1
command -v singularity >/dev/null || { echo "ERROR: singularity not on PATH" >&2; exit 1; }
STAR_SIF="{STAR_SIF}"; GATK_SIF="{GATK_SIF}"; BCFTOOLS_SIF="{BCFTOOLS_SIF}"; SAMTOOLS_SIF="{SAMTOOLS_SIF}"; PICARD_SIF="{PICARD_SIF}"
fetch_sif() {   # $1 = local path, $2 = verified URL (Step 2); downloads only when the file is missing
  [ -s "$1" ] && return 0
  mkdir -p "$(dirname "$1")"
  { wget -c -O "$1.part" "$2" && mv "$1.part" "$1"; } || { echo "ERROR: could not download $1" >&2; exit 1; }
}
BIND=""
add_bind() { case ",$BIND," in *",$1,"*) ;; *) BIND="${BIND:+$BIND,}$1" ;; esac; }
add_bind "{CWD}"; add_bind "{GENOME_DIR}"; add_bind "$(dirname "{FASTA_PATH}")"; add_bind "$(dirname "{GTF_PATH}")"
star()     { local SIF="$STAR_SIF";     singularity exec --bind "$BIND" "$SIF" STAR "$@"; }
gatk()     { local SIF="$GATK_SIF";     singularity exec --bind "$BIND" "$SIF" gatk "$@"; }
bcftools() { local SIF="$BCFTOOLS_SIF"; singularity exec --bind "$BIND" "$SIF" bcftools "$@"; }
bgzip()    { local SIF="$BCFTOOLS_SIF"; singularity exec --bind "$BIND" "$SIF" bgzip "$@"; }
tabix()    { local SIF="$BCFTOOLS_SIF"; singularity exec --bind "$BIND" "$SIF" tabix "$@"; }
samtools() { local SIF="$SAMTOOLS_SIF"; singularity exec --bind "$BIND" "$SIF" samtools "$@"; }
picard()    { local SIF="$PICARD_SIF";  singularity exec --bind "$BIND" "$SIF" picard "$@"; }
```

Defining `fetch_sif` does nothing by itself: every script must **call** it once per container it uses. The calls are not part of the code block above (a literal copy of the block must stay valid for every script): append the script's `fetch_sif` lines, from "Per-script fetch lines" below, right after the block (the calls are the lines that start with `fetch_sif "$`). The GTF directory is bound in every script (STAR `genomeGenerate` reads `{GTF_PATH}`, and with a custom GTF outside `{CWD}` and `{GENOME_DIR}` it would otherwise be invisible inside the container).

**Per-script fetch lines** (appended after block C, outside its code block; URLs from the Step 2 table):

- `prep_f1_reference.sh`, `prep_genotypes.sh`: five lines, "Prep container fetch" below.
- `align_count_f1.sh`, `align_wasp_count.sh`: four check-only calls (STAR, GATK, SAMTOOLS, PICARD), Step 11.
- `extract_mgp_parental_vcf.sh`, `concat_mgp_parental_vcf.sh`: the one line `fetch_sif "$BCFTOOLS_SIF" "https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0"`.

**Prep container fetch (both prep scripts).** The prep job (`prep_f1_reference.sh` or `prep_genotypes.sh`) fetches all five containers, not only the ones it runs itself: the per-sample array job of Step 11 only checks that its containers exist (so parallel tasks never download the same file) and needs Picard, which no prep step uses. In both prep scripts the `fetch_sif` lines appended after block C are therefore exactly:

```bash
fetch_sif "$STAR_SIF" "https://depot.galaxyproject.org/singularity/star:2.7.10b--h9ee0642_0"
fetch_sif "$GATK_SIF" "https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0"
fetch_sif "$BCFTOOLS_SIF" "https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0"
fetch_sif "$SAMTOOLS_SIF" "https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0"
fetch_sif "$PICARD_SIF" "https://depot.galaxyproject.org/singularity/picard:3.1.1--hdfd78af_0"
```

Add `add_bind "$(dirname <path>)"` for every other input directory the script reads (parental or genotype VCF, FASTQ directories). If `{FASTA_PATH}` or `{GTF_PATH}` do not exist yet (Step 4 option 1), put the Step 4 download-and-decompress commands (URLs shown to the user first) right after block C, before any command that reads the FASTA or the GTF (in `prep_f1_reference.sh` that is before its step 1, in `prep_genotypes.sh` before block R), with the same `wget -c -O FILE.part URL && mv FILE.part FILE || { echo "ERROR: ..." >&2; exit 1; }` pattern.

### Shared block R — reference files for ASEReadCounter (both modes; `FASTA="{FASTA_PATH}"` is the ORIGINAL, unmasked FASTA)

`gatk ASEReadCounter -R` needs `<fasta>.fai` and `<fasta>.dict` next to the FASTA, and counting must use the original FASTA so its sequence names match the BAM header. If the genome folder is read-only, the job fails with a clear message; copy the FASTA to a writable folder and point `{FASTA_PATH}` there.

```bash
FASTA="{FASTA_PATH}"
[ -s "$FASTA.fai" ] || samtools faidx "$FASTA" || { echo "ERROR: cannot create $FASTA.fai (read-only genome folder?)" >&2; exit 1; }
DICT="${FASTA%.*}.dict"
[ -s "$DICT" ] || gatk CreateSequenceDictionary -R "$FASTA" -O "$DICT" || { echo "ERROR: cannot create $DICT" >&2; exit 1; }
```

### Shared block I — STAR index (`INDEX_FASTA` = `$REF_DIR/masked.fa` in F1 mode, `$FASTA` in outbred mode; `STAR_INDEX` and `INDEX_KEY` set by the prep script, from `{STAR_INDEX}` of Step 5, part 2)

The prep script sets two shell variables before this block. Outbred: `STAR_INDEX="{STAR_INDEX}"; INDEX_KEY="unmasked"`. F1: `STAR_INDEX="{STAR_INDEX}_$MASK_KEY"; INDEX_KEY="$MASK_KEY"` (`MASK_KEY` from step 3 of `prep_f1_reference.sh`).

```bash
[ -s "$INDEX_FASTA.fai" ] || samtools faidx "$INDEX_FASTA" || { echo "ERROR: faidx failed for $INDEX_FASTA" >&2; exit 1; }
SA_INDEX_NBASES=$(awk '{L+=$2} END{n=int(log(L)/log(2)/2-1); if(n>14)n=14; if(n<4)n=4; print n}' "$INDEX_FASTA.fai")
# reuse only a complete index built for the same key (index_key.txt is written last, after a successful build)
if [ -s "$STAR_INDEX/SA" ] && [ -s "$STAR_INDEX/Genome" ] && [ -s "$STAR_INDEX/sjdbList.out.tab" ] \
   && [ "$(cat "$STAR_INDEX/index_key.txt" 2>/dev/null)" = "$INDEX_KEY" ]; then
  echo "Reusing STAR index $STAR_INDEX (key $INDEX_KEY)"
else
  echo "Building STAR index $STAR_INDEX (key $INDEX_KEY)"
  mkdir -p "$STAR_INDEX" || exit 1
  rm -f "$STAR_INDEX/index_key.txt"
  star --runMode genomeGenerate --genomeDir "$STAR_INDEX" --genomeFastaFiles "$INDEX_FASTA" \
       --sjdbGTFfile "{GTF_PATH}" --sjdbOverhang {SJDB_OVERHANG} --genomeSAindexNbases "$SA_INDEX_NBASES" \
       --runThreadN "${SLURM_NTASKS:-4}" --outFileNamePrefix "{RESULTS_DIR}/logs/star_index_" \
    || { echo "ERROR: STAR genomeGenerate failed" >&2; exit 1; }
  [ -s "$STAR_INDEX/SA" ] || { echo "ERROR: STAR index incomplete in $STAR_INDEX" >&2; exit 1; }
  echo "$INDEX_KEY" > "$STAR_INDEX/index_key.txt" || { echo "ERROR: cannot write $STAR_INDEX/index_key.txt" >&2; exit 1; }
fi
mkdir -p "{RESULTS_DIR}/reference" && echo "$STAR_INDEX" > "{RESULTS_DIR}/reference/star_index_path.txt" \
  || { echo "ERROR: cannot write star_index_path.txt" >&2; exit 1; }
```

An index from an older run without `index_key.txt` is rebuilt once. Two projects that share a genome folder share an index only when the key matches (outbred: same read length; F1: same strain pair, read length and masked sites).

STAR prints "Could not move Log.out" for the index build when `--outFileNamePrefix` is set; it is harmless (the log stays in `{RESULTS_DIR}/logs`). The suffix-array parameter is computed in shell from the `.fai` (light; never read the FASTA on the login node), so the wizard never substitutes it.

### `prep_f1_reference.sh` (F1 mode)

Header: `#!/bin/bash`, `#SBATCH -N 1 -p bcc`, `#SBATCH {PREP_RESOURCES}`, `#SBATCH --mail-type=END,FAIL`, `#SBATCH --mail-user={USER_EMAIL}`, `#SBATCH -o {RESULTS_DIR}/logs/prep_f1_reference_%j.out`, `set -uo pipefail`, block C with the prep container fetch below (all five wrappers are defined; plus `add_bind "$(dirname "{PARENTAL_VCF}")"`), then:

```bash
REF_DIR="{RESULTS_DIR}/reference"; mkdir -p "$REF_DIR" || exit 1
FASTA="{FASTA_PATH}"; PARENTAL="{PARENTAL_VCF}"
A="{STRAIN_A}"; B="{STRAIN_B}"
MIN_PARENTAL_SITES=1000   # stop if fewer sites remain (whole genomes have far more); use 20 for small test data
A_ALT_SITES=stop          # strain A 1/1 with strain B 0/0: stop (default) or drop (strain A is a substrain of the assembly strain; Step 6)
[ -s "$FASTA.fai" ] || samtools faidx "$FASTA" || { echo "ERROR: cannot create $FASTA.fai" >&2; exit 1; }

# 1. biallelic SNPs of the parental VCF; contig guard; with sample columns keep only strain A 0/0 + strain B 1/1
bcftools view -m2 -M2 -v snps -O v -o "$REF_DIR/parental_biallelic.vcf" "$PARENTAL" \
  || { echo "ERROR: bcftools view failed on $PARENTAL" >&2; exit 1; }
N_BI=$(grep -vc '^#' "$REF_DIR/parental_biallelic.vcf")
[ "$N_BI" -gt 0 ] || { echo "ERROR: no biallelic SNP sites in $PARENTAL" >&2; exit 1; }
# contig guard (compute node, before any alignment): every site must lie on a contig of the FASTA
N_OFF=$(grep -v '^#' "$REF_DIR/parental_biallelic.vcf" | cut -f1 | awk 'NR==FNR {c[$1]=1; next} !($1 in c)' "$FASTA.fai" - | wc -l)
if [ "$N_OFF" -gt 0 ]; then
  echo "ERROR: contig mismatch: $N_OFF of $N_BI parental SNPs lie on contigs that are not in $FASTA" >&2
  echo "  VCF contigs (first 5): $(grep -v '^#' "$REF_DIR/parental_biallelic.vcf" | cut -f1 | uniq | head -5 | paste -sd' ')" >&2
  echo "  FASTA contigs (first 5): $(cut -f1 "$FASTA.fai" | head -5 | paste -sd' ')" >&2
  echo "  Use a VCF and a FASTA with the same contig names (for example both '1' or both 'chr1'); nothing was aligned." >&2
  exit 1
fi
NS=$(bcftools query -l "$REF_DIR/parental_biallelic.vcf" | wc -l)
if [ "$NS" -eq 0 ]; then
  # sites-only VCF: REF = strain A, ALT = strain B is taken on trust (the user guarantees it, Step 6)
  bcftools view -G -O v -o "$REF_DIR/parental_snps.sites.vcf" "$REF_DIR/parental_biallelic.vcf" \
    || { echo "ERROR: bcftools view -G failed" >&2; exit 1; }
  echo "Parental VCF: $N_BI biallelic SNPs, no sample columns: all taken as REF = $A, ALT = $B (not checkable)"
else
  for S in "$A" "$B"; do
    bcftools query -l "$REF_DIR/parental_biallelic.vcf" | grep -qxF -- "$S" || { echo "ERROR: strain '$S' is not a sample of $PARENTAL (samples: $(bcftools query -l "$REF_DIR/parental_biallelic.vcf" | head -10 | paste -sd' '))" >&2; exit 1; }
  done
  # one line per site: CHROM POS REF ALT class; keep = A 0/0 and B 1/1; a_alt = A 1/1 and B 0/0;
  # same = both 0/0 or both 1/1 (the strains do not differ); het_or_missing = any other genotype
  bcftools query -s "$A,$B" -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' "$REF_DIR/parental_biallelic.vcf" \
    | awk -F'\t' -v A="$A" -v B="$B" 'BEGIN{OFS="\t"} {
        for (i = 5; i <= NF; i++) { j = index($i, "="); s = substr($i, 1, j - 1); g = substr($i, j + 1); gsub(/\|/, "/", g); gt[s] = g }
        ga = gt[A]; gb = gt[B]
        if (ga == "0/0" && gb == "1/1") c = "keep"; else if (ga == "1/1" && gb == "0/0") c = "a_alt"
        else if (ga == gb && (ga == "0/0" || ga == "1/1")) c = "same"; else c = "het_or_missing"
        print $1, $2, $3, $4, c }' > "$REF_DIR/parental_site_classes.tsv" \
    || { echo "ERROR: genotype classification failed" >&2; exit 1; }
  N_KEEP=$(awk -F'\t' '$5=="keep"' "$REF_DIR/parental_site_classes.tsv" | wc -l)
  N_SAME=$(awk -F'\t' '$5=="same"' "$REF_DIR/parental_site_classes.tsv" | wc -l)
  N_HET=$(awk -F'\t' '$5=="het_or_missing"' "$REF_DIR/parental_site_classes.tsv" | wc -l)
  N_AALT=$(awk -F'\t' '$5=="a_alt"' "$REF_DIR/parental_site_classes.tsv" | wc -l)
  echo "Parental VCF: $N_BI biallelic SNPs; kept $N_KEEP ($A 0/0, $B 1/1); dropped $N_SAME (same genotype in both strains), $N_HET (heterozygous or missing), $N_AALT ($A 1/1, $B 0/0; A_ALT_SITES=$A_ALT_SITES)"
  if [ "$N_AALT" -gt 0 ] && [ "$A_ALT_SITES" != "drop" ]; then
    echo "ERROR: $N_AALT sites have $A = 1/1 and $B = 0/0: strain A carries the non-reference allele there, so REF = strain A does not hold." >&2
    echo "  Check that $A is the reference strain and that the strains are not swapped. If $A is a substrain of the assembly strain and these sites are expected, set A_ALT_SITES=drop (Step 6)." >&2
    exit 1
  fi
  bcftools view -G -O v "$REF_DIR/parental_biallelic.vcf" \
    | awk -F'\t' 'NR==FNR { if ($5 == "keep") k[$1 FS $2 FS $3 FS $4] = 1; next } /^#/ || (($1 FS $2 FS $4 FS $5) in k)' \
        "$REF_DIR/parental_site_classes.tsv" - > "$REF_DIR/parental_snps.sites.vcf" \
    || { echo "ERROR: cannot write parental_snps.sites.vcf" >&2; exit 1; }
fi
N_SITES=$(grep -vc '^#' "$REF_DIR/parental_snps.sites.vcf")
echo "Parental SNP sites used: $N_SITES (minimum $MIN_PARENTAL_SITES)"
[ "$N_SITES" -ge "$MIN_PARENTAL_SITES" ] || { echo "ERROR: only $N_SITES parental SNP sites remain (minimum $MIN_PARENTAL_SITES); check the strain names, the genotypes and the VCF" >&2; exit 1; }

# 2. the VCF REF must equal the FASTA base at every site (right assembly). This does NOT show which strain is
#    the reference strain: that is the genotype check of step 1 (VCF with sample columns) or the user's guarantee.
bcftools norm -f "$FASTA" -c e -o /dev/null "$REF_DIR/parental_snps.sites.vcf" \
  || { echo "ERROR: the REF allele of $PARENTAL differs from the FASTA base at some sites: the VCF was not called against this assembly" >&2; exit 1; }

# 3. third-allele masking: at each site the masked genome carries the first of A, C, G, T that is
#    neither the REF nor the ALT allele (same rule everywhere), so neither strain is favoured in mapping.
#    One genotype column (sample MASK, 1/1) is written because bcftools consensus applies a sites-only VCF unreliably.
awk 'BEGIN{OFS="\t"}
  function third(r,a,  i,c,s){ s="ACGT"; for(i=1;i<=4;i++){ c=substr(s,i,1); if(c!=r && c!=a) return c } }
  /^##/ {print; next}
  /^#/  {print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"; print $1,$2,$3,$4,$5,$6,$7,$8,"FORMAT","MASK"; next}
  { print $1,$2,$3,$4,third(toupper($4),toupper($5)),".",".",".","GT","1/1" }' \
  "$REF_DIR/parental_snps.sites.vcf" > "$REF_DIR/masked_sites.vcf" || { echo "ERROR: masked sites failed" >&2; exit 1; }
# mask key: md5 of the masked site records (headers excluded), first 8 characters; names the STAR index (Step 5, part 2)
MASK_KEY=$(grep -v '^#' "$REF_DIR/masked_sites.vcf" | cut -f1,2,4,5 | md5sum | cut -c1-8)
[ -n "$MASK_KEY" ] || { echo "ERROR: cannot compute the mask key" >&2; exit 1; }
echo "Mask key: $MASK_KEY"
awk 'BEGIN{OFS="\t"} !/^#/ {print $1,$2-1,$2,$4 ">" $5}' "$REF_DIR/masked_sites.vcf" > "$REF_DIR/masked_positions.bed"
bgzip -c "$REF_DIR/masked_sites.vcf" > "$REF_DIR/masked_sites.vcf.gz" && tabix -f -p vcf "$REF_DIR/masked_sites.vcf.gz" \
  || { echo "ERROR: bgzip/tabix failed for masked_sites.vcf" >&2; exit 1; }
bcftools consensus -s MASK -f "$FASTA" "$REF_DIR/masked_sites.vcf.gz" > "$REF_DIR/masked.fa" \
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

`masked_sites.vcf` carries one sample, `MASK`, with genotype `1/1` at every site, and `bcftools consensus -s MASK` applies its ALT (third) allele. Do not use a sites-only VCF here: with bcftools 1.20 a sites-only VCF is applied unreliably (verified on the synthetic data: with `-H A` all 49 sites were applied in only 4 of 10 identical runs and none in the other 6; without `-H` none were applied in any run), while the genotyped VCF with `-s MASK` applied all 49 sites in 10 of 10 runs. The check in step 4 stops the job whenever the masking is incomplete. Then block R (`.fai` and `.dict` of the original FASTA), then block I with `INDEX_FASTA="$REF_DIR/masked.fa"; STAR_INDEX="{STAR_INDEX}_$MASK_KEY"; INDEX_KEY="$MASK_KEY"` (the masked index is keyed on the mask, Step 5, part 2), then `echo "prep_f1_reference done"`. The array job later counts against the ORIGINAL FASTA with `f1_het_sites.vcf.gz`. Third-allele masking example: REF `A`, ALT `G` gives `C`; REF `A`, ALT `C` gives `G`. Sites are the parental SNPs only (no indels), and `masked.fa` is only used to build the STAR index.

### `prep_genotypes.sh` (outbred mode)

Header as above (`prep_genotypes_%j.out`), block C with the prep container fetch below (plus `add_bind` for each genotype VCF directory), then **block R first** (the contig guard below needs `$FASTA.fai`), then one heterozygous single-sample VCF per individual. The wizard writes the `MAP` lines (individual, path, sample name inside that VCF) from `{GENOTYPE_VCFS}` and shows the two filter settings for editing: rnavar VCFs use `FILTER_EXPR='FMT/DP>=10'`; external genotypes use a genotype-quality filter such as `FILTER_EXPR='FMT/GQ>=20'`; an empty value skips the expression filter (test data). `-f PASS,.` keeps PASS and unfiltered records.

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
  # contig guard (compute node, before any alignment): every het site must lie on a contig of the FASTA
  # (chr1 in the VCF vs 1 in an Ensembl FASTA would make WASP a silent no-op after hours of alignment)
  N_ALL=$(bcftools view -H "$OUT.vcf.gz" | wc -l)
  N_OFF=$(bcftools query -f '%CHROM\n' "$OUT.vcf.gz" | awk 'NR==FNR {c[$1]=1; next} !($1 in c)' "$FASTA.fai" - | wc -l)
  if [ "$N_OFF" -gt 0 ]; then
    echo "ERROR: contig mismatch for individual $IND: $N_OFF of $N_ALL heterozygous sites lie on contigs that are not in $FASTA" >&2
    echo "  VCF contigs (first 5): $(bcftools query -f '%CHROM\n' "$OUT.vcf.gz" | uniq | head -5 | paste -sd' ')" >&2
    echo "  FASTA contigs (first 5): $(cut -f1 "$FASTA.fai" | head -5 | paste -sd' ')" >&2
    echo "  Supply VCFs with the FASTA's contig names (for example renamed with bcftools annotate --rename-chrs, Step 6); nothing was aligned." >&2
    exit 1
  fi
  tabix -f -p vcf "$OUT.vcf.gz" || { echo "ERROR: tabix failed for $OUT.vcf.gz" >&2; exit 1; }
  # plain-text copy for STAR --varVCFfile (heterozygous-only, single sample: STAR uses only the first sample)
  bcftools view -O v -o "$OUT.vcf" "$OUT.vcf.gz" || { echo "ERROR: cannot write $OUT.vcf" >&2; exit 1; }
  N=$(bcftools view -H "$OUT.vcf.gz" | wc -l)
  echo "Individual $IND: $N heterozygous SNP sites"
  [ "$N" -gt 0 ] || { echo "ERROR: individual $IND has no heterozygous SNPs after filtering" >&2; exit 1; }
  [ "$N" -ge "$MIN_SITES_WARN" ] || echo "WARNING: individual $IND has only $N heterozygous sites (fewer than $MIN_SITES_WARN)"
done < "$GENO_DIR/genotype_map.tsv"
```

The array job reads `{RESULTS_DIR}/genotypes/{individual}.het.vcf` (STAR) and `{individual}.het.vcf.gz` (ASEReadCounter and phASER; a plain VCF fails there). Block R has already run before the loop; then the gene-spans block below; then block I with `INDEX_FASTA="$FASTA"; STAR_INDEX="{STAR_INDEX}"; INDEX_KEY="unmasked"` (the unmasked genome), then `echo "prep_genotypes done"`.

### Gene spans for phASER (end of `prep_genotypes.sh`, before block I)

Written in every outbred project (a few seconds; phASER, if selected, reads it in Step 18). One row per gene: the span of the gene's exons on its contig, 0-based start, 4 columns (contig, start, stop, `gene_id`), as `phaser_gene_ae` needs. Genes whose exons lie on more than one contig or strand are left out and counted. Built from the exons, not from `gene` lines, so it works for GTFs without `gene` lines and uses the same `gene_id`s as Rmd 02.

```bash
SPAN_BED="{RESULTS_DIR}/reference/genes_span.bed"; mkdir -p "{RESULTS_DIR}/reference" || exit 1
awk -F'\t' '$3 == "exon" && match($9, /gene_id "[^"]+"/) {
    g = substr($9, RSTART + 9, RLENGTH - 10); k = $1 SUBSEP $7
    if (!(g in key)) { key[g] = k; lo[g] = $4; hi[g] = $5 } else if (key[g] != k) bad[g] = 1
    else { if ($4 < lo[g]) lo[g] = $4; if ($5 > hi[g]) hi[g] = $5 } }
  END { for (g in key) if (!(g in bad)) { split(key[g], c, SUBSEP); print c[1] "\t" lo[g] - 1 "\t" hi[g] "\t" g }
        n = 0; for (g in bad) n++; if (n > 0) print n " genes on more than one contig or strand left out" > "/dev/stderr" }' "{GTF_PATH}" \
  | sort -k1,1 -k2,2n > "$SPAN_BED" || { echo "ERROR: cannot write $SPAN_BED" >&2; exit 1; }
N_SPAN=$(wc -l < "$SPAN_BED")
N_SPAN_ON=$(awk 'NR == FNR {c[$1] = 1; next} ($1 in c)' "$FASTA.fai" "$SPAN_BED" | wc -l)
[ "$N_SPAN" -gt 0 ] && [ "$N_SPAN_ON" -gt 0 ] \
  || { echo "ERROR: no gene span of {GTF_PATH} lies on a contig of $FASTA (GTF contigs: $(cut -f1 "$SPAN_BED" | uniq | head -3 | paste -sd' '); FASTA contigs: $(cut -f1 "$FASTA.fai" | head -3 | paste -sd' '))" >&2; exit 1; }
echo "Gene spans for phASER: $N_SPAN genes, $N_SPAN_ON on contigs of the FASTA"
```

### `extract_mgp_parental_vcf.sh` (F1 mode, mouse helper; optional)

For the two named strains (`{STRAIN_A}` = the reference strain, for example `C57BL_6NJ`, and `{STRAIN_B}`, for example `A_J`; names exactly as in the VCF header) the helper queries the Mouse Genomes Project VCF remotely. Remote access costs about 80 s fixed overhead plus about 80 s per 10 Mb of region, so a whole-genome serial run would take about 6 hours: it is therefore **one array task per chromosome** (each well under 4 h) followed by a final dependent `bcftools concat` job. The wizard must verify the URL with a HEAD request in the session before writing it (this was verified on 2026-09-29: HTTP 200, `ase-pipeline/tests/fixtures/verified_urls.txt`):

```bash
curl -sI "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz" | head -1
```

The release is GRCm39 and its contigs are `1`, `2`, ... (no `chr` prefix), so the FASTA must be GRCm39 with the same names. **The helper needs neither the FASTA nor its `.fai`**: it runs before the prep job, which is the job that downloads or decompresses the FASTA and writes the `.fai`, so the helper must never read them (and the wizard never reads the FASTA on the login node). The wizard writes `{RESULTS_DIR}/mgp/chromosomes.txt` from the fixed GRCm39 list, one name per line: `1` ... `19`, `X` (20 lines; Y and MT are not queried), and sets `{N_CHROM}` = 20. Each array task checks on its compute node that its chromosome is a contig of the Mouse Genomes Project VCF header; the comparison of the parental VCF's contigs with the FASTA happens in the prep job (contig guard, step 1 of `prep_f1_reference.sh`), on a compute node, once the FASTA and `.fai` exist and before any alignment. Header: `#SBATCH -N 1 -p bcc`, `#SBATCH --array=1-{N_CHROM}`, `#SBATCH -n 2 --mem=8G -t 4:00:00`, mail lines, `#SBATCH -o {RESULTS_DIR}/logs/extract_mgp_%A_%a.out`, `set -uo pipefail`, block C with the helper's bind line (see the table in block C: `add_bind "{CWD}"` only, never the genome folder), plus `add_bind "{RESULTS_DIR}"`, then:

```bash
URL="https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz"
A="{STRAIN_A}"; B="{STRAIN_B}"
MGP_DIR="{RESULTS_DIR}/mgp"; OUT_VCF="{PARENTAL_VCF}"      # {RESULTS_DIR}/mgp/parental_{STRAIN_A}_{STRAIN_B}.vcf.gz (written by the concat job)
CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MGP_DIR/chromosomes.txt")
[ -n "$CHR" ] || { echo "ERROR: no chromosome for array task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
# this task's chromosome must be a contig of the MGP header (the FASTA is not read here; the prep job compares it)
bcftools view -h "$URL" > "$MGP_DIR/header_$SLURM_ARRAY_TASK_ID.txt" || { echo "ERROR: cannot read the MGP VCF header from $URL" >&2; exit 1; }
grep -q "^##contig=<ID=${CHR%%:*}," "$MGP_DIR/header_$SLURM_ARRAY_TASK_ID.txt" \
  || { echo "ERROR: chromosome '${CHR%%:*}' is not a contig of the MGP VCF (its contigs are 1, 2, ..., X without chr); fix chromosomes.txt" >&2; exit 1; }
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

**`concat_mgp_parental_vcf.sh`** header: `#SBATCH -N 1 -p bcc`, `#SBATCH -n 1 --mem=4G -t 1:00:00`, mail lines, `#SBATCH -o {RESULTS_DIR}/logs/concat_mgp_%j.out`, `set -uo pipefail`, block C with the helper's bind line (`add_bind "{CWD}"` only), plus `add_bind "{RESULTS_DIR}"`, then:

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

Both scripts are SLURM array jobs (`#SBATCH --array=1-{ARRAY_N}`, `#SBATCH {ARRAY_RESOURCES}`), take row `SLURM_ARRAY_TASK_ID` of `{SAMPLES_CSV}` (data row 1 = task 1), and write the outputs below. Header: `#!/bin/bash`, `#SBATCH -N 1 -p bcc`, `#SBATCH --array=1-{ARRAY_N}`, `#SBATCH {ARRAY_RESOURCES}`, `#SBATCH --mail-type=END,FAIL`, `#SBATCH --mail-user={USER_EMAIL}`, and one log per array task, `%A` = array job id and `%a` = task id: `#SBATCH -o {RESULTS_DIR}/logs/align_count_f1_%A_%a.out` (F1) or `#SBATCH -o {RESULTS_DIR}/logs/align_wasp_count_%A_%a.out` (outbred). They use block C, but the containers must already be cached by the prep job (it fetches all five, see "Prep container fetch" in Step 10), so parallel tasks never download the same file: in these scripts define `fetch_sif` as a check only, `fetch_sif() { [ -s "$1" ] || { echo "ERROR: missing container $1 (run the prep job first)" >&2; exit 1; }; }`, and **call** it for each container the script uses, on the lines right after block C (see "Per-script fetch lines" in Step 10): `fetch_sif "$STAR_SIF"; fetch_sif "$GATK_SIF"; fetch_sif "$SAMTOOLS_SIF"; fetch_sif "$PICARD_SIF"`. The wrapper functions are `star`, `gatk`, `samtools` and `picard` (the table in Step 10, block C). `{MIN_MAPQ}` and `{MIN_BASEQ}` are the ASEReadCounter thresholds recorded in Step 8.

| Output | Content |
|---|---|
| `{RESULTS_DIR}/bam/{sample}.bam` (+ `.bai`) | coordinate-sorted, duplicates marked (not removed: ASEReadCounter skips them), read group kept |
| `{RESULTS_DIR}/ase_counts/{sample}.table` | ASEReadCounter table: `contig position variantID refAllele altAllele refCount altCount totalCount lowMAPQDepth lowBaseQDepth rawDepth otherBases improperPairs` |
| `{RESULTS_DIR}/ase_counts/{sample}.wasp_stats.tsv` | outbred only: alignments by `vW` and first `vA` value before filtering |
| `{RESULTS_DIR}/ase_counts/{sample}.unfiltered.table` | outbred only: the same ASEReadCounter count on the STAR BAM BEFORE the WASP filter (duplicates marked the same way), for the before/after REF fraction in Rmd 01; reported for information, not a gate. Rmd 01 reads it per sample from `{SAMPLES_CSV}` and never mistakes it for a sample table |

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
# the prep job records the index it built or reused (F1: keyed on the mask, Step 5, part 2)
STAR_INDEX=$(cat "$R/reference/star_index_path.txt" 2>/dev/null)
[ -n "$STAR_INDEX" ] && [ -s "$STAR_INDEX/SA" ] || die "no STAR index recorded in $R/reference/star_index_path.txt (run the prep job first)"
STAR_DIR="$R/star/$SAMPLE"; TMPD="$R/tmp/$SAMPLE"; BAM="$R/bam/$SAMPLE.bam"; TABLE="$R/ase_counts/$SAMPLE.table"; STATS="$R/ase_counts/$SAMPLE.wasp_stats.tsv"
UNF_TABLE="$R/ase_counts/$SAMPLE.unfiltered.table"   # outbred only (pre-WASP counts); never written in F1 mode
mkdir -p "$STAR_DIR" "$TMPD" "$R/bam" "$R/ase_counts" || die "cannot create output directories"
rm -rf "$STAR_DIR/_STARtmp"
rm -f "$TABLE" "$STATS"      # never keep a table or stats file from an earlier run: a failed sample must have none
rm -f "$UNF_TABLE"
[ -s "$FQ1" ] || die "missing $FQ1"
```

### `align_count_f1.sh` (F1 mode)

After the common start (`--readFilesIn` takes one file for single-end data and two for paired-end data; the wizard writes `"$FQ1" "$FQ2"` or only `"$FQ1"` accordingly):

```bash
star --runThreadN "${SLURM_NTASKS:-4}" --genomeDir "$STAR_INDEX" \
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
star --runThreadN "${SLURM_NTASKS:-4}" --genomeDir "$STAR_INDEX" \
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

# the same count on the UNFILTERED STAR BAM (before WASP), same flags, duplicates marked the same way, so the
# two tables differ only by the WASP filter; Rmd 01 reports the before/after REF fraction for information, not as a gate
UNF_BAM="$STAR_DIR/unfiltered.markdup.bam"
picard MarkDuplicates I="$RAW" O="$UNF_BAM" M="$STAR_DIR/unfiltered.markdup_metrics.txt" \
       TMP_DIR="$TMPD" REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT || die "MarkDuplicates (unfiltered BAM) failed"
samtools index "$UNF_BAM" || die "samtools index (unfiltered BAM) failed"
gatk ASEReadCounter -R "$FASTA" -I "$UNF_BAM" -V "$HET_GZ" \
     --min-mapping-quality {MIN_MAPQ} --min-base-quality {MIN_BASEQ} \
     --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE \
     --tmp-dir "$TMPD" -O "$UNF_TABLE" || die "ASEReadCounter (unfiltered BAM) failed"
```

Never use `--outSAMattributes` values other than the explicit list above: WASP needs `vA vG vW` spelled out. `vG` is 0-based. WASP does not promise less bias; the stats file lets Rmd 01 show how many alignments were removed and from which allele, and `{sample}.unfiltered.table` lets it show the site-level REF fraction before and after the filter (for information, not a gate).

### Empty-table guard (last block of both scripts)

The tool exits 0 even when it counted nothing, so the script itself fails, and an empty table can never reach Rmd 01. The same guard covers the outbred unfiltered table:

```bash
# ase_counts/$SAMPLE.table needs non-empty data rows (header only = ASEReadCounter counted nothing)
check_table() {   # $1 = count table; a missing or header-only table is deleted and the task fails
  local N_ROWS=0
  [ -s "$1" ] && N_ROWS=$(awk 'NR>1' "$1" | wc -l)
  if [ "$N_ROWS" -lt 1 ]; then
    echo "ERROR: sample $SAMPLE: $1 is missing or has no data rows (check the read group, the het-sites VCF and the contig names)" >&2
    rm -f "$1"            # a header-only table must not stay on disk
    exit 1
  fi
  echo "Sample $SAMPLE: $N_ROWS sites counted in $(basename "$1")"
}
check_table "$TABLE"
check_table "$UNF_TABLE"     # outbred (align_wasp_count.sh) only; omit this line in align_count_f1.sh
```

### Submission order

Submission and waiting are done by `submit_chain.sh` and `wait_chain.sh` of Step 15 (job ids checked, `afterok` dependencies, bounded wait); see Step 15 and do not submit these scripts by hand. Rmd 01 (Step 12) reads `{RESULTS_DIR}/ase_counts/{sample}.table` (and `{sample}.wasp_stats.tsv` and `{sample}.unfiltered.table` in outbred mode) for every sample in `{SAMPLES_CSV}`, by name.

---

## Step 12 — Rmd 01: import, site filters and reference-bias QC

Write `{CWD}/{TODAY}_{WD_NAME}_01_import_qc.Rmd` (ask once for `{AUTHOR}` and `{PROJECT_TITLE}`; the same two values go into Rmd 02). The Rmd is self-contained (it sources no helper files), uses `knitr::opts_chunk$set(cache = FALSE)` and `options(scipen = 9)`, loads the Bioconductor packages before `tidyverse`, and calls every dplyr verb with the `dplyr::` prefix. Substitute the placeholders from Steps 0-9 (`{MODE}`, `{STRAIN_A}`, `{STRAIN_B}`, `{SAMPLES_CSV}`, `{RESULTS_DIR}`, `{GTF_PATH}`, and the four constants of Step 8 as bare numbers); for outbred mode leave `{STRAIN_A}` and `{STRAIN_B}` as the words `REF` and `ALT`. The constants block below is the only place the four thresholds are defined; every table and figure uses these objects.

**It reads `ase_counts/{sample}.table` for every sample listed in `samples.csv`, never by globbing `*.table`.** A stray or stale table from an earlier run is therefore never accepted, and the Rmd stops with an error that names every sample whose table is missing, header-only or without data rows (the Step 11 scripts delete failed tables; this is the second line of defence). In outbred mode the same rule applies to `ase_counts/{sample}.wasp_stats.tsv` and `ase_counts/{sample}.unfiltered.table` (the pre-WASP counts). Because tables are read by name, a `{sample}.unfiltered.table` is never picked up as a sample table: the stray-table message skips it, and the Rmd checks that the samples read are exactly the rows of `samples.csv`.

**Every table is read with explicit column classes** (`colClasses`: `contig`, `variantID`, `refAllele` and `altAllele` as character in the count tables, all columns as character in `samples.csv`, fixed classes in `wasp_stats.tsv`). With Ensembl names a table whose contigs are all `1`..`22` would otherwise be read as integer and a table with `X` or `MT` as character, and binding the samples would stop with "Can't combine `contig` <integer> and <character>" (for example a male without expressed X heterozygous sites next to females with them).

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
read_counts <- function(s, suffix = ".table") {   # suffix ".unfiltered.table" = outbred pre-WASP counts
  f <- file.path(COUNTS_DIR, paste0(s, suffix))
  if (!file.exists(f)) return(list(data = NULL, problem = "table missing"))
  # text columns are read as character in EVERY table: a table whose contigs are all 1..22 would otherwise be read as
  # integer and fail to bind with one that has X or MT (and an all-"T" allele column would become logical TRUE)
  d <- tryCatch(read.delim(f, stringsAsFactors = FALSE, check.names = FALSE,
                           colClasses = c(contig = "character", variantID = "character",
                                          refAllele = "character", altAllele = "character")), error = function(e) NULL)
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
  ures <- lapply(samples$sample, read_counts, suffix = ".unfiltered.table"); names(ures) <- samples$sample
  uprob <- vapply(ures, function(r) if (is.na(r$problem)) NA_character_ else paste("unfiltered", r$problem), character(1))
  wprob <- ifelse(is.na(wprob), uprob, ifelse(is.na(uprob), wprob, paste(wprob, uprob, sep = "; ")))
  problems <- ifelse(is.na(problems), wprob, ifelse(is.na(wprob), problems, paste(problems, wprob, sep = "; ")))
}
problems <- problems[!is.na(problems)]
if (length(problems) > 0)
  stop(length(problems), " sample(s) have no usable ASE count table in ", COUNTS_DIR, ":\n",
       paste0("  - ", names(problems), ": ", problems, collapse = "\n"),
       "\nRe-run the per-sample array script for these samples; tables of samples not listed in samples.csv are ignored.",
       call. = FALSE)
tbl_files <- list.files(COUNTS_DIR, pattern = "\\.table$")   # only to report stray tables; nothing is read from this list
tbl_files <- tbl_files[!grepl("\\.unfiltered\\.table$", tbl_files)]   # pre-WASP tables are not sample tables
extra <- setdiff(sub("\\.table$", "", tbl_files), samples$sample)
if (length(extra) > 0) message("Ignoring tables of samples not in samples.csv: ", paste(extra, collapse = ", "))
sites_raw <- dplyr::bind_rows(lapply(res, function(r) r$data))
stopifnot(setequal(unique(sites_raw$sample), samples$sample))   # one table per samples.csv row, read by name, never globbed
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

For every sample the mean REF fraction is REF reads divided by REF + ALT reads summed over the filtered sites. A sample is flagged when its deviation from 0.5 exceeds `max(BIAS_TOL, 3 * SE)` with `SE = sqrt(0.25 / total reads)`. **This is a screen, not a test.** On data with planted or real allelic imbalance the all-sites mean is not expected to be 0.5, so a flag says "look at the ratios with caution", not "the sample is wrong"; a benign spread of about +0.02 to -0.03 was seen on synthetic data with no imbalance at about 2000 reads per sample. In outbred mode the table also shows the site-level REF fraction before (`ref_frac_before_wasp`, counted on the unfiltered STAR BAM) and after (`ref_frac_after_wasp`) the WASP filter, on the sites kept in both tables; this is reported for information, not a gate (the flag and every test use the WASP-filtered counts only).

```{r bias}
bias <- sites %>% dplyr::group_by(sample) %>%
  dplyr::summarise(n_sites = dplyr::n(), total_reads = sum(total), ref_reads = sum(ref_n), .groups = "drop") %>%
  dplyr::mutate(mean_ref_frac = ref_reads / total_reads, deviation = mean_ref_frac - 0.5,
                SE = sqrt(0.25 / total_reads), threshold = pmax(BIAS_TOL, 3 * SE),
                flagged = abs(deviation) > threshold)
# outbred: site-level REF fraction before (unfiltered STAR BAM) and after the WASP filter, on the sites that pass the
# same site filters in both tables, so the two numbers differ only by the filter. Reported for information, not a gate.
bias$sites_before_after <- NA_integer_; bias$ref_frac_before_wasp <- NA_real_; bias$ref_frac_after_wasp <- NA_real_
if (MODE == "outbred") {
  site_filter <- function(d) dplyr::filter(d, totalCount >= MIN_DEPTH, (lowMAPQDepth + lowBaseQDepth + otherBases) <= 0.1 * rawDepth, improperPairs == 0)
  unf <- site_filter(dplyr::bind_rows(lapply(ures, function(r) r$data)))
  ba <- dplyr::inner_join(dplyr::select(unf, sample, contig, position, ref_u = refCount, alt_u = altCount),
                          dplyr::select(s3, sample, contig, position, ref_w = refCount, alt_w = altCount),
                          by = c("sample", "contig", "position")) %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(sites_before_after = dplyr::n(), ref_frac_before_wasp = sum(ref_u) / sum(ref_u + alt_u),
                     ref_frac_after_wasp = sum(ref_w) / sum(ref_w + alt_w), .groups = "drop")
  bias <- dplyr::left_join(dplyr::select(bias, -sites_before_after, -ref_frac_before_wasp, -ref_frac_after_wasp), ba, by = "sample")
}
knitr::kable(bias, digits = 4, caption = paste0("Reference-bias diagnostic (REF = ", ref_label, "). ",
  "flagged uses mean_ref_frac only. Outbred: ref_frac_before_wasp / ref_frac_after_wasp = REF fraction on the sites_before_after ",
  "sites kept in both the unfiltered and the WASP-filtered table, reported for information, not a gate"))
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

Render each Rmd from its own `sbatch -p bcc` script, written to `{RESULTS_DIR}/scripts/` (like every other script of this skill). `run_01_import_qc.sh` is the template; `run_02_imbalance.sh` (Step 13) is identical except for the job name, the log name and the Rmd file. Requests: `-n 1 --mem=16G -t 1:00:00` (both Rmds run single-threaded). The log goes to the shared filesystem under `{RESULTS_DIR}/logs` (never `/tmp`, which is node-local). `{GTF_DIR}` is the directory of `{GTF_PATH}` (`dirname`); bind each directory once. **`--bind` line, two cases** (use exactly one; comma, no space, never a path glued to another):

- `{GTF_DIR}` inside `{CWD}`: `singularity exec --bind {CWD} {R_SIF} \` (for example `--bind /data/proj` when the GTF is `/data/proj/genome/genes.gtf`);
- `{GTF_DIR}` outside `{CWD}`: `singularity exec --bind {CWD},{GTF_DIR} {R_SIF} \` (for example `--bind /data/proj,/refs/mouse/mm39_ens112`).

The template below shows the second case:

```bash
#!/bin/bash
#SBATCH -J ase_01_import_qc
#SBATCH -N 1 -p bcc
#SBATCH -n 1 --mem=16G -t 1:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={USER_EMAIL}
#SBATCH -o {RESULTS_DIR}/logs/run_01_import_qc_%j.out
set -uo pipefail
module add singularity/3.10.4 || exit 1
singularity exec --bind {CWD},{GTF_DIR} {R_SIF} \
  Rscript -e "rmarkdown::render('{CWD}/{TODAY}_{WD_NAME}_01_import_qc.Rmd', output_dir = '{RESULTS_DIR}', knit_root_dir = '{CWD}')" \
  || { echo "ERROR: Rmd 01 failed" >&2; exit 1; }
```

`{RESULTS_DIR}` is inside `{CWD}`, so one bind covers the Rmd, the count tables and the output. The stage HTML is written to `{RESULTS_DIR}/{TODAY}_{WD_NAME}_01_import_qc.html`. The full submission order (with dependencies) is in Step 15.

If the job stops with "sample(s) have no usable ASE count table", show the listed samples to the user; do not delete anything, and do not work around it by globbing. Either fix and re-run the failed samples, or drop them as described in Step 15 ("Dropping a failed sample").

---

## Step 13 — Rmd 02: per-sample allelic imbalance

Write `{CWD}/{TODAY}_{WD_NAME}_02_imbalance.Rmd` with the same header conventions as Step 12. It loads `ase_checkpoint.rds`, defines its own constants block (same names, same values; the Rmd stops if they differ from the checkpoint, so Rmd 01 and Rmd 02 can never disagree), pastes the statistics functions of Step 14 (the code between the two marker lines, without the marker lines) into the chunk marked below, and writes `{TODAY}_{WD_NAME}_ASE_imbalance.xlsx` (sheets `SNP`, `Gene`, `Summary`) and `ase_imbalance_checkpoint.rds` to `{RESULTS_DIR}`.

- Per sample, `rho_h0 <- bb_estimate_rho(alt, total)` is estimated under H0 (p = 0.5) from all filtered sites. Every truly imbalanced site inflates it (on the synthetic acceptance data it was 0.08-0.10 against a binomial truth and detected 0 of 16 planted genes), so no test uses it except the F1 fallback below; it is reported as the diagnostic `rho_h0_naive`. **Outbred** uses `bb_estimate_rho_trim(alt, total)`: the H0 fit restricted to the central sites (counts folded around n/2; each site's central region holds 90 percent of the H0 probability), with the truncation corrected in the likelihood (Step 14); the tests use `rho_used = max(rho_trim, RHO_MIN)`. Unphased data allow no free mean per gene (the SNP orientation is unknown), which is why the outbred estimator trims instead. **F1** uses `bb_estimate_rho_gene(alt, total, gene)` on all the sample's SNPs in GTF genes, estimated with a free mean per gene and bias-corrected (`rho_corrected`), because real imbalance inflates `rho_h0`; the SNP-level tests use `rho_used = max(rho_corrected, RHO_MIN)`; if it is `NA` (fewer than 5 genes with 2 or more SNPs) the Rmd prints a WARNING, falls back to `rho_h0` and records that in the Summary (`rho_source`). The F1 gene LRT has its own `rho_gene`, with this fallback chain: (1) `max(rho_corrected, RHO_MIN)` of `bb_estimate_rho_gene` on the thinned SNPs (the SNPs the gene test uses); (2) if that is `NA` (fewer than 5 genes with 2 or more thinned SNPs), the unthinned free-mean value `max(rho_corrected, RHO_MIN)` used by the SNP tests; (3) only if both are `NA`, `max(rho_h0, RHO_MIN)`. The Summary records which one was used (`rho_gene_source`); when it is the unthinned fallback (2), note that the fallback with the unthinned dispersion has not been verified when the true overdispersion is above RHO_MIN; because SNPs that share reads push that estimate low, the gene test may be anti-conservative in that case. The Summary shows all the estimates, and any estimator boundary warning per sample. With fewer than 20 sites `rho_h0` is `NA` and the sample is not tested.
- Per SNP, `bb_pvalue`, then Benjamini-Hochberg within each sample. The single column `sig` is `padj < FDR_SIG` and `|ALT fraction - 0.5| >= ABS_DEV_SIG` (F1: ALT is strain B). The Summary table and every plot use this column and nothing else; the Rmd checks that the counts drawn in the figures equal the Summary counts.
- Gene level, both modes: SNP positions are overlapped with the GTF exons by `GenomicRanges::findOverlaps` (`gene_id` from the GTF). **F1:** counts are not summed before testing (summing and applying the per-SNP `rho` to the total inflates the variance by about `1 + (n - 1) * rho` and destroys power); each gene's SNPs are first thinned with `thin_snps` (at most one SNP per `THIN_BP` window in exon coordinates, deepest SNP first, with the depth pooled over all samples, so every sample uses the same choice of SNPs), so no read pair is counted twice; the thinned SNPs share one strain-B fraction and `bb_gene_lrt` tests it against 0.5 with `rho_gene`; the table reports `n_snps` (the gene's SNPs in the sample), `n_snps_used` (the thinned SNPs the test used; the read sums are over these), `phat` (strain-B fraction), p, BH within the sample, and `sig` from `FDR_SIG` and `ABS_DEV_SIG` on `|phat - 0.5|`, plus a `note` column saying that the test is thinned: at most one SNP per `THIN_BP` window (exon coordinates) is used, so no read pair is counted twice, which is approximate because alternative splicing can bring distant exons into one fragment (Step 14). **Outbred:** SNPs cannot be pooled without phasing, so the gene p-value is the `acat` combination of the SNP p-values, labelled "unphased, no direction" (no direction column); the gene is `sig` when its BH-adjusted `acat` p-value is below `FDR_SIG` and at least one of its SNPs deviates by `ABS_DEV_SIG` or more.
- **X chromosome (F1).** Rmd 01 and Rmd 02 do not set X aside (Rmd 03 and Rmd 04 do, Step 14): male F1 animals carry one X, so its genes are monoallelic (the mother's strain) and can look strongly imbalanced in the gene and SNP tables, and the X sites shift `mean_ref_frac` (up in AxB males, down in BxA males), which can trip the reference-bias flag without any bias. Tell the user to read X rows with care; say this in the summary page when a male sample is flagged.
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
DATE_TAG    <- "{TODAY}_{WD_NAME}"
MIN_DEPTH   <- {MIN_DEPTH}
FDR_SIG     <- {FDR_SIG}
ABS_DEV_SIG <- {ABS_DEV_SIG}
RHO_MIN     <- {RHO_MIN}             # both modes: floor for the overdispersion used by the tests
THIN_BP     <- {THIN_BP}             # F1 gene test: SNPs closer than this (exon coordinates) can share a read pair: one kept per window
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
# F1 gene test: at most one SNP per THIN_BP window of each gene (exon coordinates), deepest first, depth pooled over all samples.
# The exon-coordinate table (sg_exon) is used only for thinning; snp_gene above keeps its rows and order (per-SNP rho labels, SNP sheet).
exd <- as.data.frame(GenomicRanges::reduce(split(ex, ex$gene_id)))              # exon union per gene: group_name = gene_id
exd <- exd[order(exd$group_name, exd$start), ]
exd$offset <- ave(exd$width, exd$group_name, FUN = function(w) cumsum(w) - w)    # exonic bases of the gene before this exon
exr <- GenomicRanges::GRanges(as.character(exd$seqnames), IRanges::IRanges(exd$start, exd$end))
eh <- GenomicRanges::findOverlaps(gr, exr); qh <- S4Vectors::queryHits(eh); sh <- S4Vectors::subjectHits(eh)
sg_exon <- data.frame(contig = pos$contig[qh], position = pos$position[qh], gene_id = exd$group_name[sh],
                      exon_pos = exd$offset[sh] + pos$position[qh] - exd$start[sh] + 1, stringsAsFactors = FALSE) %>%
  dplyr::distinct(contig, position, gene_id, .keep_all = TRUE)                  # one row per SNP and gene, so the join below keeps snp_gene's rows
thin_depth <- dplyr::inner_join(sites, sg_exon, by = c("contig", "position")) %>%
  dplyr::group_by(gene_id, contig, position, exon_pos) %>% dplyr::summarise(depth = sum(total), .groups = "drop") %>%
  dplyr::group_by(gene_id) %>% dplyr::mutate(thin_keep = thin_snps(exon_pos, depth, THIN_BP)) %>% dplyr::ungroup()
snp_gene <- dplyr::left_join(snp_gene, dplyr::select(thin_depth, gene_id, contig, position, thin_keep), by = c("gene_id", "contig", "position"))
snp_gene$thin_keep[is.na(snp_gene$thin_keep)] <- FALSE
cat(nrow(pos), "SNP positions,", length(unique(snp_gene$gene_id)), "genes with at least one SNP;", sum(snp_gene$thin_keep),
    "of", nrow(snp_gene), "SNP-gene pairs kept after thinning (THIN_BP", THIN_BP, ")\n")
```

## Per-SNP beta-binomial test

The overdispersion `rho` is estimated per sample. The plain H0 fit on all filtered sites (`rho_h0`, `bb_estimate_rho`) absorbs the spread of every truly imbalanced site; it is shown for comparison only (Summary column `rho_h0_naive`), because testing with it removes almost all power (0 of 16 planted genes on the synthetic acceptance data). **Outbred:** `bb_estimate_rho_trim` fits the H0 beta-binomial on the central sites only: each site's central region is the set of counts closest to n/2 that holds 90 percent of the H0 probability under the current `rho`, sites whose count falls outside it are set aside, and each kept site's likelihood is divided by the H0 probability of its region (a truncated likelihood, so trimming the tails does not shrink `rho`); region and `rho` are iterated upward from a small start to a fixed point. The tests use `rho_used = max(rho_trim, RHO_MIN)`, and `central_frac` (the fraction of sites kept) is reported. On simulated unphased data (Step 14) this keeps the null size of the SNP and ACAT gene tests at or below about 0.05 and restores the power; with real overdispersion and many imbalanced sites `rho_trim` stays somewhat inflated, which lowers power but does not add false positives. **F1:** the H0-based value (`rho_h0`) is inflated by every truly imbalanced gene, so `rho` is instead estimated with a free mean per gene (`bb_estimate_rho_gene`, genes with at least 2 SNPs). That estimate (`rho_free`) is biased low, because a free mean per gene absorbs part of the variance (by about (k - 1) / k for genes with k SNPs), so it is bias-corrected (`rho_corrected`, see Step 14). The SNP-level tests use `rho_used = max(rho_corrected, RHO_MIN)` (`RHO_MIN` from Step 8), so a near-boundary corrected estimate can never make them anti-conservative. If the free-mean estimate is not available (fewer than 5 usable genes) the H0-based value is used and a WARNING is printed and written to the Summary (`rho_source`). The F1 gene LRT uses its own `rho_gene`: the same estimator on the thinned SNPs of the gene map above (the SNPs the gene test uses), floored at `RHO_MIN`; if fewer than 5 genes keep 2 or more thinned SNPs it falls back to the unthinned value used by the SNP tests (a NOTE is printed), and only if that is also `NA` to `max(rho_h0, RHO_MIN)`; `rho_gene_source` in the Summary says which. Any estimator warning (for example an estimate at the upper boundary) is printed per sample and written to `rho_warning`.

```{r snp}
bb_p_safe <- function(x, n, rho) if (is.na(rho)) NA_real_ else bb_pvalue(x, n, rho = rho)
collect_warnings <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") })
  list(value = v, warn = paste(unique(w), collapse = "; "))
}
snp <- dplyr::bind_rows(lapply(split(sites, sites$sample), function(d) {
  r0 <- collect_warnings(bb_estimate_rho(d$alt_n, d$total)); rho_h0 <- r0$value; warn <- r0$warn
  rho_free <- NA_real_; rho_corrected <- NA_real_; rho_trim <- NA_real_; central_frac <- NA_real_; rho <- rho_h0; src <- "H0-based"
  if (MODE == "outbred") {
    r2 <- collect_warnings(bb_estimate_rho_trim(d$alt_n, d$total))
    warn <- paste(c(warn, r2$warn)[nzchar(c(warn, r2$warn))], collapse = "; ")
    rho_trim <- r2$value[["rho"]]; central_frac <- r2$value[["central_frac"]]
    rho <- if (is.na(rho_trim)) NA_real_ else max(rho_trim, RHO_MIN)   # NA (fewer than 20 sites): sample not tested
    src <- "central sites, truncation-corrected (trimmed H0 fit), floored at RHO_MIN"
  }
  rho_gene <- NA_real_; gsrc <- NA_character_
  if (MODE == "f1") {
    g <- snp_gene$gene_id[match(paste(d$contig, d$position), paste(snp_gene$contig, snp_gene$position))]
    r1 <- collect_warnings(bb_estimate_rho_gene(d$alt_n, d$total, g))
    warn <- paste(c(warn, r1$warn)[nzchar(c(warn, r1$warn))], collapse = "; ")
    if (!is.na(r1$value[["corrected"]])) {
      rho_free <- r1$value[["free"]]; rho_corrected <- r1$value[["corrected"]]
      rho <- max(rho_corrected, RHO_MIN)   # floor at SNP level (F1)
      src <- "free-mean per gene, bias-corrected"
    } else {
      cat("WARNING: sample", d$sample[1], "has fewer than 5 genes with 2 or more SNPs; using the H0-based rho\n"); src <- "H0-based (fallback)" }
    # gene LRT: its own rho_gene from the THINNED SNPs; fallback chain thinned -> unthinned free-mean -> H0-based
    sg_keep <- snp_gene[snp_gene$thin_keep, ]
    g_thin <- sg_keep$gene_id[match(paste(d$contig, d$position), paste(sg_keep$contig, sg_keep$position))]
    rt <- collect_warnings(bb_estimate_rho_gene(d$alt_n, d$total, g_thin))
    if (!is.na(rt$value[["corrected"]])) {
      rho_gene <- max(rt$value[["corrected"]], RHO_MIN); gsrc <- "thinned SNPs, free-mean per gene, bias-corrected"
    } else if (!is.na(rho_corrected)) {
      rho_gene <- max(rho_corrected, RHO_MIN); gsrc <- "unthinned SNPs, free-mean per gene, bias-corrected (fallback: fewer than 5 genes with 2 or more thinned SNPs)"
      cat("NOTE: sample", d$sample[1], "has fewer than 5 genes with 2 or more thinned SNPs; the gene test uses the unthinned free-mean rho;",
          "the fallback with the unthinned dispersion has not been verified when the true overdispersion is above RHO_MIN; because SNPs that share reads push that estimate low, the gene test may be anti-conservative in that case\n")
    } else {
      rho_gene <- if (is.na(rho_h0)) NA_real_ else max(rho_h0, RHO_MIN); gsrc <- "H0-based (fallback)"
    }
    if (nzchar(rt$warn)) { cat("WARNING (thinned rho_gene), sample", d$sample[1], ":", rt$warn, "\n"); gsrc <- paste0(gsrc, "; warning: ", rt$warn) }
  }
  if (nzchar(warn)) cat("WARNING (rho estimation), sample", d$sample[1], ":", warn, "\n")
  d$rho_h0 <- rho_h0; d$rho_free <- rho_free; d$rho_corrected <- rho_corrected; d$rho_trim <- rho_trim; d$central_frac <- central_frac; d$rho <- rho
  d$rho_gene <- rho_gene   # F1 gene LRT only (NA for outbred)
  d$rho_source <- src; d$rho_warning <- warn; d$rho_gene_source <- gsrc
  d$p <- mapply(bb_p_safe, d$alt_n, d$total, MoreArgs = list(rho = rho))
  d$padj <- p.adjust(d$p, method = "BH")
  d
})) %>%
  dplyr::mutate(dev = abs(alt_frac - 0.5),
                sig = !is.na(padj) & padj < FDR_SIG & dev >= ABS_DEV_SIG,
                direction = dplyr::case_when(!sig ~ "none", alt_frac > 0.5 ~ paste(alt_label, "higher"),
                                             TRUE ~ paste(ref_label, "higher")))
knitr::kable(dplyr::distinct(snp, sample, rho_h0, rho_trim, central_frac, rho_free, rho_corrected, rho, rho_gene, rho_source, rho_gene_source), digits = 4,
             caption = paste("Overdispersion per sample: rho_h0 (p = 0.5 at every site, inflated by real imbalance; diagnostic only, used only as the F1 fallback);",
                             "outbred: rho_trim (H0 fit on the central sites, truncation-corrected) and central_frac (fraction of sites kept), rho = max(rho_trim, RHO_MIN) (used for the SNP tests);",
                             "F1: rho_free (free mean per gene, all SNPs, biased low), rho_corrected (bias-corrected), rho = max(rho_corrected, RHO_MIN) (used for the SNP tests);",
                             "rho_gene (F1 gene LRT): max(bias-corrected free-mean rho of the thinned SNPs, RHO_MIN); fallbacks: unthinned free-mean rho, then rho_h0 (rho_gene_source).",
                             "NA = fewer than 20 sites, sample not tested"))
```

## Gene level

**F1: counts are not summed before testing.** Summing the two strains' reads over a gene's SNPs and applying the per-SNP overdispersion to the total would inflate the variance by `1 + (n - 1) * rho` on the summed depth `n` (about 45 at n = 450), which destroys power. Instead each gene's SNPs share one strain-B fraction `p`, and `bb_gene_lrt` runs a likelihood-ratio test of `p = 0.5` against `p` free on the SNP-level counts with the sample's `rho_gene`; `phat` is the strain-B fraction estimate (`alt_frac` in the table). The summed reads are shown for information only. **Thinned SNPs (F1):** ASEReadCounter counts a read pair at every SNP it covers, so the test uses only the SNPs kept by `thin_snps` in the gene map (`n_snps_used`, against `n_snps` in the sample). F1 gene test: at most one SNP per THIN_BP window (exon coordinates) is used, so no read pair is counted twice (approximate: alternative splicing can bring distant exons into one fragment). The F1 gene table carries a `note` column saying so, and the Rmd prints the same note (Step 14). **Outbred:** `acat` over the SNP p-values (valid under dependence).

```{r gene}
snp_by_gene <- dplyr::inner_join(snp, snp_gene, by = c("contig", "position"))

if (MODE == "f1") {
  n_all <- dplyr::count(snp_by_gene, sample, gene_id, name = "n_snps")   # all SNPs of the gene in the sample (before thinning)
  snp_by_gene_f1 <- dplyr::inner_join(snp, dplyr::filter(snp_gene, thin_keep), by = c("contig", "position"))
  lost <- setdiff(unique(snp_by_gene$sample), unique(snp_by_gene_f1$sample))
  if (length(lost) > 0)
    stop("sample(s) ", paste(lost, collapse = ", "), " have SNPs in GTF genes but none of the SNPs kept by the thinning (one per THIN_BP window ",
         "per gene, chosen on the depth pooled over all samples): every kept SNP is below MIN_DEPTH or missing in these samples, so no F1 ",
         "gene test can be run for them. Check their depth, or drop them (Step 15, \"Dropping a failed sample\")", call. = FALSE)
  gene <- dplyr::bind_rows(lapply(split(snp_by_gene_f1, snp_by_gene_f1$sample), function(d) {
    dplyr::bind_rows(lapply(split(d, d$gene_id), function(g) {
      r <- if (is.na(g$rho_gene[1])) list(p = NA_real_, phat = NA_real_) else bb_gene_lrt(g$alt_n, g$total, g$rho_gene[1])
      data.frame(sample = g$sample[1], gene_id = g$gene_id[1], n_snps_used = nrow(g), ref_n = sum(g$ref_n), alt_n = sum(g$alt_n),
                 total = sum(g$total), rho = g$rho_gene[1], alt_frac = r$phat, p = r$p, stringsAsFactors = FALSE)
    }))
  })) %>%
    dplyr::left_join(n_all, by = c("sample", "gene_id")) %>% dplyr::select(sample, gene_id, n_snps, dplyr::everything()) %>%
    dplyr::mutate(dev = abs(alt_frac - 0.5)) %>%
    dplyr::group_by(sample) %>% dplyr::mutate(padj = p.adjust(p, method = "BH")) %>% dplyr::ungroup() %>%
    dplyr::mutate(sig = !is.na(padj) & padj < FDR_SIG & dev >= ABS_DEV_SIG,
                  direction = dplyr::case_when(!sig ~ "none", alt_frac > 0.5 ~ paste(alt_label, "higher"),
                                               TRUE ~ paste(ref_label, "higher")),
                  note = "thinned: at most one SNP per THIN_BP window (exon coordinates), so no read pair is counted twice; approximate (ignores alternative splicing)")
  f1_note <- "F1 gene test: at most one SNP per THIN_BP window (exon coordinates) is used, so no read pair is counted twice (approximate: alternative splicing can bring distant exons into one fragment)."
  cat("NOTE:", f1_note, "THIN_BP =", THIN_BP, "\n")
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
knitr::kable(head(dplyr::filter(gene, sig), 30), digits = 4,
             caption = paste0("Significant genes (first 30)", if (MODE == "f1") paste0(". ", f1_note) else ""))
```

## Summary

```{r summary}
snp_plot <- dplyr::filter(snp, !is.na(p))      # data behind the SNP figure
gene_plot <- dplyr::filter(gene, if (MODE == "f1") !is.na(p) else !is.na(acat_p))   # data behind the gene figure
summary_tbl <- snp %>% dplyr::group_by(sample) %>%
  dplyr::summarise(rho_used = dplyr::first(rho), rho_robust = dplyr::first(if (MODE == "f1") rho_corrected else rho_trim), rho_h0_naive = dplyr::first(rho_h0), central_frac = dplyr::first(central_frac), rho_free = dplyr::first(rho_free), rho_corrected = dplyr::first(rho_corrected), rho_gene = dplyr::first(rho_gene), rho_source = dplyr::first(rho_source), rho_gene_source = dplyr::first(rho_gene_source), rho_warning = dplyr::first(rho_warning), sites_tested = sum(!is.na(p)), sig_sites = sum(sig),
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
  dplyr::left_join(dplyr::select(bias, sample, mean_ref_frac, ref_frac_before_wasp, ref_frac_after_wasp, ref_bias_flagged = flagged), by = "sample") %>%
  dplyr::mutate(ref_bias_note = bias_note[sample],
                gene_level = if (MODE == "f1") "per-gene LRT on thinned SNP-level counts (at most one SNP per THIN_BP window; shared strain fraction)" else "unphased, no direction (acat of SNP p-values)")
if (MODE != "f1") summary_tbl$rho_gene_source <- NULL   # F1 only (the outbred gene test has no rho_gene)
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
# headline numbers per sample for the summary page (Step 15); every value comes from summary_tbl and bias above
summary_numbers <- summary_tbl %>%
  dplyr::left_join(dplyr::select(bias, sample, filtered_sites = n_sites), by = "sample") %>%
  dplyr::left_join(unique(ck$samples[, c("sample", "condition")]), by = "sample") %>%
  dplyr::transmute(sample, condition, filtered_sites, sig_snps = sig_sites, genes_tested, sig_genes,
                   rho_used, rho_robust, rho_h0_naive, mean_ref_frac, ref_frac_before_wasp, ref_frac_after_wasp,
                   bias_flag = ref_bias_flagged)
write.table(summary_numbers, file.path(RESULTS_DIR, "summary_numbers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", file.path(RESULTS_DIR, "summary_numbers.tsv"), "\n")
sessionInfo()
```
````

Render it with `{RESULTS_DIR}/scripts/run_02_imbalance.sh`: the same script as `run_01_import_qc.sh` (Step 12) with job name `ase_02_imbalance`, log `run_02_imbalance_%j.out`, the Rmd 02 file name and the same `-n 1 --mem=16G -t 1:00:00`. If Rmd 02 later moves to `mclapply`, raise the CPU request together with that change. Submission order and dependencies: Step 15.

Besides the xlsx, Rmd 02 writes `{RESULTS_DIR}/summary_numbers.tsv` (one row per sample, tab-separated, taken from the Summary table above), so the summary page of Step 15 can be written without starting R.

---

## Step 14 — Statistics functions

These base-R functions are pasted verbatim into Rmd 02, Rmd 03, Rmd 04 and Rmd 05 (`{TODAY}_{WD_NAME}_02_imbalance.Rmd`, `..._03_reciprocal.Rmd`, `..._04_differential.Rmd`, `..._05_phaser.Rmd`); the marker lines delimit the code that `ase-pipeline/tests/r/test_ase_stats.R` extracts from this file and unit-tests, so the tested code is the shipped code. They need no packages beyond base R.

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
bb_estimate_rho_trim <- function(x, n, keep = 0.9, max_iter = 100) {   # c(rho, central_frac): outbred rho from the central sites, truncation-corrected
  use <- n > 0 & !is.na(x) & !is.na(n); x <- x[use]; n <- n[use]
  if (length(x) < 20) return(c(rho = NA_real_, central_frac = NA_real_))
  rho <- bb_estimate_rho(x, n)                      # naive all-sites H0 fit (inflated by real imbalance)
  if (rho == 0) return(c(rho = 0, central_frac = 1))
  hi <- qlogis(min(rho, 0.88))                      # upper bound: trimming the tails never makes rho larger than the all-sites fit
  rho <- min(rho, 1e-3)                             # iterate up from a small start: the lowest self-consistent value, least contaminated
  ldb <- function(xx, nn, r) if (r < 1e-8) dbinom(xx, nn, 0.5, log = TRUE) else {
    a <- 0.5 * (1 - r) / r; lchoose(nn, xx) + lbeta(xx + a, nn - xx + a) - lbeta(a, a) }
  s <- rep(seq_along(n), n + 1); k <- sequence(n + 1) - 1; d <- abs(k - n[s] / 2)   # every possible count of every site, folded
  o <- order(s, d); so <- s[o]; ko <- k[o]; do <- d[o]
  hist <- rho; fit <- NULL
  for (it in seq_len(max_iter)) {
    # central region of each site: the counts closest to n/2 holding at least `keep` of the H0 probability under the current rho
    cp <- ave(exp(ldb(ko, n[so], rho)), so, FUN = cumsum)
    t <- as.numeric(tapply(ifelse(cp >= keep - 1e-10, do, Inf), so, min))
    inside <- abs(x - n / 2) <= t                     # sites whose observed count lies in its central region
    reg <- d <= t[s] & inside[s]
    nll <- function(par) {                            # truncated likelihood: each kept site's density divided by the H0 mass of its region
      r <- plogis(par)
      -sum(ldb(x[inside], n[inside], r)) + sum(log(rowsum(exp(ldb(k[reg], n[s[reg]], r)), s[reg])))
    }
    fit <- optimize(nll, interval = c(-12, hi))
    rn <- if (fit$minimum < -12 + 0.05 || plogis(fit$minimum) < 1e-4) 0 else plogis(fit$minimum)
    tol <- 1e-3 * max(rn, rho) + 1e-6
    if (abs(rn - rho) <= tol) { rho <- rn; break }
    cyc <- which(abs(hist - rn) <= tol)               # discrete regions can make the iteration cycle: keep the largest value of the cycle
    if (length(cyc) > 0) { rho <- max(hist[min(cyc):length(hist)], rn); break }
    rho <- rn; hist <- c(hist, rn)
    if (it == max_iter) warning("trimmed rho estimate did not converge in ", max_iter, " iterations; the last value is used")
  }
  c(rho = rho, central_frac = mean(inside))
}
chrom_class <- function(contig) {   # "autosome", "X", "Y" or "MT" (Ensembl and chr-prefixed names); Rmd 03/04 test autosomes only
  s <- toupper(sub("^chr", "", as.character(contig), ignore.case = TRUE))
  ifelse(s == "X", "X", ifelse(s == "Y", "Y", ifelse(s %in% c("MT", "M"), "MT", "autosome")))
}
thin_snps <- function(pos, depth, window) {   # logical keep: deepest SNP first (ties: lower position); no two kept SNPs closer than window
  bad <- is.na(pos) | is.na(depth)
  if (any(bad)) stop("thin_snps: ", sum(bad), " SNP(s) with a missing position or depth; remove them before thinning")
  keep <- rep(FALSE, length(pos))
  for (i in order(-depth, pos)) if (!any(keep & abs(pos - pos[i]) < window)) keep[i] <- TRUE
  keep
}
bb_ll_rows <- function(x, n, p, rho) {   # per-row beta-binomial log-likelihood, p may differ by row; rho < 1e-8 -> binomial
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
  if (rho < 1e-8) return(dbinom(x, n, p, log = TRUE))
  a <- p * (1 - rho) / rho; b <- (1 - p) * (1 - rho) / rho
  lchoose(n, x) + lbeta(x + a, n - x + b) - lbeta(a, b)
}
bb_glm_fit <- function(y, n, X, rho, bound = 15) {   # beta-binomial GLM, logit link, dispersion rho FIXED; coefficients boxed to [-bound, bound]
  X <- as.matrix(X)
  if (ncol(X) == 0) return(list(beta = numeric(0), loglik = sum(bb_ll_rows(y, n, rep(0.5, length(y)), rho)), conv = 0L, at_bound = FALSE))
  nll <- function(b) -sum(bb_ll_rows(y, n, plogis(drop(X %*% b)), rho))
  gr <- function(b) {   # analytic gradient (finite differences stall near the bounds of monoallelic genes)
    p <- pmin(pmax(plogis(drop(X %*% b)), 1e-10), 1 - 1e-10)
    s <- if (rho < 1e-8) y - n * p else {
      k <- (1 - rho) / rho
      k * p * (1 - p) * (digamma(y + p * k) - digamma(n - y + (1 - p) * k) - digamma(p * k) + digamma((1 - p) * k))
    }
    -drop(crossprod(X, s))
  }
  st <- tryCatch(as.numeric(qr.solve(X, qlogis((y + 0.5) / (n + 1)))), error = function(e) rep(0, ncol(X)))
  st[!is.finite(st)] <- 0; st <- pmin(pmax(st, -bound + 1), bound - 1)
  o <- optim(st, nll, gr, method = "L-BFGS-B", lower = -bound, upper = bound)
  conv <- o$convergence
  if (conv != 0) {   # L-BFGS-B often stops with code 52 (line search) AT the optimum, where the log-likelihood is flat to rounding:
    g <- gr(o$par)   # accept the fit when the projected gradient vanishes (KKT: at a bound only an outward gradient remains)
    g[(o$par <= -bound + 1e-8 & g >= 0) | (o$par >= bound - 1e-8 & g <= 0)] <- 0
    if (all(is.finite(g)) && max(abs(g)) < 1e-3) conv <- 0L
  }
  list(beta = setNames(o$par, colnames(X)), loglik = -o$value, conv = conv, at_bound = any(abs(o$par) > bound - 1e-3))
}
bb_glm_lrt <- function(y, n, X, drop_cols, rho, bound = 15) {   # LRT of H0: coefficients of drop_cols = 0, rho fixed; df = length(drop_cols)
  X <- as.matrix(X)
  miss <- setdiff(drop_cols, colnames(X))   # never a silent stat 0 / p = 1 for a column the design does not have
  if (length(miss) > 0) stop("bb_glm_lrt: tested column(s) not in the design: ", paste(miss, collapse = ", "))
  f1 <- bb_glm_fit(y, n, X, rho, bound)
  f0 <- bb_glm_fit(y, n, X[, setdiff(colnames(X), drop_cols), drop = FALSE], rho, bound)
  stat <- max(0, 2 * (f1$loglik - f0$loglik))
  list(beta = f1$beta, stat = stat, df = length(drop_cols), p = pchisq(stat, length(drop_cols), lower.tail = FALSE),
       conv = max(f1$conv, f0$conv), at_bound = f1$at_bound)
}
bb_moment_phi <- function(units, max_iter = 25, tol = 1e-6, bound = 15) {
  # pooled moment (Williams-type) estimate of the beta-binomial dispersion over units (each list(y, n, X), rows with n > 0):
  # solves sum(Pearson residual^2 / (1 + (n - 1) phi)) = residual df, the means refitted at each phi. Dividing by the residual
  # df (rows minus coefficients) is what avoids the Neyman-Scott shrinkage of a maximum-likelihood fit with free means.
  phi <- 0
  for (it in seq_len(max_iter)) {
    parts <- lapply(units, function(u) tryCatch({
      f <- bb_glm_fit(u$y, u$n, u$X, phi, bound)
      p <- pmin(pmax(plogis(drop(as.matrix(u$X) %*% f$beta)), 1e-6), 1 - 1e-6)
      r2 <- (u$y - u$n * p)^2 / (u$n * p * (1 - p)); dfu <- as.numeric(length(u$y) - qr(as.matrix(u$X))$rank)
      if (dfu > 0 && all(is.finite(r2))) list(r2 = r2, n = u$n, df = dfu) else NULL
    }, error = function(e) NULL))
    parts <- parts[!vapply(parts, is.null, logical(1))]
    if (length(parts) == 0) return(NA_real_)
    r2 <- unlist(lapply(parts, `[[`, "r2")); nn <- unlist(lapply(parts, `[[`, "n")); df <- sum(vapply(parts, `[[`, numeric(1), "df"))
    g <- function(ph) sum(r2 / (1 + (nn - 1) * ph)) - df
    new <- if (g(0) <= 0) 0 else if (g(0.99) > 0) 0.99 else uniroot(g, c(0, 0.99), tol = 1e-12)$root
    if (abs(new - phi) < tol) return(new)
    phi <- new
  }
  warning("dispersion estimate did not converge in ", max_iter, " iterations; the last value is used")
  phi
}
ase_glm_test <- function(units, tests, rho_min, min_df_unit = 2, bound = 15) {
  # units: named list of list(y, n, X), one row per sample (rows with n > 0 only), X with named columns
  # tests: named list, test name -> the X columns set to 0 under H0 (likelihood-ratio test, df = number of columns)
  # dispersion: phi_common = pooled moment estimate over every estimable unit; per unit phi_used = max(phi_common, rho_min,
  # phi_unit when the unit has at least min_df_unit residual df): the larger of the shared and the unit's own estimate
  nr <- vapply(units, function(u) length(u$y), integer(1))
  rk <- vapply(units, function(u) if (length(u$y) == 0) 0L else qr(as.matrix(u$X))$rank, integer(1))
  est <- nr > 0 & rk == vapply(units, function(u) ncol(as.matrix(u$X)), integer(1)) & nr > rk
  common <- if (any(est)) bb_moment_phi(units[est], bound = bound) else NA_real_
  cn <- unique(unlist(lapply(units, function(u) colnames(u$X))))
  one <- function(i) {
    u <- units[[i]]
    row <- c(list(unit = names(units)[i], n_rows = nr[i], df_resid = nr[i] - rk[i], depth = sum(u$n),
                  phi_common = common, phi_unit = NA_real_, phi_used = NA_real_, status = "not_estimable"),
             setNames(as.list(rep(NA_real_, length(cn))), paste0("beta_", cn)),
             setNames(as.list(rep(NA_real_, length(tests))), paste0("stat_", names(tests))),
             setNames(as.list(rep(NA_real_, length(tests))), paste0("p_", names(tests))))
    if (est[i] && !is.na(common)) {
      res <- tryCatch({
        pu <- if (nr[i] - rk[i] >= min_df_unit) bb_moment_phi(list(u), bound = bound) else NA_real_
        ph <- max(c(common, rho_min, pu), na.rm = TRUE)
        list(pu = pu, ph = ph, tt = lapply(tests, function(cols) bb_glm_lrt(u$y, u$n, u$X, cols, ph, bound)))
      }, error = function(e) NULL)
      if (is.null(res)) row$status <- "fit_error" else {
        row$phi_unit <- res$pu; row$phi_used <- res$ph
        b <- res$tt[[1]]$beta; row[paste0("beta_", names(b))] <- as.list(unname(b))
        conv <- max(vapply(res$tt, function(t) t$conv, numeric(1)))
        for (nm in names(tests)) {
          row[[paste0("stat_", nm)]] <- res$tt[[nm]]$stat
          row[[paste0("p_", nm)]] <- if (conv == 0) res$tt[[nm]]$p else NA_real_
        }
        row$status <- if (conv != 0) "not_converged" else if (any(vapply(res$tt, function(t) t$at_bound, logical(1)))) "ok_at_bound" else "ok"
      }
    }
    as.data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, lapply(seq_along(units), one)); rownames(out) <- NULL; out
}
bb_pair_phi_trim <- function(y, n, cell, keep = 0.9, max_iter = 50) {
  # c(phi, phi_all, central_frac): pair dispersion of ONE individual, robust to real condition changes.
  # cell = SNP id; the rows of a cell are that SNP's samples (both conditions). Per cell, under H0 (one REF fraction
  # for all its rows), the Williams statistic X2 = sum((y - n p)^2 / (n p (1 - p) (1 + (n - 1) phi))), p the
  # weighted pooled fraction, is about chi-square with df = rows - 1. phi_all solves sum(X2) = sum(df) over all
  # cells (the plain moment estimate, inflated by every SNP that really changes). The trimmed estimate keeps the
  # central cells (pchisq(X2, df) <= keep, so the most-changed cells are set aside) and solves
  # sum(X2) = sum(E[chi2_df | chi2_df <= qchisq(keep, df)]) over them: the truncation-corrected moment equation.
  # Kept set and phi are iterated DOWN from phi_all to a fixed point (if the iteration cycles, the largest value of the
  # cycle is kept); the result is capped at phi_all (trimming never makes the data more overdispersed).
  ok <- is.finite(y) & is.finite(n) & n > 0; y <- y[ok]; n <- n[ok]; cell <- as.character(cell[ok])
  rows <- table(cell); cell <- factor(cell, levels = names(rows)[rows >= 2]); use <- !is.na(cell)
  y <- y[use]; n <- n[use]; cell <- droplevels(cell)
  if (nlevels(cell) < 20) return(c(phi = NA_real_, phi_all = NA_real_, central_frac = NA_real_))
  df <- as.numeric(table(cell)) - 1
  x2 <- function(ph) {
    w <- n / (1 + (n - 1) * ph)
    p <- (rowsum(w * y / n, cell) / rowsum(w, cell))[as.integer(cell)]
    p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
    as.numeric(rowsum((y - n * p)^2 / (n * p * (1 - p) * (1 + (n - 1) * ph)), cell))
  }
  solve_phi <- function(k, target, hi) {   # phi in [0, hi] with sum(x2[k]) = target
    g <- function(ph) sum(x2(ph)[k]) - target
    if (g(0) <= 0) 0 else if (g(hi) >= 0) hi else uniroot(g, c(0, hi), tol = 1e-10)$root
  }
  phi_all <- solve_phi(rep(TRUE, length(df)), sum(df), 0.99)
  if (phi_all == 0) return(c(phi = 0, phi_all = 0, central_frac = 1))
  m <- df * pchisq(qchisq(keep, df), df + 2) / keep   # E[chi2_df | chi2_df <= its keep quantile]
  phi <- phi_all; hist <- phi; k <- rep(TRUE, length(df))
  for (it in seq_len(max_iter)) {
    k <- pchisq(x2(phi), df) <= keep
    pn <- solve_phi(k, sum(m[k]), phi_all)
    tol <- 1e-4 * max(pn, phi) + 1e-8
    if (abs(pn - phi) <= tol) { phi <- pn; break }
    cyc <- which(abs(hist - pn) <= tol)
    if (length(cyc) > 0) { phi <- max(hist[min(cyc):length(hist)], pn); break }
    phi <- pn; hist <- c(hist, pn)
    if (it == max_iter) warning("trimmed pair dispersion did not converge in ", max_iter, " iterations; the last value is used")
  }
  c(phi = phi, phi_all = phi_all, central_frac = mean(k))
}
ase_paired_test <- function(d, rho_min, min_individuals = 2, bound = 15) {
  # d: data.frame(snp, individual, cond, y, n); cond 0 = reference condition, 1 = tested condition; y = REF count, n = REF + ALT
  # 1. per individual: pair dispersion phi_i = bb_pair_phi_trim (trimmed moment estimate, one REF fraction per SNP; the
  #    SNPs that change most between conditions are set aside, so real changes do not inflate it), floored at rho_min
  # 2. per SNP and individual: LRT of one shared REF fraction against one per condition (1 df) with phi_i
  # 3. per SNP: the statistics summed over the informative individuals (df = their number). Direction-free, because
  #    which allele carries a regulatory variant differs between individuals (a shared slope would cancel real changes)
  empty <- data.frame(snp = character(0), n_individuals = integer(0), n_failed = integer(0), stat = numeric(0), df = integer(0),
                      p = numeric(0), mean_delta = numeric(0), max_abs_delta = numeric(0), n_up = integer(0), n_down = integer(0),
                      status = character(0), stringsAsFactors = FALSE)
  d <- d[is.finite(d$y) & is.finite(d$n) & d$n > 0, c("snp", "individual", "cond", "y", "n"), drop = FALSE]
  key <- paste(d$individual, d$snp, sep = "\t")
  paired <- tapply(d$cond, key, function(v) any(v == 0) & any(v == 1))
  d <- d[as.logical(paired[key]), , drop = FALSE]
  if (nrow(d) == 0) return(list(snp = empty, phi = data.frame(individual = character(0), phi_pair = numeric(0))))
  key <- paste(d$individual, d$snp, sep = "\t")
  cells <- split(seq_len(nrow(d)), key); first <- vapply(cells, `[`, integer(1), 1)
  cell_ind <- d$individual[first]; cell_snp <- d$snp[first]
  inds <- sort(unique(cell_ind))
  phi <- vapply(inds, function(ind) {   # NA (individual not tested) below 20 paired SNPs
    i <- d$individual == ind
    ph <- bb_pair_phi_trim(d$y[i], d$n[i], d$snp[i])[["phi"]]
    if (is.na(ph)) NA_real_ else max(ph, rho_min)
  }, numeric(1))
  names(phi) <- inds
  cr <- do.call(rbind, lapply(seq_along(cells), function(j) {
    i <- cells[[j]]; ph <- phi[[cell_ind[j]]]
    r <- if (is.na(ph)) NULL else tryCatch(bb_glm_lrt(d$y[i], d$n[i], cbind("(Intercept)" = 1, cond = d$cond[i]), "cond", ph, bound),
                                           error = function(e) NULL)
    if (is.null(r) || r$conv != 0) return(data.frame(snp = cell_snp[j], stat = NA_real_, delta = NA_real_, stringsAsFactors = FALSE))
    data.frame(snp = cell_snp[j], stat = r$stat, delta = plogis(sum(r$beta)) - plogis(r$beta[[1]]), stringsAsFactors = FALSE)
  }))
  snp <- do.call(rbind, lapply(split(cr, cr$snp), function(s) {
    k <- sum(is.finite(s$stat)); so <- s[is.finite(s$stat), , drop = FALSE]; enough <- k >= min_individuals
    data.frame(snp = s$snp[1], n_individuals = k, n_failed = nrow(s) - k,
               stat = if (enough) sum(so$stat) else NA_real_, df = k,
               p = if (enough) pchisq(sum(so$stat), k, lower.tail = FALSE) else NA_real_,
               mean_delta = if (k > 0) mean(so$delta) else NA_real_, max_abs_delta = if (k > 0) max(abs(so$delta)) else NA_real_,
               n_up = sum(so$delta > 0), n_down = sum(so$delta < 0),
               status = if (enough) "ok" else "too_few_individuals", stringsAsFactors = FALSE)
  }))
  rownames(snp) <- NULL
  list(snp = snp, phi = data.frame(individual = inds, phi_pair = unname(phi), stringsAsFactors = FALSE))
}
hap_gene_test <- function(a, b, sample, rho_min, min_genes = 20) {
  # phASER gene-level haplotype counts, one row per sample and gene (a = aCount, b = bCount of phaser_gene_ae). Symmetric in a
  # and b: without phased genotypes the labels A and B are arbitrary per gene and sample, so only |a / (a + b) - 0.5| carries
  # information, and every quantity below is unchanged when a and b are swapped in any row.
  # Dispersion: rho_own = bb_estimate_rho_trim over the sample's usable genes (folded central-region H0 fit, truncation-
  # corrected, so genes with real imbalance inflate it far less than an all-genes fit; NA below min_genes genes);
  # rho_cohort = the median of rho_own over the samples; per sample rho_used = max(rho_cohort, rho_own, rho_min).
  # Test per row: bb_pvalue against 0.5 (exact, two-sided). Rows with a or b missing, infinite, negative or non-integer,
  # or a + b = 0, are not used (not in n_genes, not in the fit) and get p = NA, as does every row when no sample has
  # min_genes usable rows. MIN_DEPTH is not applied here: the caller sets rows below MIN_DEPTH to NA first. A missing
  # sample label or vectors of different lengths stop.
  if (length(a) != length(b) || length(b) != length(sample)) stop("hap_gene_test: a, b and sample must have the same length")
  if (anyNA(sample)) stop("hap_gene_test: sample labels must not be NA")
  n <- a + b
  use <- is.finite(a) & is.finite(b) & a >= 0 & b >= 0 & n > 0 & a == round(a) & b == round(b)
  sample <- as.character(sample); smp <- sort(unique(sample))
  own <- vapply(smp, function(s) {
    i <- use & sample == s
    if (sum(i) >= min_genes) bb_estimate_rho_trim(a[i], n[i]) else c(rho = NA_real_, central_frac = NA_real_)
  }, c(rho = 0, central_frac = 0))
  coh <- if (all(is.na(own["rho", ]))) NA_real_ else median(own["rho", ], na.rm = TRUE)
  ru <- if (is.na(coh)) rep(NA_real_, length(smp)) else pmax(coh, ifelse(is.na(own["rho", ]), 0, own["rho", ]), rho_min)
  per <- data.frame(sample = smp, n_genes = unname(vapply(smp, function(s) sum(use & sample == s), integer(1))),
                    rho_own = unname(own["rho", ]), central_frac = unname(own["central_frac", ]), rho_cohort = coh,
                    rho_used = unname(ru), stringsAsFactors = FALSE)
  rr <- per$rho_used[match(sample, per$sample)]
  p <- rep(NA_real_, length(n))
  for (j in which(use & !is.na(rr))) p[j] <- bb_pvalue(a[j], n[j], rr[j])
  list(p = p, rho_used = rr, samples = per)
}
# --- ase-stats-end
```

- `bb_pvalue` is the two-sided exact beta-binomial test against `p0`: the p-value is the total probability of all outcomes that are no more likely than the observed count. With `rho = 0` it is the exact binomial test.
- **`bb_estimate_rho` is the naive H0 fit (p = 0.5) on all filtered sites.** When real imbalance exists, the between-site spread it causes is absorbed into `rho`, which inflates it; on the synthetic acceptance data (binomial truth, about a third of the sites imbalanced) it gave 0.08-0.10, and testing with it detected 0 of 16 planted outbred genes. It is therefore only a diagnostic (`rho_h0_naive`), the upper bound of the trimmed fit, and the F1 fallback. It is estimated by maximum likelihood on the logistic scale (`optimize`), and a value below 1e-4 is reported as 0 (pure binomial). Fewer than 20 sites gives `NA` and the sample is not tested.
- **Outbred `rho`: `bb_estimate_rho_trim` (central sites, truncation-corrected).** Unphased data give no free mean per gene (the SNP orientation is unknown), so imbalanced sites are handled by trimming. Each site's counts are folded around n/2; its central region is the set of counts closest to n/2 that holds `keep` = 90 percent of the H0 beta-binomial probability under the current `rho`. Sites whose observed count falls outside their region are set aside, and `rho` is re-fitted on the kept sites with a **truncated likelihood**: each kept site's H0 density is divided by the H0 probability of its region. That division is the correction for the truncation: without it, cutting the tails would shrink `rho` and make the tests anti-conservative. Region and `rho` are iterated upward from a small start (0.001) to a fixed point, the lowest self-consistent value and so the one least contaminated by imbalanced sites (if the discrete regions make the iteration cycle, the largest value of the cycle is kept), and `rho_trim` is bounded above by the all-sites fit (trimming tails cannot make the data more overdispersed, and an unbounded fit on very few sites can run to a U-shaped distribution under which central counts look significant). `central_frac` is the fraction of sites kept. Rmd 02 tests with `max(rho_trim, RHO_MIN)`. Unit-test results: 10 seeds x 300 genes of 1-4 SNPs, unphased, depth 30-200, `bb_pvalue` per SNP and `acat` per gene. With a third of the genes imbalanced (0.85/0.70/0.30) on binomial data, as in the acceptance data, `rho_trim` is 0.0007 (naive fit 0.093), the null size is 0.004 (SNP) / 0.003 (gene) and the power is 1.00 for 0.85 genes and 0.88 for 0.70/0.30 genes, against 0 with the naive fit; with a true `rho` of 0.01 the power is 1.00 / 0.74 (null gene size 0.016). Without imbalance `rho_trim` recovers the truth (0.0196 for 0.02, 0.0497 for 0.05) and the pooled null gene size is 0.040 / 0.053 (worst seed 0.083). On samples the size of the synthetic acceptance data (about 30 sites, 4 of 12 genes planted, 100 seeds) 3.3 of 4 planted genes are detected with a null false-positive rate of 0.003; in about one sample in ten the fit cannot separate the few imbalanced sites and stays near the naive value (conservative). On the acceptance count tables themselves Rmd 02 now detects 10 of 16 planted gene x sample tests (before: 0) with 0 of 28 null false positives. The remaining limitation: with real overdispersion and many imbalanced sites, imbalanced sites that fall inside the central regions keep `rho_trim` above the truth (0.039 for a true 0.02 with a third of the genes imbalanced), which lowers the power for moderate imbalance (0.23 for 0.70/0.30 genes) but adds no false positives.
- `acat` is the Cauchy combination test with equal weights; it is valid for correlated p-values, which is why it is used for the SNPs of one gene. It carries no direction. P-values are capped at 0.99 before the transform: `bb_pvalue` returns exactly 1 at the modal count, and `tan((0.5 - 1) * pi)` is about -1.6e16, which would swamp any signal and make the combined p-value 1. Values are also clipped at 1e-15 from below.
- The unit tests check the null size and uniformity of the p-values, power, recovery of `rho` (including `rho` near 0), the edge cases (`n = 0`, `x = 0`, `x = n`), `acat`, and the full outbred path (trimmed `rho`, floor, per-SNP test, ACAT per gene) on simulated unphased data with planted 0.85/0.70/0.30 genes. Run them from the repository root with `singularity exec --bind /net/bmc-lab3 <bulkrnaseq sif> Rscript ase-pipeline/tests/r/test_ase_stats.R ase-pipeline/ase-pipeline.md` inside an `sbatch -p bcc` job.
- **F1 gene level uses `bb_estimate_rho_gene` and `bb_gene_lrt`, not summed counts.** `bb_estimate_rho` (H0, p = 0.5 everywhere) treats the between-SNP spread caused by true imbalance as overdispersion, so with many imbalanced genes it is inflated (in the unit test about four times the true value), and applying it to counts summed over a gene multiplies the variance by about `1 + (n - 1) * rho`. `bb_estimate_rho_gene` maximises the likelihood with a free mean for every gene (genes with at least 2 SNPs; `NA` below 5 such genes) and returns `c(free, corrected)`. The free estimate is biased low: a free mean per gene shrinks the whole variance factor `1 + (n - 1) * rho` by `(k - 1) / k` (Neyman-Scott), which for 4-SNP genes gives 0.0115 for a true 0.02. The corrected value undoes that shrinkage, `rho_corrected = cf * rho_free + (cf - 1) / m` with `cf = sum(k) / sum(k - 1)` and `m` the mean of `n - 1` (0.021 in the same test). Rmd 02 then uses `rho_gene = max(rho_corrected, RHO_MIN)` for the gene tests (`RHO_MIN`, default 0.01, from Step 8), so a noisy or zero estimate can never make them anti-conservative. `bb_gene_lrt` tests `p = 0.5` against a shared free `p` per gene on the SNP-level counts and returns `phat`. Outbred mode uses `bb_estimate_rho_trim`, `bb_pvalue` and `acat`.
- **Thinned SNPs (F1 gene LRT and `bb_estimate_rho_gene`).** ASEReadCounter counts a fragment at every SNP it overlaps, so when several SNPs of a gene lie within one fragment (common in divergent F1 haplotype blocks) their counts share read pairs; a likelihood that treats them as independent draws lets the LRT statistic grow with the number of SNPs per fragment, and the free-mean `rho` estimate cannot compensate (SNPs that share reads agree with each other). The Rmd 02 F1 gene LRT therefore uses only the SNPs kept by `thin_snps`: at most one SNP per `THIN_BP` window of each gene in exon coordinates (exon union of the gene), the deepest first (ties: lower position), so no read pair is counted twice; this is approximate, because alternative splicing can bring distant exons into one fragment. The gene table reports `n_snps` (the gene's SNPs in the sample) and `n_snps_used` (the thinned SNPs the test used). **Fallback chain of `rho_gene`.** The gene LRT has its own `rho_gene`: (1) `bb_estimate_rho_gene` on the thinned SNPs, bias-corrected and floored at `RHO_MIN`; (2) if that is `NA` (fewer than 5 genes with 2 or more thinned SNPs, as on the synthetic acceptance data, where 2 of 12 genes keep 2 SNPs), the unthinned free-mean value that the SNP tests use, `max(rho_corrected, RHO_MIN)`; (3) only if both are `NA`, `max(rho_h0, RHO_MIN)`; the Summary column `rho_gene_source` records which. Measured in the fragment-level simulation of `test_ase_stats_stage2.R` (see "SNPs sharing read pairs" below; 10 seeds x 1500 null genes, 6 samples, 30 percent SNP-dense genes), the size of the Rmd 02 F1 gene LRT at a nominal 0.05 was 0.1450 unthinned and 0.0497 thinned (genes with 2 or more SNPs left; worst seed 0.0551); in the SNP-dense genes it was 0.2084 unthinned, and thinned these genes keep a single SNP, so they have no multi-SNP gene test (a gene with one SNP left is tested on that SNP alone). The same thinned test with the unthinned `rho` gave 0.0497 (worst seed 0.0551). In a fallback scenario (10 seeds x 300 null genes, 6 samples, all but 3 genes SNP-dense, so step (2) is taken in every sample) the size was 0.0467 (worst seed 0.0556), 0.0462 in the SNP-dense genes (worst seed 0.0550). In both simulations the free-mean estimate is 0 and `rho_gene` sits at the `RHO_MIN` floor (the simulated animal-level spread, 0.005, is below it), so they do not show how an unthinned estimate above the floor would behave: the fallback with the unthinned dispersion has not been verified when the true overdispersion is above RHO_MIN; because SNPs that share reads push that estimate low, the gene test may be anti-conservative in that case. A correct gene test costs power in genes whose SNPs are close together: on the synthetic acceptance F1 data most genes keep a single SNP, the thinned test finds 21 of 24 planted gene x sample tests with 0 of 42 null false positives (the three misses are 0.70 genes whose single kept SNP gives padj 0.051-0.19), and the Stage 1 figure of 24 of 24 planted gene x sample tests no longer applies (it relied on counting shared read pairs several times). Stated property of the thinning: to choose the kept SNP, the depth is pooled over all samples of the project, so every sample uses the same choice of SNPs, a kept SNP that is below `MIN_DEPTH` in one sample leaves that sample's gene test with fewer SNPs (or none), and adding samples or conditions to the project can change which SNP is kept in a dense gene (and so the gene p-values of the samples already analysed). The SNP-level tests are unchanged: the same `bb_pvalue` per SNP with the Stage 1 `rho` (free-mean estimate on all SNPs, fallback `rho_h0`), and so is the outbred path (per-SNP tests and the ACAT gene test, which is valid under dependence).
- Boundary handling: the `rho` estimators return exactly 0 when the optimum sits at the lower boundary (or below 1e-4), `bb_estimate_rho` and `bb_estimate_rho_gene` emit a warning when it sits at the upper boundary, and `bb_estimate_rho_trim` warns if it does not converge; Rmd 02 prints such warnings per sample and writes them to the Summary (`rho_warning`). `bb_pvalue` returns `NA` for `x` below 0, above `n`, non-integer or missing.
- **Stage 2 models (Rmd 03 and 04): `ase_glm_test`.** Each gene (or SNP) is one unit, with one row per sample. For F1 genes, a row is the sum of strain-A and total counts over the gene's SNPs after `thin_snps` (one SNP per `THIN_BP` window in exon coordinates, deepest first), so a read pair is not counted twice. The beta-binomial GLM (`bb_glm_fit`: logit link, dispersion fixed, coefficients boxed to ±15 so that monoallelic genes give a finite likelihood-ratio statistic) is tested by `bb_glm_lrt`. L-BFGS-B often stops with code 52 (line search) exactly at the optimum; such a fit is accepted when the projected gradient of the log-likelihood is below 1e-3 (in a development run, all 124 code-52 stops had at most 9.3e-6; the regression test in `test_ase_stats_stage2.R` reproduces one such stop), and a fit that truly stopped early stays `not_converged` with p `NA`. The dispersion is `phi_common` from `bb_moment_phi`: the spread between replicate samples, divided by the residual df, so there is no Neyman-Scott shrinkage (it recovers 0.0051 / 0.0207 / 0.0514 for 0.005 / 0.02 / 0.05). Each unit uses `phi_used = max(phi_common, RHO_MIN, phi_unit)`, and `phi_unit` counts only with at least 2 residual df. Measured (all Stage 2 numbers below come from one run of `ase-pipeline/tests/r/test_ase_stats_stage2.R`, except where a bullet names another source), 10 seeds x 2000 null genes per design, true dispersion 0.02: the reciprocal null sizes (strain / parent of origin) are 0.0370 / 0.0367 (3+3), 0.0362 / 0.0346 (2+2), 0.0356 / 0.0341 (3+2) and 0.0380 / 0.0365 (4+2); leakage (strain genes in the parent-of-origin test and the reverse) is at most 0.0387; with heterogeneous dispersion (lognormal sd 0.7) 0.0309 / 0.0291; at a true dispersion of 0.005 (below `RHO_MIN`) 0.0132 / 0.0129. The tests are therefore conservative (about 0.036 at a nominal 0.05, because of the maximum rule and the floor). Power for 3+3 is 1.000 (strain 0.7) / 1.000 (maternal 0.8), and signs are correct in every call. The F1 differential null size is 0.0369 (3 vs 3), 0.0363 (2 vs 2), 0.0365 (3 vs 2) and 0.0372 (unbalanced directions with the cross-direction term); imprinted genes with no condition effect give 0.0420 / 0.0420 / 0.0430 / 0.0383. Power for a 0.5 -> 0.7 change is 0.7813 (3 vs 3), 0.6853 (3 vs 2) and 0.5983 (2 vs 2). Without the cross-direction term, the imprinted-gene size in the unbalanced design is 0.0000 (conservative, not inflated, in this design); the term is kept so that the condition effect is not mixed with parent of origin.
- **Outbred differential: `ase_paired_test`.** Each individual gets its own pair dispersion from `bb_pair_phi_trim`, floored at `RHO_MIN`. Under H0 each SNP has one REF fraction for all its samples, so the Williams statistic of a SNP is about chi-square with df = rows - 1. The plain moment estimate (`phi_all`, sum of the statistics = sum of the df) counts every SNP that really changes between conditions as noise: with 10 percent of the SNPs changed 0.5 -> 0.8 it gives 0.0365 (phase-heterogeneous) / 0.0394 (consistent) for a true 0.02, and a power of 0.7625 / 0.7977. The trimmed estimate keeps the central SNPs (`pchisq(X2, df) <= keep`, `keep` = 0.9, so the most-changed SNPs are set aside), solves the truncation-corrected moment equation on them (each kept statistic matched to `E[chi2_df | chi2_df <= qchisq(keep, df)]`), iterates down from `phi_all` to a fixed point and is capped at `phi_all`; fewer than 20 paired SNPs give `NA` and the individual is not tested. It recovers 0.0097 / 0.0201 / 0.0507 for 0.01 / 0.02 / 0.05 without changes, and 0.0248 for 0.02 with 10 percent of the SNPs changed. A 1-df LRT per SNP and individual uses that dispersion; the statistics are summed over individuals (df = number of individuals), which is direction-free because the allele carrying a regulatory variant differs between individuals. Genes use `acat`. A SNP is tested only when at least 2 individuals are heterozygous and paired for it; the others are reported as `too_few_individuals` with p `NA`. With 60 percent of individuals heterozygous per SNP, the tested fraction of the reported SNPs is 0.4244 (I = 2), 0.6961 (I = 3), 0.8434 (I = 4) and 0.9646 (I = 6). The paired SNP sizes and power are over the tested SNPs only. Measured, 10 seeds x 3000 SNPs, true dispersion 0.02: the SNP / gene null sizes are 0.0485 / 0.0490 (I = 2), 0.0511 / 0.0505 (I = 3), 0.0515 / 0.0526 (I = 4) and 0.0490 / 0.0514 (I = 6), so the paired test is nominal (not conservative); at a true dispersion of 0.005 (below `RHO_MIN`) they are at most 0.0205 (SNP) / 0.0191 (gene). Power for a 0.5 -> 0.8 change in 10 percent of the SNPs with 4 individuals is 0.8653 (phase-heterogeneous) / 0.9174 (consistent), against an oracle power of 0.9022 / 0.9424 with the dispersion fixed at its true value (the ceiling of this design: 4 individuals, 60 percent heterozygous, depth 30-200), while the unchanged SNPs keep a size of 0.0287 / 0.0304. When 30 percent of the SNPs change, the trimmed dispersion is still inflated (0.0454 for 0.02) and power drops to 0.6698 against an oracle 0.9163 (untrimmed: 0.3675), which is conservative (no false positives added).
- **X chromosome in F1 mode.** The functions above do not know chromosomes: Rmd 03 and Rmd 04 leave X, Y and MT out with `chrom_class`, but Rmd 01 and Rmd 02 do not set X aside (they keep it in every table and in the reference-bias mean). Male F1 animals carry one X, so X genes can look imbalanced and X sites can trip the reference-bias flag; read X rows with care (Step 13).
- **SNPs sharing read pairs.** In a fragment-level simulation (ASEReadCounter-style counting, 30 percent SNP-dense genes), the thinned gene tests have null sizes of 0.0314 (strain) / 0.0288 (parent of origin) (dense genes 0.0331), against 0.0354 without thinning (dense genes 0.0578). The Stage 1 Rmd 02 F1 gene LRT has a size of 0.1450 unthinned (0.2084 in SNP-dense genes): counting shared read pairs at every SNP makes it anti-conservative. Thinned, its size is 0.0497; that figure excludes the SNP-dense genes, which keep a single SNP after thinning and so have no multi-SNP gene test.
- **phASER gene test (Rmd 05): `hap_gene_test`.**
  - **Input:** one row per sample and gene with phASER's haplotype counts (`aCount`, `bCount` of `phaser_gene_ae`). phASER counts a read pair once per haplotype block, so a block's counts are not inflated by SNPs that share reads.
  - **One block per gene without genome-wide phasing:** with unphased genotypes `phaser_gene_ae` does not add up a gene's haplotype blocks; it keeps only the most-covered block (or single variant) of the gene (`phaser_gene_ae.py`, lines 136-164 at the pinned commit). The counts therefore cover the SNPs that read pairs link into that one block, not every SNP of the gene. Only with phased genotypes (phASER's genome-wide phasing option) are the blocks of a gene summed.
  - **Rows the function uses:** `hap_gene_test` does not apply `MIN_DEPTH`; Rmd 05 sets rows below `MIN_DEPTH` to `NA` before the call. Rows with a missing, infinite, negative or non-integer count, or no reads, get p `NA` and count neither in `n_genes` nor in the dispersion fit. A missing sample label, or `a`, `b` and `sample` of different lengths, stop with an error.
  - **Test:** each row is tested against 0.5 with `bb_pvalue`. The test is symmetric in the two haplotypes: without phased genotypes the labels are arbitrary per gene and sample. Swapping the labels in any row changes no p-value and no dispersion; the unit tests check this.
  - **Dispersion:** `rho_own` is the trimmed central-region fit (`bb_estimate_rho_trim`) over the sample's genes. `rho_cohort` is the median of `rho_own` over the samples. `rho_used = max(rho_cohort, rho_own, RHO_MIN)`. A sample with fewer than 20 genes at `MIN_DEPTH` is tested with `max(rho_cohort, RHO_MIN)`. Because `rho_cohort` is a median over the samples of the run, adding or removing samples changes `rho_cohort`, and so can change the p-values of the samples already analysed.
  - **Null sizes** (measured; 10 seeds; cohorts of 6 samples × 2000 genes):
    - true rho 0.005 / 0.02 / 0.05 without imbalance: 0.0141 / 0.0454 / 0.0462;
    - with 20 % imbalanced genes: 0.0147 / 0.0357 / 0.0238;
    - with 40 %: 0.0145 / 0.0188 / 0.0008;
    - few genes (6 samples × 20 / 30 / 60): 0.0233 / 0.0210 / 0.0265;
    - one sample with 15 genes: 0.0357;
    - heterogeneous cohort: 0.0223 overall and 0.0225 for the noisiest sample.
  - **Power:** for 0.75/0.25 genes at rho 0.02 the power is 0.8560, against an oracle with the true rho of 0.8751.
  - **Limit: many imbalanced genes inflate the dispersion, which is conservative but costs power.** Real imbalance and overdispersion both add a depth-independent spread to a / (a + b), and one sample's unphased data separate the two only by the shape of the distribution: imbalanced genes that fall inside the central region of the trimmed fit raise `rho_own`, and so `rho_cohort` and `rho_used`. Measured (same simulation, 0.75/0.25 genes):
    - true rho 0.05 with 40 % imbalanced genes: `rho_cohort` 0.1324, power 0.1460 against an oracle 0.5961;
    - true rho 0.05 with 20 %: power 0.4882 against 0.5957;
    - true rho 0.02 with 40 %: power 0.7843 against 0.8748;
    - 20 genes per sample (rho 0.02, 20 %): power 0.6330 against 0.8765, because noisy per-sample estimates enter the maximum rule.
    - These losses add no false positives: the null sizes stay below 0.05, down to 0.0008 in the worst case (rho 0.05, 40 %).
  - **Comparison with the unphased path** (fragment-level simulation in which read pairs are shared between SNPs; genes with 4 or more SNPs within 600 transcript bases):
    - haplotype-test power 0.7877 against 0.5873 for the per-SNP tests combined with ACAT (Rmd 02) at 300 fragments per gene, and 0.1495 against 0.0357 at 40;
    - at 40 fragments per gene, 0.9962 of 8-SNP genes are tested by the haplotype test and 0.4169 by the unphased path. That path needs one SNP with at least `MIN_DEPTH` reads.
    - The haplotype test's null size is 0.0184 (300) and 0.0042 (40).
    - The simulation puts all SNPs of a gene into one haplotype block (every fragment that covers any of its SNPs counts). With unphased genotypes, SNPs that no read pair links fall into separate blocks and `phaser_gene_ae` keeps only the most-covered one, so the haplotype-test power and tested fractions above are optimistic for that case.
  - **Runtime:** 55 s and 2027 MB for one sample of 20,000 genes on one core.
- **Excel row limit.** Rmd 05 writes its full gene tables as `.tsv.gz` because an xlsx sheet holds at most 1,048,576 rows. The `SNP` and `Gene` sheets of Rmd 02 and Rmd 04 share the same Excel limit of 1,048,576 rows per sheet, which a very large cohort can exceed; their checkpoints (`.rds`) always hold the full tables.

---

## Step 15 — Submission order and summary report

### Submission order (one block, every job `sbatch -p bcc`)

Write all scripts to `{RESULTS_DIR}/scripts/` (Steps 10-13, and Steps 16 and 17 for the selected analyses), run `mkdir -p {RESULTS_DIR}/logs {RESULTS_DIR}/tmp`, then submit the whole chain with `--parsable` job ids and `--dependency=afterok` so that no job starts before its inputs exist and none starts after a failure. Write the block below to `{RESULTS_DIR}/scripts/submit_chain.sh`, substitute the placeholders, keep or delete the optional mouse helper lines, show it to the user, and run it once from `{CWD}` with `bash {RESULTS_DIR}/scripts/submit_chain.sh` (it only calls `sbatch`, so it may run on the login node). It stops with an error, and submits nothing more, when a script it needs is missing or empty (`NEED`) or when `sbatch` returns no job id (`got`); it also stops when `{ANALYSES}` is empty or holds any word other than `per-sample`, `reciprocal`, `differential` (for example menu numbers), and when `{RESULTS_DIR}/logs/chain_job_ids.tsv` already exists (a chain has been submitted from this results directory). It records every job id in `{RESULTS_DIR}/logs/chain_job_ids.tsv` (script name, job id), which `wait_chain.sh` below reads.

**Running the wizard again into an existing results directory.** `{RESULTS_DIR}` is `results/{TODAY}_{WD_NAME}`, so a second run on the same day uses the same directory. Step 0 and Step 10 never delete an existing results directory: the scripts and Rmds the wizard writes again replace the old ones of the same name, and everything else (count tables, BAMs, logs, `chain_job_ids.tsv`) stays. `submit_chain.sh` then refuses to run because `chain_job_ids.tsv` exists; look at `squeue -u $USER` first, and only when the old jobs are finished or cancelled move the file away (`mv {RESULTS_DIR}/logs/chain_job_ids.tsv {RESULTS_DIR}/logs/chain_job_ids_previous.tsv`) and run the script again.

```bash
#!/bin/bash
# submit_chain.sh: run once from {CWD} with `bash {RESULTS_DIR}/scripts/submit_chain.sh`; it only calls sbatch.
set -u
S={RESULTS_DIR}/scripts
IDS={RESULTS_DIR}/logs/chain_job_ids.tsv
die() { echo "ERROR: $*" >&2; exit 1; }
got() { [ -n "$2" ] || die "sbatch returned no job id for $1; nothing after it was submitted (cancel the jobs already listed in $IDS)"; printf '%s\t%s\n' "$1" "$2" >> "$IDS"; }
MODE="{MODE}"; ANALYSES="{ANALYSES}"   # MODE: f1 or outbred; ANALYSES: the words selected in Step 7, e.g. "per-sample reciprocal differential"
case "$MODE" in
  f1)      PREP_SCRIPT=prep_f1_reference.sh; ARRAY_SCRIPT=align_count_f1.sh ;;
  outbred) PREP_SCRIPT=prep_genotypes.sh;    ARRAY_SCRIPT=align_wasp_count.sh ;;
  *)       die "MODE must be f1 or outbred, got '$MODE'" ;;
esac
case "$ANALYSES" in "") die "ANALYSES is empty: write the Step 7 words (per-sample, reciprocal, differential), separated by spaces" ;; esac
for w in $ANALYSES; do case "$w" in per-sample|reciprocal|differential) ;; *) die "ANALYSES contains '$w'; the accepted words are per-sample, reciprocal and differential, separated by spaces (not menu numbers or other spellings)" ;; esac; done
NEED="$PREP_SCRIPT $ARRAY_SCRIPT run_01_import_qc.sh run_02_imbalance.sh"
case " $ANALYSES " in *reciprocal*) [ "$MODE" = f1 ] || die "reciprocal analysis is F1 only"; NEED="$NEED run_03_reciprocal.sh" ;; esac
case " $ANALYSES " in *differential*) NEED="$NEED run_04_differential.sh" ;; esac
MISSING=$(for f in $NEED; do [ -s "$S/$f" ] || echo "$f"; done)
[ -z "$MISSING" ] || die "MISSING in $S: $(echo $MISSING) - nothing submitted"
[ ! -e "$IDS" ] || die "a chain has been submitted from this results directory ($IDS exists): check squeue, then remove or rename that file to submit again"
: > "$IDS"
DEP=""; R3="not selected"; R4="not selected"   # DEP stays empty unless the mouse helper below is used
# ONLY IF the mouse helper is used (F1, no parental VCF): array of chromosomes -> concat, then prep depends on the concat job.
# The helper reads neither the FASTA nor its .fai (fixed GRCm39 chromosome list), so it can run before the prep job,
# which downloads/decompresses the FASTA, writes the .fai and compares the parental VCF's contigs with it.
# Without the helper, delete the next four lines; the scripts do not exist and sbatch would fail.
NEED="$NEED extract_mgp_parental_vcf.sh concat_mgp_parental_vcf.sh"; MISSING=$(for f in $NEED; do [ -s "$S/$f" ] || echo "$f"; done); [ -z "$MISSING" ] || die "MISSING in $S: $(echo $MISSING) - nothing submitted"
X=$(sbatch -p bcc --parsable $S/extract_mgp_parental_vcf.sh); got extract_mgp_parental_vcf.sh "$X"
C=$(sbatch -p bcc --parsable --dependency=afterok:$X $S/concat_mgp_parental_vcf.sh); got concat_mgp_parental_vcf.sh "$C"
DEP="--dependency=afterok:$C"
# 1. prep job (PREP_SCRIPT was set from MODE above)
P=$(sbatch -p bcc --parsable $DEP $S/$PREP_SCRIPT); got $PREP_SCRIPT "$P"
# 2. per-sample array job (ARRAY_SCRIPT was set from MODE above)
A=$(sbatch -p bcc --parsable --dependency=afterok:$P $S/$ARRAY_SCRIPT); got $ARRAY_SCRIPT "$A"
# 3. Rmd 01 (import, filters, reference-bias QC) after every array task succeeded
R1=$(sbatch -p bcc --parsable --dependency=afterok:$A $S/run_01_import_qc.sh); got run_01_import_qc.sh "$R1"
# 4. Rmd 02 (per-sample imbalance, xlsx, summary_numbers.tsv) after Rmd 01
R2=$(sbatch -p bcc --parsable --dependency=afterok:$R1 $S/run_02_imbalance.sh); got run_02_imbalance.sh "$R2"
# 5. Rmd 03 only when reciprocal F1 was selected (Step 7): after Rmd 01 (it reads only the Rmd 01 checkpoint; runs beside Rmd 02)
case " $ANALYSES " in *reciprocal*)
  R3=$(sbatch -p bcc --parsable --dependency=afterok:$R1 $S/run_03_reciprocal.sh); got run_03_reciprocal.sh "$R3" ;;
esac
# 6. Rmd 04 only when differential ASE was selected (Step 7): after Rmd 01 (runs beside Rmd 02 and Rmd 03)
case " $ANALYSES " in *differential*)
  R4=$(sbatch -p bcc --parsable --dependency=afterok:$R1 $S/run_04_differential.sh); got run_04_differential.sh "$R4" ;;
esac
echo "submitted (script, job id):"; cat "$IDS"
echo "Rmd 03: $R3; Rmd 04: $R4"
```

Do not edit the block except for the mouse helper: substitute `{RESULTS_DIR}`, `{CWD}`, `{MODE}` (`f1` or `outbred`) and `{ANALYSES}` (the analyses selected in Step 7 written as the words `per-sample`, `reciprocal`, `differential`, separated by spaces, for example `per-sample reciprocal differential`; never the menu numbers). `PREP_SCRIPT` and `ARRAY_SCRIPT` are set from `MODE` (F1: `prep_f1_reference.sh`, `align_count_f1.sh`; outbred: `prep_genotypes.sh`, `align_wasp_count.sh`), and the Rmd 03 and Rmd 04 submissions and their names in `NEED` are conditional on `ANALYSES`; a job that is not selected leaves `R3` or `R4` as "not selected", which the last line prints (distinct from a failed submission, which stops with an error). When the mouse helper is used, keep its four lines; without it, delete the four lines after the "ONLY IF" comment (`NEED=`, `X=`, `C=`, `DEP=`) and `DEP` stays empty. Rmd 03 and Rmd 04 read only the Rmd 01 checkpoint (`ase_checkpoint.rds`) and check its constants themselves, so they depend on Rmd 01 and run beside Rmd 02; a failure of Rmd 02 does not stop them.

**Waiting.** Write `{RESULTS_DIR}/scripts/wait_chain.sh` and run it with `bash` (in the background; it only calls `squeue` and `scancel`). It waits until no job of `chain_job_ids.tsv` is left in the queue, checking every 60 s for at most 4 h. When every job still queued is pending with the reason `DependencyNeverSatisfied` (an upstream job failed, so they can never start), it cancels them and stops waiting. It then prints one line per job: `ok` for an Rmd whose result file is newer than `chain_job_ids.tsv`, `FAILED` with the first error line of its log, `CANCELLED` with the upstream job that failed (the job just before the first cancelled one in chain order), and `finished` for the prep, helper and array jobs (with the first error line of their log, if there is one). If `submit_chain.sh` itself stopped with an error, do not wait: cancel the jobs it lists in `chain_job_ids.tsv` (`cut -f2 {RESULTS_DIR}/logs/chain_job_ids.tsv | xargs -r scancel`), fix the cause, move the old `chain_job_ids.tsv` away (`mv {RESULTS_DIR}/logs/chain_job_ids.tsv {RESULTS_DIR}/logs/chain_job_ids_previous.tsv`, otherwise `submit_chain.sh` refuses) and run it again.

**On TIMEOUT (exit 2)** (the chain is still in the queue after 4 h, for example a long index build plus a large array): do not write the summary page. Report the job ids that `wait_chain.sh` printed as still queued or running (script name and id, from `chain_job_ids.tsv`); do not cancel them, they are working. Tell the user that the chain keeps running and that waiting is resumed by running `bash {RESULTS_DIR}/scripts/wait_chain.sh` again (it re-reads `chain_job_ids.tsv` and waits for another 4 h at most); then continue only when it returns 0 or 1.

```bash
#!/bin/bash
# wait_chain.sh: waits for every job of chain_job_ids.tsv (at most 4 h), cancels jobs that can never start, reports each job.
set -u
RD={RESULTS_DIR}
IDS=$RD/logs/chain_job_ids.tsv
JOBS=$(cut -f2 "$IDS" | paste -sd, -)
STUCK=""; n=0
PAT='error|halted|due to time limit|cancelled|out of memory|oom-kill'   # first matching line of a job log (a time-limit kill has no "error" word)
while :; do
  Q=$(squeue -h -j "$JOBS" -o '%i %r' 2>/dev/null)
  [ -z "$Q" ] && break
  if ! printf '%s\n' "$Q" | grep -qv ' DependencyNeverSatisfied$'; then   # only never-satisfiable jobs are left
    STUCK=$(printf '%s\n' "$Q" | awk '{sub(/_.*/, "", $1); print $1}' | sort -u | paste -sd' ' -)
    scancel $STUCK; break
  fi
  n=$((n + 1)); if [ $n -ge 240 ]; then echo "TIMEOUT after 4 h; still in the queue:"; printf '%s\n' "$Q"; exit 2; fi
  sleep 60
done
result_of() { case $1 in run_01_*) echo ase_checkpoint.rds ;; run_02_*) echo summary_numbers.tsv ;;
                         run_03_*) echo summary_numbers_reciprocal.tsv ;; run_04_*) echo summary_numbers_differential.tsv ;; esac; }
prev=""; up=""; bad=0
while IFS=$'\t' read -r name id; do
  if [[ " $STUCK " == *" $id "* ]]; then
    [ -n "$up" ] || up=${prev:-unknown}
    echo "$name $id CANCELLED: never satisfiable, upstream failure: $up"; bad=1
  elif [ -n "$(result_of "$name")" ]; then
    if [ "$RD/$(result_of "$name")" -nt "$IDS" ]; then echo "$name $id ok"
    else echo "$name $id FAILED: $(grep -h -m1 -iE "$PAT" $RD/logs/*_${id}*.out 2>/dev/null | head -1)"; bad=1; fi
  else e=$(grep -h -m1 -iE "$PAT" $RD/logs/*_${id}*.out 2>/dev/null | head -1); echo "$name $id finished${e:+; first error line of its log: $e}"; fi
  prev="$name $id"
done < "$IDS"
exit $bad
```

Write the summary page only after `wait_chain.sh` has returned. A failure of the prep job, the array, Rmd 01 or Rmd 02 stops the run: show the error line, fix the cause and re-submit (below). A failure of Rmd 03 or Rmd 04 does not invalidate the per-sample results: write the summary page anyway and say there which analysis failed, with its error line. If only Rmd 02 failed, Rmd 03 and 04 are still valid: fix it and re-submit Rmd 02 alone.

**Re-submitting part of the chain.** Copy `submit_chain.sh` to `submit_rmds.sh`, delete the lines of the jobs that need not run again (and their names in `NEED`), and remove `--dependency=afterok:$A` from the Rmd 01 line (or, when Rmd 01 is not re-run, `--dependency=afterok:$R1` from the lines of the Rmds that are); first move the old job list away (`mv {RESULTS_DIR}/logs/chain_job_ids.tsv {RESULTS_DIR}/logs/chain_job_ids_previous.tsv`; the `IDS` guard of the copied script refuses otherwise), run `submit_rmds.sh`, then `wait_chain.sh`: the new `submit_rmds.sh` writes a new `chain_job_ids.tsv`, so the wait covers the new jobs only. `afterok` on the array job id means every array task must succeed; if one fails, Rmd 01 and every Rmd after it stay pending with `DependencyNeverSatisfied` (`wait_chain.sh` cancels them; if you do it by hand, cancel them all: `awk -F'\t' '$1 ~ /^run_0/ {print $2}' {RESULTS_DIR}/logs/chain_job_ids.tsv | xargs -r scancel`). Fix the failed sample (its log is under `{RESULTS_DIR}/logs`), re-submit the failed task, wait until it has succeeded, then submit Rmd 01, Rmd 02 and the selected Rmd 03 / Rmd 04 again with `submit_rmds.sh`.

**Dropping a failed sample.** If a sample cannot be fixed (for example it has no counted sites, or no sites left after filtering in Rmd 01), the supported way out is to remove it from the analysis, never to edit the Rmds: with the user's confirmation, remove that sample's row from `{SAMPLES_CSV}` (keep a copy of the original sheet next to it), cancel every pending Rmd job (see "Re-submitting part of the chain"), and re-run from Rmd 01 with `submit_rmds.sh` (Rmd 01 without a dependency on the old array job, then Rmd 02, and the selected Rmd 03 / Rmd 04 after Rmd 01), whether an array task or only an Rmd failed (only Rmd 01 reads `{SAMPLES_CSV}`; Rmd 02, 03, 04 and 05 read its checkpoint, so they see the dropped sample until Rmd 01 has run again). The Rmds read the tables by the sample names in `{SAMPLES_CSV}`, so the dropped sample's files are simply ignored (if its table exists, Rmd 01 lists it as a table of a sample not in the sheet). Re-submitting the array job is not needed, because the other samples' tables already exist; if the array job is re-run anyway, `{ARRAY_N}` must be updated to the new number of rows (task numbers follow the rows of the edited sheet). Rmd 03 and Rmd 04 are re-run after Rmd 01 in the same way (a dropped sample can change whether a design is still valid, so re-check Step 7's preconditions first). Rmd 02 is always re-run too, even when it had succeeded: the thinning pools the depth over all samples of the project, so removing a sample can change which SNP a dense gene keeps, and with it the gene p-values of the remaining samples. When phASER was selected, Rmd 05 (Step 19) is re-run after Rmd 02 for a similar reason: its dispersion uses the median over all samples (`rho_cohort`), so removing a sample can change the p-values of the remaining samples. Say in the summary page which sample was dropped and why.

### Summary report `{WD_NAME}_summary_report.html`

After `wait_chain.sh` has returned (every submitted job has left the queue: Rmd 02 and, when selected, Rmd 03 and Rmd 04), write the standalone page `{RESULTS_DIR}/{WD_NAME}_summary_report.html`. It needs **no R** and no external dependency: inline CSS only, no scripts, no fonts, no images from a URL. Read `{RESULTS_DIR}/summary_numbers.tsv` with `awk -F'\t'` (for example `awk -F'\t' 'NR>1 {print $1, $3, $4}'`, so no number is transcribed by hand; columns: `sample`, `condition`, `filtered_sites`, `sig_snps`, `genes_tested`, `sig_genes`, `rho_used`, `rho_robust`, `rho_h0_naive`, `mean_ref_frac`, `ref_frac_before_wasp`, `ref_frac_after_wasp`, `bias_flag`; `rho_used` is the overdispersion the SNP tests used (the F1 gene tests use `rho_gene`, see Step 13), `rho_robust` the estimate before the `RHO_MIN` floor (F1: bias-corrected free-mean fit; outbred: trimmed central-sites fit), `rho_h0_naive` the all-sites H0 fit shown for comparison only, and the two `ref_frac_*_wasp` columns are `NA` in F1 mode), and write the HTML yourself; never open R on the login node to produce it. The page contains:

- a header with the project title, `{MODE}`, the strain names (F1: REF = `{STRAIN_A}`, ALT = `{STRAIN_B}`) or REF/ALT, the constants of Step 8 and the date;
- the **reference-bias flags, shown prominently** at the top: one box that lists every sample with `bias_flag` = `TRUE` (with its `mean_ref_frac`) in a warning colour and the sentence "read the allelic ratios of these samples with caution", or "no sample is flagged" when none is; the per-sample table repeats the flag in its own column, with the same colour;
- a per-sample table from `summary_numbers.tsv`: sample, condition, filtered sites, significant SNPs, genes tested, significant genes, `rho_used`, `rho_robust`, `rho_h0_naive` (labelled "naive, not used"), mean REF fraction, bias flag;
- outbred only: the WASP note, "Alignments whose vW tag is 2-7 were removed before counting; Rmd 01 shows how many were removed and from which allele. WASP does not promise less bias." Outbred gene calls are also labelled "unphased, no direction". The per-sample table also shows `ref_frac_before_wasp` and `ref_frac_after_wasp` under the heading "REF fraction before / after WASP (reported for information, not a gate)";
- **relative links** (bare file names, no URL prefix and no absolute path, so the page works from the local filesystem and over a web mount as long as it sits in `{RESULTS_DIR}` beside the files) as cards to: `{TODAY}_{WD_NAME}_01_import_qc.html`, `{TODAY}_{WD_NAME}_02_imbalance.html`, `{TODAY}_{WD_NAME}_ASE_imbalance.xlsx`, and the two figure PDFs `{TODAY}_{WD_NAME}_ASE_sites_vs_depth.pdf` and `{TODAY}_{WD_NAME}_ASE_genes.pdf`.

If Rmd 03 or Rmd 04 was selected but failed, its section below is replaced by a box that names the analysis, quotes the error line of its log (`{RESULTS_DIR}/logs/run_03_reciprocal_<jobid>.out` or `run_04_differential_<jobid>.out`) and says that the per-sample results above are not affected; the page has no link to its missing files.

### Reciprocal F1 section (only when Rmd 03 ran)

Read `{RESULTS_DIR}/summary_numbers_reciprocal.tsv` with `awk -F'\t'` (columns `test, genes_tested, genes_not_tested, sig_genes, direction_pos, n_pos, direction_neg, n_neg, phi_common, samples_AxB, samples_BxA`) and show one table row per test ("strain effect" and "parent-of-origin effect") with the tested and significant genes and the counts per direction. Add the sentences: "p is the `{STRAIN_A}` fraction; d = +1 when `{STRAIN_A}` is the mother (AxB). The strain effect is sensitive to reference mapping bias (see the flags above); the parent-of-origin effect is not. X, Y and MT are not tested." Add link cards to `{TODAY}_{WD_NAME}_03_reciprocal.html`, `{TODAY}_{WD_NAME}_ASE_reciprocal.xlsx` and `{TODAY}_{WD_NAME}_ASE_reciprocal.pdf`.

### Differential ASE section (only when Rmd 04 ran)

Read `{RESULTS_DIR}/summary_numbers_differential.tsv` with `awk -F'\t'` (columns `contrast, mode, level, tested, not_tested, sig, n_up, n_down, dispersion`) and show one row per contrast and level. Add, in F1 mode: "n_up / n_down: `{STRAIN_A}` fraction higher / lower in the tested condition than in `{REF_CONDITION}`". In outbred mode: "unphased: SNP directions are counted only when all tested individuals agree; gene calls have no direction" (the gene rows have `NA` in `n_up` and `n_down`). Mark a contrast with `tested` = 0 as "nothing tested" (the reasons are in the Rmd 04 HTML). A condition other than `{REF_CONDITION}` that has no row at all was not tested at the design check (fewer than 2 samples in a condition, or, outbred, fewer than 2 individuals sampled in both conditions); list it with the reason Rmd 04 printed for it. Take the reason from the rendered output line, not from the echoed code (which also contains the words): `grep -o '## Contrasts not tested:[^<]*' {RESULTS_DIR}/{TODAY}_{WD_NAME}_04_differential.html` (knitr prefixes printed output with `## `). Add: "X, Y and MT are not tested in either mode." Add link cards to `{TODAY}_{WD_NAME}_04_differential.html`, `{TODAY}_{WD_NAME}_ASE_differential.xlsx` and `{TODAY}_{WD_NAME}_ASE_differential.pdf`.

**Verify before finishing.** Extract every `href` of the page and check with a shell loop that each target exists next to the page (`grep -o 'href="[^"]*"' {RESULTS_DIR}/{WD_NAME}_summary_report.html | sed 's/href="//;s/"$//' | while read -r f; do [ -s "{RESULTS_DIR}/$f" ] || echo "MISSING $f"; done`); also check that no `href` or `src` starts with `http` or `/`. Fix the page (or report the missing file) until nothing is printed. Then tell the user the paths of the summary page, the two stage HTML files, the xlsx and `summary_numbers.tsv`, and, when they ran, the Rmd 03/04 HTML files, xlsx files and summary TSVs.

---

## Step 16 — Rmd 03: reciprocal F1 (strain and parent-of-origin effects)

Only when "Reciprocal F1 analysis" was selected in Step 7 (`{MODE}` = `f1`). Write `{CWD}/{TODAY}_{WD_NAME}_03_reciprocal.Rmd` with the same conventions as Step 12, and the same `{AUTHOR}` and `{PROJECT_TITLE}`. It loads `ase_checkpoint.rds` (Rmd 01), checks the constants against it, pastes the statistics functions of Step 14 into the chunk marked below, and writes these files to `{RESULTS_DIR}`: `{TODAY}_{WD_NAME}_ASE_reciprocal.xlsx` (sheets `Gene`, `SNPs_used`, `Design`, `Excluded`, `Summary`), `{TODAY}_{WD_NAME}_ASE_reciprocal.pdf`, `ase_reciprocal_checkpoint.rds` and `summary_numbers_reciprocal.tsv`.

- **Model.** Per gene, `logit(p) = b0 + b1 * d`, where p is the **strain-A (`{STRAIN_A}`, REF) fraction** and d = +1 for `AxB` (strain A is the mother) and −1 for `BxA`. When several conditions exist, sum-to-zero condition terms are added: b0 is then the unweighted average over the conditions, and b1 is one parent-of-origin effect common to all conditions (the model has no direction x condition term). At least one condition must contain samples of both directions; if every condition holds a single direction, condition and direction cannot be told apart and the Rmd stops. **b0 is the strain effect** (cis-regulatory divergence): b0 > 0 means "`{STRAIN_A}` higher". **b1 is the parent-of-origin effect** (imprinting): b1 > 0 means "maternal higher", and `maternal_frac = plogis(b1)`.
- **Rows.** One row per sample and gene: the strain-A and total counts summed over the gene's SNPs after `thin_snps` (one SNP per `THIN_BP` window in exon coordinates, deepest SNP first, the same SNPs in every sample). ASEReadCounter counts a read pair at every SNP it covers, so the gene test uses no SNP pair that one read pair can span. A gene with a single SNP is tested normally, because the replication comes from the samples. SNPs that lie in two genes are left out.
- **Dispersion and tests.** `ase_glm_test` (Step 14): the pooled dispersion between replicate animals (`phi_common`), each gene's own estimate when it is larger and the gene has at least 2 residual df, floored at `RHO_MIN`. b0 and b1 are tested by likelihood-ratio tests (1 df each), with BH within each test. `sig_strain` requires `padj_strain < FDR_SIG` and `|frac_A - 0.5| >= ABS_DEV_SIG`; `sig_parent_of_origin` requires `padj_parent_of_origin < FDR_SIG` and `|maternal_frac - 0.5| >= ABS_DEV_SIG`. Monoallelic genes (complete imprinting) have their coefficients at the ±15 bound (`status` `ok_at_bound`) with a valid p-value.
- **Requirements.** Every sample's `cross_direction` must be exactly `AxB` or `BxA`, and each direction needs at least 2 samples; otherwise the Rmd stops. Per gene, at least 2 samples with coverage in each direction; genes below that are listed in `not_tested`.
- **Excluded chromosomes.** X, Y and MT (and `chrX`, `chrY`, `chrM`) are not tested. A hemizygous X in males and the maternally inherited MT would look like maternal imprinting. They are counted in the `Excluded` sheet. Only the contigs present in the GTF gene map are tested (a SNP must lie in an exon of a GTF gene); alt, random and unplaced contigs are treated as autosomes (`chrom_class` recognises only the X, Y and MT names), so a gene placed on such a contig is tested like an autosomal gene.
- **Reference bias.** b0 is exactly what a mapping bias toward strain A would imitate, so the Rmd 01 bias flags are printed first and repeated next to the strain results. b1 is not affected by a constant mapping bias toward one strain (it cancels between the two directions).

````rmd
---
title: "{PROJECT_TITLE} - ASE reciprocal F1"
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
DATE_TAG    <- "{TODAY}_{WD_NAME}"
MIN_DEPTH   <- {MIN_DEPTH}
FDR_SIG     <- {FDR_SIG}
ABS_DEV_SIG <- {ABS_DEV_SIG}
RHO_MIN     <- {RHO_MIN}             # floor for the dispersion used by the tests
BIAS_TOL    <- {BIAS_TOL}
THIN_BP     <- {THIN_BP}             # SNPs closer than this (exon coordinates) can share a read pair: one kept per window
if (MODE != "f1") stop("Rmd 03 (reciprocal F1) needs MODE = f1", call. = FALSE)
ck <- readRDS(file.path(RESULTS_DIR, "ase_checkpoint.rds"))
same <- c(MODE = identical(ck$constants$MODE, MODE),   # swapped strain names would invert every strain label
          STRAIN_A = identical(ck$constants$STRAIN_A, STRAIN_A), STRAIN_B = identical(ck$constants$STRAIN_B, STRAIN_B),
          vapply(c("MIN_DEPTH", "FDR_SIG", "ABS_DEV_SIG", "BIAS_TOL"),
                 function(k) isTRUE(all.equal(ck$constants[[k]], get(k))), logical(1)))
if (!all(same)) stop("constants differ from the Rmd 01 checkpoint (", paste(names(same)[!same], collapse = ", "),
                     "); re-render Rmd 01 with the same values", call. = FALSE)
sites <- ck$sites; bias <- ck$bias
knitr::kable(bias, digits = 4, caption = paste0("Reference-bias diagnostic from Rmd 01 (REF = ", STRAIN_A, "); a screen, not a test. ",
                                                 "The strain effect b0 is sensitive to it; the parent-of-origin effect b1 is not"))
if (any(bias$flagged)) cat("FLAGGED samples:", paste(bias$sample[bias$flagged], collapse = ", "), "- read the strain effects with caution\n")
```

## Statistics functions

```{r stats}
# <<< paste here the code of Step 14 between the two marker lines (the marker lines themselves are not pasted) >>>
```

## Design

```{r design}
dirs <- unique(ck$samples[, c("sample", "condition", "cross_direction")])
bad <- dirs$sample[!dirs$cross_direction %in% c("AxB", "BxA")]
if (length(bad) > 0) stop("cross_direction must be AxB (", STRAIN_A, " mother) or BxA (", STRAIN_B, " mother); other values for: ",
                          paste(bad, collapse = ", "), call. = FALSE)
n_dir <- table(factor(dirs$cross_direction, levels = c("AxB", "BxA")))
if (any(n_dir < 2)) stop("the reciprocal analysis needs at least 2 samples in each cross direction (AxB: ", n_dir[["AxB"]],
                         ", BxA: ", n_dir[["BxA"]], ")", call. = FALSE)
dirs$d <- ifelse(dirs$cross_direction == "AxB", 1, -1)   # d = +1: strain A is the mother (maternal allele = STRAIN_A)
cl_all <- sort(unique(dirs$condition))   # the model's design over all samples: intercept, d, sum-to-zero condition terms
Xd <- cbind(1, dirs$d)
if (length(cl_all) > 1) Xd <- cbind(Xd, contr.sum(length(cl_all))[match(dirs$condition, cl_all), , drop = FALSE])
if (qr(Xd)$rank < ncol(Xd))
  stop("condition is confounded with cross direction (every condition holds a single direction), so the parent-of-origin ",
       "effect cannot be told apart from the condition effect. Samples per condition and direction:\n",
       paste(capture.output(print(table(condition = dirs$condition, direction = dirs$cross_direction))), collapse = "\n"),
       "\nAt least one condition needs samples of both directions", call. = FALSE)
knitr::kable(dirs, caption = "Samples, condition and cross direction (d = +1: AxB, strain A maternal)")
```

## Gene map in exon coordinates and SNP thinning

Only SNPs in exons of GTF genes are used: only the contigs present in the GTF gene map are tested. X, Y and MT (`chrX`, `chrY`, `chrM`) are set aside and counted below; alt, random and unplaced contigs are treated as autosomes.

```{r genemap}
gtf <- rtracklayer::import(GTF_PATH)
if (is.null(gtf$gene_id)) stop("the GTF has no gene_id attribute", call. = FALSE)
ex <- gtf[gtf$type == "exon" & !is.na(gtf$gene_id)]
exd <- as.data.frame(GenomicRanges::reduce(split(ex, ex$gene_id)))              # exon union per gene: group_name = gene_id
exd <- exd[order(exd$group_name, exd$start), ]
exd$offset <- ave(exd$width, exd$group_name, FUN = function(w) cumsum(w) - w)    # exonic bases of the gene before this exon
pos <- unique(sites[, c("contig", "position")])
gr <- GenomicRanges::GRanges(pos$contig, IRanges::IRanges(pos$position, width = 1))
exr <- GenomicRanges::GRanges(as.character(exd$seqnames), IRanges::IRanges(exd$start, exd$end))
hits <- GenomicRanges::findOverlaps(gr, exr); qh <- S4Vectors::queryHits(hits); sh <- S4Vectors::subjectHits(hits)
if (length(qh) == 0) stop("no SNP lies in an exon of the GTF (contig names: counts ", paste(head(unique(pos$contig), 3), collapse = ", "),
                          "; GTF ", paste(head(unique(as.character(exd$seqnames)), 3), collapse = ", "), ")", call. = FALSE)
snp_gene <- data.frame(contig = pos$contig[qh], position = pos$position[qh], gene_id = exd$group_name[sh],
                       exon_pos = exd$offset[sh] + pos$position[qh] - exd$start[sh] + 1, stringsAsFactors = FALSE)
ng <- table(paste(snp_gene$contig, snp_gene$position))
snp_gene$multi_gene <- as.vector(ng[paste(snp_gene$contig, snp_gene$position)]) > 1
use <- sites %>% dplyr::inner_join(dplyr::filter(snp_gene, !multi_gene), by = c("contig", "position")) %>%
  dplyr::mutate(chrom = chrom_class(contig))
excluded_chrom <- dplyr::count(dplyr::filter(use, chrom != "autosome"), chrom, name = "sample_site_rows")
use <- dplyr::filter(use, chrom == "autosome")
kept <- use %>% dplyr::group_by(gene_id, contig, position, exon_pos) %>% dplyr::summarise(depth = sum(total), .groups = "drop") %>%
  dplyr::group_by(gene_id) %>% dplyr::mutate(keep = thin_snps(exon_pos, depth, THIN_BP)) %>% dplyr::ungroup() %>% dplyr::filter(keep)
gene_rows <- use %>% dplyr::semi_join(kept, by = c("gene_id", "contig", "position")) %>%
  dplyr::group_by(gene_id, sample, condition, cross_direction) %>%
  dplyr::summarise(y = sum(ref_n), n = sum(total), .groups = "drop") %>% dplyr::filter(n > 0)
if (nrow(kept) == 0)
  stop("no autosomal SNP in a single GTF gene is left after the gene map: ", nrow(unique(snp_gene[, c("contig", "position")])),
       " SNP positions lie in GTF exons, ", sum(snp_gene$multi_gene), " SNP-gene pairs are in overlapping genes, ",
       sum(excluded_chrom$sample_site_rows), " sample-site rows are on X, Y or MT. Check the contigs and genes of the GTF (",
       GTF_PATH, ") against the count tables", call. = FALSE)
cat(sum(snp_gene$multi_gene), "SNP-gene pairs in overlapping genes left out;", nrow(kept), "SNPs kept after thinning in",
    length(unique(kept$gene_id)), "genes\n")
if (nrow(excluded_chrom) > 0) knitr::kable(excluded_chrom, caption = "Rows on X, Y or MT: not tested (hemizygous X and maternal MT mimic imprinting)")
```

## Per-gene model and tests

```{r test}
build <- lapply(split(gene_rows, gene_rows$gene_id), function(r) {
  nA <- sum(r$cross_direction == "AxB"); nB <- sum(r$cross_direction == "BxA")
  if (nA < 2 || nB < 2) return(list(unit = NULL, why = sprintf("fewer than 2 samples with coverage in a direction (AxB %d, BxA %d)", nA, nB)))
  X <- cbind("(Intercept)" = 1, d = ifelse(r$cross_direction == "AxB", 1, -1))
  cl <- sort(unique(r$condition))
  if (length(cl) > 1) {   # condition adjustment, sum-to-zero coding: b0 = unweighted average over the conditions, b1 common to all
    C <- contr.sum(length(cl)); colnames(C) <- paste0("cond_", cl[-length(cl)])
    X <- cbind(X, C[match(r$condition, cl), , drop = FALSE])
  }
  if (qr(X)$rank < ncol(X)) return(list(unit = NULL, why = "design not estimable (condition confounded with direction in the covered samples)"))
  list(unit = list(y = r$y, n = r$n, X = X), why = NA_character_)
})
units <- lapply(Filter(function(b) !is.null(b$unit), build), `[[`, "unit")
not_tested <- data.frame(gene_id = names(build), reason = vapply(build, `[[`, character(1), "why"), stringsAsFactors = FALSE)
not_tested <- not_tested[!is.na(not_tested$reason), , drop = FALSE]
if (length(units) == 0) stop("no gene can be tested; genes per reason: ",
                             paste(sprintf("%s: %d", names(table(not_tested$reason)), table(not_tested$reason)), collapse = "; "), call. = FALSE)
res <- ase_glm_test(units, list(strain = "(Intercept)", parent_of_origin = "d"), RHO_MIN)
st <- table(res$status)
if (!any(res$status %in% c("ok", "ok_at_bound")))
  stop("no gene could be tested (fit status: ", paste(sprintf("%s %d", names(st), st), collapse = ", "), "; phi_common ", res$phi_common[1], ")",
       call. = FALSE)
obs <- gene_rows %>% dplyr::group_by(gene_id) %>%
  dplyr::summarise(n_samples_AxB = sum(cross_direction == "AxB"), n_samples_BxA = sum(cross_direction == "BxA"),
                   frac_A_AxB = sum(y[cross_direction == "AxB"]) / sum(n[cross_direction == "AxB"]),
                   frac_A_BxA = sum(y[cross_direction == "BxA"]) / sum(n[cross_direction == "BxA"]), .groups = "drop")
gene <- res %>% dplyr::rename(gene_id = unit) %>%
  dplyr::left_join(dplyr::count(kept, gene_id, name = "n_snps_used"), by = "gene_id") %>%
  dplyr::left_join(obs, by = "gene_id") %>%
  dplyr::mutate(b0_strain = `beta_(Intercept)`, b1_parent_of_origin = beta_d,
                frac_A = plogis(b0_strain), maternal_frac = plogis(b1_parent_of_origin),
                padj_strain = p.adjust(p_strain, "BH"), padj_parent_of_origin = p.adjust(p_parent_of_origin, "BH"),
                sig_strain = !is.na(padj_strain) & padj_strain < FDR_SIG & abs(frac_A - 0.5) >= ABS_DEV_SIG,
                sig_parent_of_origin = !is.na(padj_parent_of_origin) & padj_parent_of_origin < FDR_SIG & abs(maternal_frac - 0.5) >= ABS_DEV_SIG,
                direction_strain = dplyr::case_when(!sig_strain ~ "none", b0_strain > 0 ~ paste(STRAIN_A, "higher"), TRUE ~ paste(STRAIN_B, "higher")),
                direction_parent_of_origin = dplyr::case_when(!sig_parent_of_origin ~ "none", b1_parent_of_origin > 0 ~ "maternal higher",
                                                              TRUE ~ "paternal higher"),
                category = dplyr::case_when(sig_strain & sig_parent_of_origin ~ "strain and parent-of-origin", sig_strain ~ "strain (cis)",
                                            sig_parent_of_origin ~ "parent-of-origin", TRUE ~ "none")) %>%
  dplyr::arrange(gene_id)
cat("Dispersion between replicate animals (phi_common):", signif(gene$phi_common[1], 4), "; floor RHO_MIN =", RHO_MIN, "\n")
knitr::kable(dplyr::count(gene, status), caption = "Fit status (ok_at_bound = monoallelic gene, coefficient at the bound, p valid)")
knitr::kable(head(dplyr::filter(gene, category != "none"), 40), digits = 4,
             caption = paste0("Significant genes (first 40): strain effect = ", STRAIN_A, " fraction averaged over directions; ",
                              "parent-of-origin effect = maternal fraction", if (any(bias$flagged)) "; some samples flagged for reference bias" else ""))
```

## Summary

```{r summary}
row_for <- function(test, sig, b, pos_lab, neg_lab) {
  p <- gene[[paste0("p_", test)]]
  data.frame(test = test, genes_tested = sum(!is.na(p)), genes_not_tested = nrow(not_tested) + sum(is.na(p)), sig_genes = sum(sig),
             direction_pos = pos_lab, n_pos = sum(sig & b > 0, na.rm = TRUE), direction_neg = neg_lab, n_neg = sum(sig & b < 0, na.rm = TRUE),
             phi_common = gene$phi_common[1], samples_AxB = n_dir[["AxB"]], samples_BxA = n_dir[["BxA"]], stringsAsFactors = FALSE)
}
summary_rec <- rbind(row_for("strain", gene$sig_strain, gene$b0_strain, paste(STRAIN_A, "higher"), paste(STRAIN_B, "higher")),
                     row_for("parent_of_origin", gene$sig_parent_of_origin, gene$b1_parent_of_origin, "maternal higher", "paternal higher"))
knitr::kable(summary_rec, digits = 4, caption = "Summary (FDR_SIG, ABS_DEV_SIG as in the constants block)")
```

## Figure

```{r figure}
gene_plot <- dplyr::filter(gene, !is.na(p_strain))
stopifnot(sum(gene_plot$sig_strain) == summary_rec$sig_genes[1], sum(gene_plot$sig_parent_of_origin) == summary_rec$sig_genes[2])
cat("Figure sig counts equal Summary counts: TRUE\n")
p1 <- ggplot(gene_plot, aes(frac_A, maternal_frac, colour = category)) + geom_hline(yintercept = 0.5, linetype = 2) +
  geom_vline(xintercept = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.8) +
  labs(x = paste0(STRAIN_A, " fraction (strain effect, plogis(b0))"), y = "maternal fraction (parent-of-origin effect, plogis(b1))",
       colour = NULL, caption = paste0("sig: BH < ", FDR_SIG, " and |fraction - 0.5| >= ", ABS_DEV_SIG,
                                       if (any(bias$flagged)) "; some samples flagged for reference bias" else "")) + theme_bw()
print(p1)
ggsave(file.path(RESULTS_DIR, paste0(DATE_TAG, "_ASE_reciprocal.pdf")), p1, width = 9, height = 7)
```

## Export

```{r export}
out_names <- function(d) {
  names(d) <- sub("^frac_A$", paste0("frac_", STRAIN_A), names(d))
  names(d) <- sub("^frac_A_(AxB|BxA)$", paste0("frac_", STRAIN_A, "_\\1"), names(d))
  d
}
xlsx_file <- file.path(RESULTS_DIR, paste0(DATE_TAG, "_ASE_reciprocal.xlsx"))
wb <- openxlsx::createWorkbook()
for (nm in c("Gene", "SNPs_used", "Design", "Excluded", "Summary")) openxlsx::addWorksheet(wb, nm)
openxlsx::writeData(wb, "Gene", out_names(gene))
openxlsx::writeData(wb, "SNPs_used", kept)
openxlsx::writeData(wb, "Design", dirs)
openxlsx::writeData(wb, "Excluded", rbind(data.frame(what = sprintf("chromosome %s", excluded_chrom$chrom), n = excluded_chrom$sample_site_rows),
                                           data.frame(what = sprintf("not tested: %s", not_tested$reason), n = rep(1, nrow(not_tested)))) %>%
                                       dplyr::group_by(what) %>% dplyr::summarise(n = sum(n), .groups = "drop"))
openxlsx::writeData(wb, "Summary", summary_rec)
openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)
saveRDS(list(gene = gene, not_tested = not_tested, snps_used = kept, design = dirs, excluded_chrom = excluded_chrom, summary = summary_rec,
             constants = c(ck$constants, list(RHO_MIN = RHO_MIN, THIN_BP = THIN_BP))),
        file.path(RESULTS_DIR, "ase_reciprocal_checkpoint.rds"))
write.table(summary_rec, file.path(RESULTS_DIR, "summary_numbers_reciprocal.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", xlsx_file, "and summary_numbers_reciprocal.tsv\n")
sessionInfo()
```
````

Render it with `{RESULTS_DIR}/scripts/run_03_reciprocal.sh`. This is the same script as `run_01_import_qc.sh` (Step 12), with job name `ase_03_reciprocal`, log `run_03_reciprocal_%j.out`, the Rmd 03 file name, the same `--bind` rule and `-n 1 --mem=16G -t 4:00:00`. The Rmd runs single-threaded. The time covers genome-wide gene sets (the runtime check of the Step 14 unit tests, `ase-pipeline/tests/r/test_ase_stats_stage2.R`, measured 32 s for 5000 units x 12 rows with both tests on one core and projects 0.11 h for 60,000 genes; the 4 h request is deliberately conservative, because the projection was not run at genome scale). If the job stops with "cross_direction must be AxB", correct `{SAMPLES_CSV}` (Step 5), re-render Rmd 01, then Rmd 03. Submission order: Step 15.

---

## Step 17 — Rmd 04: differential ASE between conditions

Only when "Differential ASE between conditions" was selected in Step 7. Write `{CWD}/{TODAY}_{WD_NAME}_04_differential.Rmd` with the conventions of Step 12, and the same `{AUTHOR}` and `{PROJECT_TITLE}`. It loads `ase_checkpoint.rds` (Rmd 01), checks the constants (including the strain names) against it, pastes the statistics functions of Step 14 into the chunk marked below, and writes these files to `{RESULTS_DIR}`: `{TODAY}_{WD_NAME}_ASE_differential.xlsx` (sheets `SNP`, `Gene`, `Design`, `Excluded`, `Summary`), `{TODAY}_{WD_NAME}_ASE_differential.pdf`, `ase_differential_checkpoint.rds` and `summary_numbers_differential.tsv`. Every condition other than `{REF_CONDITION}` (Step 7) is compared with it, one contrast at a time (1-df likelihood-ratio test), with BH within each contrast and level (SNP or gene). The contrast is labelled `<condition> vs {REF_CONDITION}`, for example `treat vs ctrl`.

- **F1.** Per gene (rows = samples, strain-A and total counts summed over the thinned SNPs exactly as in Rmd 03) and per SNP (rows = samples), `logit(p) = b0 + b_condition * cond`. Here p is the strain-A (`{STRAIN_A}`, REF) fraction and cond = 1 for the tested condition. When both cross directions are among the contrast's samples, the cross-direction term `d` (+1 `AxB`, −1 `BxA`) is added, so that parent-of-origin genes are not mistaken for condition effects. If condition and direction are confounded (every sample of one condition has one direction and every sample of the other the other direction), the Rmd stops, because the two effects cannot be separated. The dispersion and the tests are those of Rmd 03 (`ase_glm_test`). `delta_frac = plogis(b0 + b_condition) - plogis(b0)`: **delta_frac > 0 means the `{STRAIN_A}` fraction is higher in the tested condition**. With the `d` term, b0 is the logit at d = 0 (halfway between the two directions), so frac_A_ref and frac_A_test are direction-averaged fractions of the reference and the tested condition; with unbalanced directions they differ from the observed pooled fractions. `sig` requires `padj < FDR_SIG` and `|delta_frac| >= ABS_DEV_SIG`. Per gene or SNP, at least 2 samples with coverage in each condition are needed; a unit whose covered samples make the design not estimable (condition confounded with direction in those samples) is listed in `not_tested` with that reason. If no gene has a valid fit (status `ok` or `ok_at_bound`) in any contrast, the Rmd stops and prints the fit status.
- **Several conditions.** A contrast in which nothing can be tested (for example no gene with 2 covered samples per condition) does not stop the Rmd when another contrast can be tested: it prints "nothing tested in contrast <contrast> (<level> level): <reasons>" and still gets its `Summary` rows (one per contrast and level) with `tested` = 0.
- **Outbred.** Only individuals sampled in both conditions of a contrast are used (the others are listed as unpaired in the `Design` sheet and do not change the results of the paired ones). `ase_paired_test` works per individual. First, a pair dispersion is estimated from all its SNPs (each SNP keeps its own REF fraction, floored at `RHO_MIN`; below 20 paired SNPs the individual is not tested). Then each SNP gets a likelihood-ratio test of one REF fraction against one per condition. The per-SNP statistic is the sum over the informative individuals (df = their number), and at least 2 individuals are needed. The test has no direction, because the allele that carries a regulatory variant differs between individuals: in one person the REF allele of a SNP can rise while in another it falls. `mean_delta` (REF fraction change, tested minus reference), `n_up` and `n_down` describe the individual changes. A significant SNP is labelled "REF higher (or lower) in the tested condition (all tested individuals)" only when all tested individuals agree (the individuals with a finite statistic for that SNP; those counted in `n_failed` do not enter the label), and "mixed (phase differs)" otherwise. Genes use `acat` over their SNP p-values and are labelled "unphased, no direction". `sig` requires `padj < FDR_SIG` and `max_abs_delta >= ABS_DEV_SIG`. If no SNP can be tested in any contrast, the Rmd stops; a single contrast with nothing tested (all its SNPs `too_few_individuals`, for example when one of its 2 paired individuals has fewer than 20 paired SNPs) prints "nothing tested in contrast ..." with the SNP status and the individuals without a pair dispersion, and its dispersion is reported as `phi_pair NA`.
- **Excluded chromosomes.** X, Y and MT are excluded in both modes, as in Rmd 03. Only the contigs present in the GTF gene map are tested at the gene level (a SNP must lie in an exon of a GTF gene); alt, random and unplaced contigs are treated as autosomes (`chrom_class` recognises only the X, Y and MT names). The `Excluded` sheets of the two Rmds count different rows: Rmd 03 counts the sample-site rows of SNPs in exons of a single GTF gene on X, Y or MT (the only SNPs it could test), and Rmd 04 counts every sample-site row on X, Y or MT (it also tests SNPs outside genes). The Rmd 01 bias flags are printed first.

````rmd
---
title: "{PROJECT_TITLE} - ASE differential between conditions"
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
MODE          <- "{MODE}"
STRAIN_A      <- "{STRAIN_A}"
STRAIN_B      <- "{STRAIN_B}"
RESULTS_DIR   <- "{RESULTS_DIR}"
GTF_PATH      <- "{GTF_PATH}"
DATE_TAG      <- "{TODAY}_{WD_NAME}"
MIN_DEPTH     <- {MIN_DEPTH}
FDR_SIG       <- {FDR_SIG}
ABS_DEV_SIG   <- {ABS_DEV_SIG}
RHO_MIN       <- {RHO_MIN}             # floor for the dispersion used by the tests
BIAS_TOL      <- {BIAS_TOL}
THIN_BP       <- {THIN_BP}             # F1 genes: SNPs closer than this (exon coordinates) can share a read pair: one kept per window
REF_CONDITION <- "{REF_CONDITION}"     # every other condition is compared with this one
stopifnot(MODE %in% c("f1", "outbred"))
ck <- readRDS(file.path(RESULTS_DIR, "ase_checkpoint.rds"))
same <- c(MODE = identical(ck$constants$MODE, MODE),   # swapped strain names would invert every strain label
          STRAIN_A = identical(ck$constants$STRAIN_A, STRAIN_A), STRAIN_B = identical(ck$constants$STRAIN_B, STRAIN_B),
          vapply(c("MIN_DEPTH", "FDR_SIG", "ABS_DEV_SIG", "BIAS_TOL"),
                 function(k) isTRUE(all.equal(ck$constants[[k]], get(k))), logical(1)))
if (!all(same)) stop("constants differ from the Rmd 01 checkpoint (", paste(names(same)[!same], collapse = ", "),
                     "); re-render Rmd 01 with the same values", call. = FALSE)
sites <- ck$sites; bias <- ck$bias
knitr::kable(bias, digits = 4, caption = paste0("Reference-bias diagnostic from Rmd 01 (REF = ", STRAIN_A, "); a screen, not a test"))
if (any(bias$flagged)) cat("FLAGGED samples:", paste(bias$sample[bias$flagged], collapse = ", "), "- read the fractions with caution\n")
```

## Statistics functions

```{r stats}
# <<< paste here the code of Step 14 between the two marker lines (the marker lines themselves are not pasted) >>>
```

## Design and contrasts

```{r design}
smp <- unique(ck$samples[, c("sample", "condition", "cross_direction", "individual")])
if (!REF_CONDITION %in% smp$condition) stop("reference condition '", REF_CONDITION, "' is not a condition in samples.csv (conditions: ",
                                            paste(sort(unique(smp$condition)), collapse = ", "), ")", call. = FALSE)
if (MODE == "f1") {
  bad <- smp$sample[!smp$cross_direction %in% c("AxB", "BxA")]
  if (length(bad) > 0) stop("cross_direction must be AxB (", STRAIN_A, " mother) or BxA (", STRAIN_B, " mother); other values for: ",
                            paste(bad, collapse = ", "), call. = FALSE)
} else {
  bad <- smp$sample[is.na(smp$individual) | smp$individual %in% c("", "NA")]
  if (length(bad) > 0) stop("outbred mode needs the individual of every sample; missing for: ", paste(bad, collapse = ", "), call. = FALSE)
}
levels_test <- setdiff(sort(unique(smp$condition)), REF_CONDITION)
if (length(levels_test) == 0) stop("differential ASE needs at least two conditions", call. = FALSE)
design <- list(); skipped <- character()
for (L in levels_test) {
  s <- smp[smp$condition %in% c(REF_CONDITION, L), ]
  if (MODE == "f1") {
    nc <- table(factor(s$condition, levels = c(REF_CONDITION, L)))
    if (any(nc < 2)) { skipped[L] <- sprintf("fewer than 2 samples in a condition (%s %d, %s %d)", REF_CONDITION, nc[[1]], L, nc[[2]]); next }
    use_d <- all(c("AxB", "BxA") %in% s$cross_direction)
    if (use_d && qr(cbind(1, s$condition == L, s$cross_direction == "AxB"))$rank < 3)
      stop("condition '", L, "' versus '", REF_CONDITION, "' is confounded with the cross direction (each condition has only one ",
           "direction), so a condition effect cannot be separated from a parent-of-origin effect. Samples per condition and direction:\n",
           paste(capture.output(print(table(condition = s$condition, direction = s$cross_direction))), collapse = "\n"),
           "\nAdd samples of both directions to a condition, or restrict samples.csv to one direction and re-run Rmd 01", call. = FALSE)
    design[[L]] <- list(samples = s, use_d = use_d)
  } else {
    pairs <- intersect(s$individual[s$condition == REF_CONDITION], s$individual[s$condition == L])
    if (length(pairs) < 2) { skipped[L] <- sprintf("fewer than 2 individuals sampled in both %s and %s", REF_CONDITION, L); next }
    design[[L]] <- list(samples = s[s$individual %in% pairs, ], unpaired = sort(setdiff(s$individual, pairs)))
    if (length(design[[L]]$unpaired) > 0) cat(paste(L, "vs", REF_CONDITION), "- unpaired individuals, not used:",
                                              paste(design[[L]]$unpaired, collapse = ", "), "\n")
  }
}
if (length(skipped) > 0) cat("Contrasts not tested:", paste(names(skipped), skipped, sep = ": ", collapse = "; "), "\n")
if (length(design) == 0) stop("no contrast can be tested: ", paste(names(skipped), skipped, sep = ": ", collapse = "; "), call. = FALSE)
design_tbl <- dplyr::bind_rows(lapply(names(design), function(L) dplyr::mutate(design[[L]]$samples, contrast = paste(L, "vs", REF_CONDITION),
  cross_direction_term = if (MODE == "f1") design[[L]]$use_d else NA,
  unpaired_individuals = if (MODE == "outbred") paste(design[[L]]$unpaired, collapse = ";") else NA_character_)))
knitr::kable(design_tbl, caption = "Samples per contrast (outbred: paired individuals only; unpaired ones listed)")
```

## Gene map in exon coordinates and SNP thinning

Only SNPs in exons of GTF genes enter the gene tables: only the contigs present in the GTF gene map are tested at the gene level. X, Y and MT (`chrX`, `chrY`, `chrM`) are set aside and counted below; alt, random and unplaced contigs are treated as autosomes. The per-SNP tests use every autosomal SNP.

```{r genemap}
gtf <- rtracklayer::import(GTF_PATH)
if (is.null(gtf$gene_id)) stop("the GTF has no gene_id attribute", call. = FALSE)
ex <- gtf[gtf$type == "exon" & !is.na(gtf$gene_id)]
exd <- as.data.frame(GenomicRanges::reduce(split(ex, ex$gene_id)))              # exon union per gene: group_name = gene_id
exd <- exd[order(exd$group_name, exd$start), ]
exd$offset <- ave(exd$width, exd$group_name, FUN = function(w) cumsum(w) - w)    # exonic bases of the gene before this exon
pos <- unique(sites[, c("contig", "position")])
gr <- GenomicRanges::GRanges(pos$contig, IRanges::IRanges(pos$position, width = 1))
exr <- GenomicRanges::GRanges(as.character(exd$seqnames), IRanges::IRanges(exd$start, exd$end))
hits <- GenomicRanges::findOverlaps(gr, exr); qh <- S4Vectors::queryHits(hits); sh <- S4Vectors::subjectHits(hits)
if (length(qh) == 0) stop("no SNP lies in an exon of the GTF (contig names: counts ", paste(head(unique(pos$contig), 3), collapse = ", "),
                          "; GTF ", paste(head(unique(as.character(exd$seqnames)), 3), collapse = ", "), ")", call. = FALSE)
snp_gene <- data.frame(contig = pos$contig[qh], position = pos$position[qh], gene_id = exd$group_name[sh],
                       exon_pos = exd$offset[sh] + pos$position[qh] - exd$start[sh] + 1, stringsAsFactors = FALSE)
ng <- table(paste(snp_gene$contig, snp_gene$position))
snp_gene$multi_gene <- as.vector(ng[paste(snp_gene$contig, snp_gene$position)]) > 1
sites <- dplyr::mutate(sites, chrom = chrom_class(contig))
excluded_chrom <- dplyr::count(dplyr::filter(sites, chrom != "autosome"), chrom, name = "sample_site_rows")
auto <- dplyr::filter(sites, chrom == "autosome")                                  # every autosomal SNP: the per-SNP tests
use <- dplyr::inner_join(auto, dplyr::filter(snp_gene, !multi_gene), by = c("contig", "position"))
kept <- use %>% dplyr::group_by(gene_id, contig, position, exon_pos) %>% dplyr::summarise(depth = sum(total), .groups = "drop") %>%
  dplyr::group_by(gene_id) %>% dplyr::mutate(keep = thin_snps(exon_pos, depth, THIN_BP)) %>% dplyr::ungroup() %>% dplyr::filter(keep)
gene_rows <- use %>% dplyr::semi_join(kept, by = c("gene_id", "contig", "position")) %>%
  dplyr::group_by(gene_id, sample, condition, cross_direction) %>%
  dplyr::summarise(y = sum(ref_n), n = sum(total), .groups = "drop") %>% dplyr::filter(n > 0)
if (nrow(kept) == 0)
  stop("no autosomal SNP in a single GTF gene is left after the gene map: ", nrow(unique(snp_gene[, c("contig", "position")])),
       " SNP positions lie in GTF exons, ", sum(snp_gene$multi_gene), " SNP-gene pairs are in overlapping genes, ",
       sum(excluded_chrom$sample_site_rows), " sample-site rows are on X, Y or MT. Check the contigs and genes of the GTF (",
       GTF_PATH, ") against the count tables", call. = FALSE)
snp_to_gene <- dplyr::distinct(use, SNP, gene_id)   # outbred genes: acat over all their SNPs (valid under dependence, no thinning)
cat(sum(snp_gene$multi_gene), "SNP-gene pairs in overlapping genes left out;", nrow(kept), "SNPs kept after thinning in",
    length(unique(kept$gene_id)), "genes (F1 gene rows)\n")
if (nrow(excluded_chrom) > 0) knitr::kable(excluded_chrom, caption = "Rows on X, Y or MT: not tested")
```

## F1: per-gene and per-SNP condition effect

delta_frac > 0: the strain-A fraction is higher in the tested condition. With both cross directions in a contrast, `d` (+1 AxB, −1 BxA) is in the model, so parent-of-origin effects are not read as condition effects.

```{r f1, eval = (MODE == "f1")}
f1_units <- function(rows, id, L, use_d) {
  b <- lapply(split(rows, rows[[id]]), function(r) {
    cc <- as.numeric(r$condition == L)
    if (sum(cc == 1) < 2 || sum(cc == 0) < 2) return(list(unit = NULL, why = "fewer than 2 samples with coverage in a condition"))
    X <- cbind("(Intercept)" = 1, cond = cc)
    if (use_d && length(unique(r$cross_direction)) == 2) X <- cbind(X, d = ifelse(r$cross_direction == "AxB", 1, -1))
    if (qr(X)$rank < ncol(X)) return(list(unit = NULL, why = "design not estimable (condition confounded with direction in the covered samples)"))
    list(unit = list(y = r$y, n = r$n, X = X), why = NA_character_)
  })
  ok <- vapply(b, function(x) !is.null(x$unit), logical(1))
  list(units = lapply(b[ok], `[[`, "unit"),
       not_tested = data.frame(id = names(b)[!ok], reason = vapply(b[!ok], `[[`, character(1), "why"), stringsAsFactors = FALSE))
}
f1_table <- function(res, L, idname) {
  names(res)[names(res) == "unit"] <- idname
  res %>% dplyr::mutate(contrast = paste(L, "vs", REF_CONDITION), b0 = `beta_(Intercept)`, b_condition = beta_cond,
                        frac_A_ref = plogis(b0), frac_A_test = plogis(b0 + b_condition), delta_frac = frac_A_test - frac_A_ref,
                        padj = p.adjust(p_condition, "BH"), sig = !is.na(padj) & padj < FDR_SIG & abs(delta_frac) >= ABS_DEV_SIG,
                        direction = dplyr::case_when(!sig ~ "none", b_condition > 0 ~ paste0(STRAIN_A, " fraction higher in ", L),
                                                     TRUE ~ paste0(STRAIN_A, " fraction lower in ", L))) %>%
    dplyr::relocate(contrast, dplyr::all_of(idname), b0, b_condition, frac_A_ref, frac_A_test, delta_frac, p_condition, padj, sig,
                    direction, status, phi_common, phi_unit, phi_used, n_rows)
}
snp_rows <- dplyr::transmute(auto, SNP, sample, condition, cross_direction, y = ref_n, n = total) %>% dplyr::filter(n > 0)
gene_l <- list(); snp_l <- list()
nt_l <- list(data.frame(contrast = character(0), level = character(0), id = character(0), reason = character(0), stringsAsFactors = FALSE))
for (L in names(design)) {
  ks <- design[[L]]$samples$sample; ct <- paste(L, "vs", REF_CONDITION)
  gu <- f1_units(dplyr::filter(gene_rows, sample %in% ks), "gene_id", L, design[[L]]$use_d)
  su <- f1_units(dplyr::filter(snp_rows, sample %in% ks), "SNP", L, design[[L]]$use_d)
  if (length(gu$units) > 0) gene_l[[L]] <- f1_table(ase_glm_test(gu$units, list(condition = "cond"), RHO_MIN), L, "gene_id")
  if (length(su$units) > 0) snp_l[[L]] <- f1_table(ase_glm_test(su$units, list(condition = "cond"), RHO_MIN), L, "SNP")
  nt_l <- c(nt_l, list(dplyr::mutate(gu$not_tested, contrast = ct, level = "gene"), dplyr::mutate(su$not_tested, contrast = ct, level = "SNP")))
  for (lv in c("gene", "SNP")) {   # several conditions: one contrast can have nothing tested while another is fine
    rs <- if (lv == "gene") gene_l[[L]] else snp_l[[L]]; nt <- if (lv == "gene") gu$not_tested else su$not_tested
    if (is.null(rs) || !any(rs$status %in% c("ok", "ok_at_bound")))
      cat("nothing tested in contrast", ct, paste0("(", lv, " level): "),
          if (!is.null(rs)) paste("fit status", paste(sprintf("%s %d", names(table(rs$status)), table(rs$status)), collapse = ", "))
          else if (nrow(nt) > 0) paste(sprintf("%s: %d", names(table(nt$reason)), table(nt$reason)), collapse = "; ")
          else "no covered unit", "\n")
  }
}
gene <- dplyr::bind_rows(gene_l); snp <- dplyr::bind_rows(snp_l)
not_tested <- dplyr::bind_rows(nt_l) %>% dplyr::select(contrast, level, id, reason)
if (nrow(gene) == 0) stop("no gene can be tested in any contrast; genes per reason: ",
                          paste(sprintf("%s: %d", names(table(not_tested$reason[not_tested$level == "gene"])),
                                        table(not_tested$reason[not_tested$level == "gene"])), collapse = "; "), call. = FALSE)
st <- table(gene$status)
if (!any(gene$status %in% c("ok", "ok_at_bound")))
  stop("no gene could be tested in any contrast (fit status: ", paste(sprintf("%s %d", names(st), st), collapse = ", "),
       "; phi_common ", paste(signif(unique(gene$phi_common), 4), collapse = ", "), ")", call. = FALSE)
dispersion <- dplyr::distinct(gene, contrast, phi_common)
knitr::kable(dispersion, digits = 4, caption = paste("Dispersion between replicate animals per contrast (gene level); floor RHO_MIN =", RHO_MIN))
knitr::kable(dplyr::count(gene, contrast, status), caption = "Fit status per contrast (gene level)")
knitr::kable(head(dplyr::filter(gene, sig), 40), digits = 4, caption = paste0("Significant genes (first 40); delta_frac = change of the ", STRAIN_A,
                                                                             " fraction (tested - reference)", if (any(bias$flagged)) "; some samples flagged for reference bias" else ""))
```

## Outbred: paired per-SNP test and gene combination

Unphased: a SNP's direction is reported only when all tested individuals agree; genes carry no direction.

```{r outbred, eval = (MODE == "outbred")}
snp_l <- list(); gene_l <- list(); phi_l <- list()
nt_l <- list(data.frame(contrast = character(0), level = character(0), id = character(0), reason = character(0), stringsAsFactors = FALSE))
for (L in names(design)) {
  s <- design[[L]]$samples; ct <- paste(L, "vs", REF_CONDITION)
  dp <- auto %>% dplyr::filter(sample %in% s$sample) %>%
    dplyr::transmute(snp = SNP, individual, cond = as.numeric(condition == L), y = ref_n, n = total)
  r <- ase_paired_test(as.data.frame(dp), RHO_MIN)
  st <- r$snp %>% dplyr::rename(SNP = snp) %>%
    dplyr::mutate(contrast = ct, padj = p.adjust(p, "BH"), sig = !is.na(padj) & padj < FDR_SIG & max_abs_delta >= ABS_DEV_SIG,
                  direction = dplyr::case_when(!sig ~ "none", n_down == 0 ~ paste0("REF higher in ", L, " (all tested individuals)"),
                                               n_up == 0 ~ paste0("REF lower in ", L, " (all tested individuals)"), TRUE ~ "mixed (phase differs)")) %>%
    dplyr::select(contrast, SNP, n_individuals, n_failed, stat, df, p, padj, mean_delta, max_abs_delta, n_up, n_down, sig, direction, status)
  gt <- st %>% dplyr::inner_join(snp_to_gene, by = "SNP") %>% dplyr::filter(!is.na(p)) %>% dplyr::group_by(gene_id) %>%
    dplyr::summarise(n_snps = dplyr::n(), acat_p = acat(p), max_abs_delta = max(max_abs_delta), .groups = "drop") %>%
    dplyr::mutate(contrast = ct, padj = p.adjust(acat_p, "BH"), sig = !is.na(padj) & padj < FDR_SIG & max_abs_delta >= ABS_DEV_SIG,
                  note = "unphased, no direction") %>%
    dplyr::select(contrast, gene_id, n_snps, acat_p, padj, max_abs_delta, sig, note)
  snp_l[[L]] <- st; gene_l[[L]] <- gt; phi_l[[L]] <- dplyr::mutate(r$phi, contrast = ct)
  if (!any(st$status == "ok"))   # several conditions: one contrast can have nothing tested while another is fine
    cat("nothing tested in contrast", ct, "(SNP and gene level): SNP status",
        paste(sprintf("%s %d", names(table(st$status)), table(st$status)), collapse = ", "),
        "; individuals without a pair dispersion (fewer than 20 paired SNPs):", paste(r$phi$individual[is.na(r$phi$phi_pair)], collapse = ", "), "\n")
  bad <- st$status != "ok"
  ng <- setdiff(unique(snp_to_gene$gene_id[snp_to_gene$SNP %in% st$SNP]), gt$gene_id)   # genes without any tested SNP
  nt_l <- c(nt_l, list(data.frame(contrast = rep(ct, sum(bad)), level = rep("SNP", sum(bad)), id = st$SNP[bad], reason = st$status[bad]),
                       data.frame(contrast = rep(ct, length(ng)), level = rep("gene", length(ng)), id = ng, reason = rep("no tested SNP", length(ng)))))
}
snp <- dplyr::bind_rows(snp_l); gene <- dplyr::bind_rows(gene_l); not_tested <- dplyr::bind_rows(nt_l)
dispersion <- dplyr::bind_rows(phi_l)
if (!any(snp$status == "ok"))
  stop("no SNP could be tested in any contrast (SNP status: ", paste(sprintf("%s %d", names(table(snp$status)), table(snp$status)), collapse = ", "),
       "; individuals without a pair dispersion (fewer than 20 paired SNPs): ",
       paste(dispersion$individual[is.na(dispersion$phi_pair)], collapse = ", "), ")", call. = FALSE)
knitr::kable(dispersion, digits = 4, caption = "Pair dispersion per individual (floored at RHO_MIN; NA = fewer than 20 paired SNPs, not tested)")
knitr::kable(head(dplyr::filter(snp, sig), 40), digits = 4, caption = paste0("Significant SNPs (first 40); unphased: direction only when all tested individuals agree",
                                                                            if (any(bias$flagged)) "; some samples flagged for reference bias" else ""))
```

## Summary

```{r summary}
contrasts <- paste(names(design), "vs", REF_CONDITION)   # one Summary row per contrast and level, also when nothing was tested
cnt <- function(tab, lev) {
  rows <- list()
  for (ct in contrasts) {
    t <- if (nrow(tab) > 0) tab[tab$contrast == ct, , drop = FALSE] else tab   # 0 rows: nothing tested in this contrast
    pcol <- if ("p_condition" %in% names(t)) t$p_condition else if ("acat_p" %in% names(t)) t$acat_p else t$p
    ids <- if (lev == "gene") t$gene_id else t$SNP
    nt_ids <- not_tested$id[not_tested$contrast == ct & not_tested$level == lev]
    up <- if (MODE == "f1") sum(t$sig & t$delta_frac > 0) else if (lev == "SNP") sum(t$sig & t$n_down == 0) else NA_integer_
    dn <- if (MODE == "f1") sum(t$sig & t$delta_frac < 0) else if (lev == "SNP") sum(t$sig & t$n_up == 0) else NA_integer_
    disp <- if (MODE == "f1") {
      phc <- if (nrow(t) > 0) t$phi_common[1] else NA_real_
      if (is.na(phc)) "phi_common NA" else sprintf("phi_common %.4f", phc)
    } else {
      ph <- dispersion$phi_pair[dispersion$contrast == ct]
      if (all(is.na(ph))) "phi_pair NA" else sprintf("phi_pair %.4f-%.4f", min(ph, na.rm = TRUE), max(ph, na.rm = TRUE))
    }
    rows[[ct]] <- data.frame(contrast = ct, mode = MODE, level = lev, tested = sum(!is.na(pcol)),
                             not_tested = length(unique(c(nt_ids, ids[is.na(pcol)]))),   # a unit listed in not_tested and with p NA counts once
                             sig = sum(t$sig), n_up = up, n_down = dn, dispersion = disp, stringsAsFactors = FALSE)
  }
  dplyr::bind_rows(rows)
}
summary_diff <- dplyr::bind_rows(cnt(gene, "gene"), cnt(snp, "SNP"))
knitr::kable(summary_diff, caption = paste0("Summary (n_up / n_down: F1 = ", STRAIN_A, " fraction higher / lower in the tested condition; ",
                                            "outbred SNPs = REF higher / lower in all tested individuals; outbred genes: no direction)"))
```

## Figure

```{r figure}
fig <- if (MODE == "f1") dplyr::filter(gene, !is.na(p_condition)) %>% dplyr::transmute(contrast, x = delta_frac, p = p_condition, sig) else
  dplyr::filter(snp, !is.na(p)) %>% dplyr::transmute(contrast, x = mean_delta, p = p, sig)
lev <- if (MODE == "f1") "gene" else "SNP"
sum_n <- setNames(summary_diff$sig[summary_diff$level == lev], summary_diff$contrast[summary_diff$level == lev])
fig_n <- vapply(names(sum_n), function(ct) sum(fig$sig[fig$contrast == ct]), integer(1))   # 0 for a contrast with nothing tested
stopifnot(identical(unname(fig_n), unname(as.integer(sum_n))))
cat("Figure sig counts equal Summary counts: TRUE\n")
p1 <- ggplot(fig, aes(x, -log10(p), colour = sig)) + geom_vline(xintercept = 0, linetype = 2) + geom_point(alpha = 0.7) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick")) + facet_wrap(~contrast) +
  labs(x = if (MODE == "f1") paste0("change of the ", STRAIN_A, " fraction (tested - reference)") else "mean change of the REF fraction (unphased)",
       y = "-log10 p", caption = paste0("sig: BH < ", FDR_SIG, " and |change| >= ", ABS_DEV_SIG,
                                        if (any(bias$flagged)) "; some samples flagged for reference bias" else "")) + theme_bw()
print(p1)
ggsave(file.path(RESULTS_DIR, paste0(DATE_TAG, "_ASE_differential.pdf")), p1, width = 10, height = 7)
```

## Export

```{r export}
xlsx_file <- file.path(RESULTS_DIR, paste0(DATE_TAG, "_ASE_differential.xlsx"))
wb <- openxlsx::createWorkbook()
for (nm in c("SNP", "Gene", "Design", "Excluded", "Summary")) openxlsx::addWorksheet(wb, nm)
openxlsx::writeData(wb, "SNP", snp); openxlsx::writeData(wb, "Gene", gene); openxlsx::writeData(wb, "Design", design_tbl)
openxlsx::writeData(wb, "Excluded", dplyr::bind_rows(
  dplyr::transmute(excluded_chrom, what = sprintf("chromosome %s (sample-site rows)", chrom), n = sample_site_rows),
  dplyr::count(not_tested, contrast, level, reason, name = "n") %>% dplyr::transmute(what = sprintf("%s, %s not tested: %s", contrast, level, reason), n)))
openxlsx::writeData(wb, "Summary", summary_diff)
openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)
saveRDS(list(snp = snp, gene = gene, not_tested = not_tested, design = design_tbl, excluded_chrom = excluded_chrom, snp_gene = snp_to_gene,
             dispersion = dispersion, summary = summary_diff,
             constants = c(ck$constants, list(RHO_MIN = RHO_MIN, THIN_BP = THIN_BP, REF_CONDITION = REF_CONDITION))),
        file.path(RESULTS_DIR, "ase_differential_checkpoint.rds"))
write.table(summary_diff, file.path(RESULTS_DIR, "summary_numbers_differential.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", xlsx_file, "and summary_numbers_differential.tsv\n")
sessionInfo()
```
````

Render it with `{RESULTS_DIR}/scripts/run_04_differential.sh`. This is the same script as `run_01_import_qc.sh` (Step 12), with job name `ase_04_differential`, log `run_04_differential_%j.out`, the Rmd 04 file name, the same `--bind` rule and `-n 1 --mem=16G -t 4:00:00` (single-threaded; F1 tests every gene and every SNP). Runtime: the same unit-test benchmark as for Rmd 03 gives 6.4 ms per unit (32 s for 5000 units, both tests, one core), so the F1 per-SNP level projects to about 0.5 h for 300,000 covered SNPs; this is a projection from the unit test, not a genome-scale run. If the job stops with "is confounded with the cross direction", show the message to the user. The design cannot answer the question; never edit the Rmd. If it stops with "cross_direction must be AxB", correct `{SAMPLES_CSV}` (Step 5), re-render Rmd 01, then Rmd 04. If it stops with "reference condition '...' is not a condition in samples.csv (conditions: ...)", nothing was written: the message lists the conditions of `samples.csv`. Ask the user which of them is the reference, set `{REF_CONDITION}` (Step 7) to it, write Rmd 04 again from this step with the new value (a new generated file, not a hand edit of the old one) and submit `run_04_differential.sh` again; Rmd 01 need not be re-rendered, because its checkpoint does not depend on `{REF_CONDITION}`. Only if the condition names in `samples.csv` themselves are wrong, correct `{SAMPLES_CSV}` and re-render Rmd 01 before Rmd 04. Submission order: Step 15.

---

## Step 18 — phASER: installation and per-sample haplotype counts (outbred, optional)

Only when phASER was selected in Step 7 (`{MODE}` = `outbred`, paired-end data).

phASER phases each sample's heterozygous SNPs from the reads themselves: a read pair that covers two SNPs shows which alleles lie on one molecule. It then counts the read pairs of each haplotype, once per read pair, over the SNPs of a haplotype block, and `phaser_gene_ae` turns the blocks into one haplotype-A and one haplotype-B count per gene. These counts are tested in Rmd 05 (Step 19), next to the unphased per-SNP and ACAT results of Rmd 02, never instead of them.

**Which phASER.** Use the repository https://github.com/secastel/phaser (maintained by PEJ Lab as a "Fast Beta"; Python 3; no releases or tags), pinned to commit `aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301` (phASER v1.2.0, 2026-03-21). Use this repository, never the bioconda package `phaser`, which is version 0.1.1 for Python 2.7.

**Environment.** Installed once per user by `setup_phaser_env.sh` into `{PHASER_HOME}` (Step 7) and reused by every project:
- `src/phaser`: the clone at the pinned commit;
- `env`: a conda environment, about 630 MB;
- `environment.phaser.yml`: the pinned package list;
- `conda_list.txt`;
- `.ase_pipeline_phaser`: an empty marker that the directory belongs to this installation;
- `install_ok.txt`: written last, after the checks pass.

`{PHASER_HOME}` must be an absolute path to a new or empty directory that holds nothing but phASER. The script deletes and rebuilds `env`, `src/phaser` and `pkgs` inside it, so it refuses a directory that is not empty and has no `.ase_pipeline_phaser` marker. The installation needs about 2.6 GB free while it runs: the conda package cache `pkgs` (about 1.9 GB) is removed at the end, and about 650 MB stay (the environment and the 22 MB clone).

The installation is a compute job: conda and `git` need the internet, which the compute nodes have, and the job takes 10 to 20 minutes, mostly conda's classic solver. `module add miniconda3/v4` in this job is the only module besides `singularity/3.10.4` that the skill ever loads. The phASER scripts of a project load no module and do not activate the environment: they put `{PHASER_HOME}/env/bin` first on `PATH`, and `samtools`, `bcftools`, `tabix` and `python` come from it. Every script that runs the environment's Python exports `PYTHONNOUSERSITE=1` and unsets `PYTHONPATH` and `PYTHONHOME`. Otherwise the user site-packages directory of the same Python version (`~/.local/lib/python3.14`, filled by pip's per-user installs) would shadow the environment's packages (it fails with `GLIBC_2.27 not found` on this cluster), and an inherited Python search path would do the same.

### `setup_phaser_env.sh`

Written to `{RESULTS_DIR}/scripts/setup_phaser_env.sh` (substitute the placeholders). `install_ok.txt` starts with two key lines: line 1 is exactly `commit <PHASER_COMMIT>` and line 2 is exactly `env_sha256 <sha256 of the package list>`, the hash of `environment.phaser.yml` as the script writes it (`sha256sum environment.phaser.yml`). `submit_chain.sh` (Step 15) submits the job only when these two lines do not both match the script; otherwise the chain skips it. The script makes the same test first and, when both lines match, prints "phASER already installed" and exits 0. Two projects that install into the same `{PHASER_HOME}` at the same time are stopped by the lock directory `{PHASER_HOME}/.installing`. If a job was killed, remove that directory by hand.

```bash
#!/bin/bash
#SBATCH -J ase_phaser_setup
#SBATCH -N 1 -p bcc
#SBATCH -n 4 --mem=16G -t 2:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={USER_EMAIL}
#SBATCH -o {RESULTS_DIR}/logs/setup_phaser_env_%j.out
set -uo pipefail
# One-time installation of phASER (pinned commit) and its pinned conda environment under {PHASER_HOME}; reused by every project.
# never the bioconda package "phaser" (version 0.1.1, Python 2.7). Runs on a compute node (internet for git and conda).
PHASER_HOME="{PHASER_HOME}"
PHASER_REPO=https://github.com/secastel/phaser.git
PHASER_COMMIT=aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301   # no releases or tags upstream: pinned by SHA (Step 18, bumping)
MARK=.ase_pipeline_phaser   # ownership marker: env, src/phaser and pkgs are deleted and rebuilt only in a directory that carries it
die() { echo "ERROR: $*" >&2; exit 1; }
case "$PHASER_HOME" in /*) ;; *) die "PHASER_HOME must be an absolute path, not '$PHASER_HOME'" ;; esac
ENV_YML=$(cat <<'YML'
name: phaser-pinned
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.14.7
  - numpy=2.5.3
  - scipy=1.18.1
  - pandas=3.0.6
  - pysam=0.24.1
  - intervaltree=3.2.1
  - samtools=1.24
  - bcftools=1.24
  - htslib=1.24
  - bedtools=2.31.1
  - pip
YML
)
ENV_SHA=$(printf '%s\n' "$ENV_YML" | sha256sum | cut -d ' ' -f 1)   # line 2 of install_ok.txt: a changed package list reinstalls
if [ "$(sed -n 1p "$PHASER_HOME/install_ok.txt" 2>/dev/null)" = "commit $PHASER_COMMIT" ] && [ "$(sed -n 2p "$PHASER_HOME/install_ok.txt" 2>/dev/null)" = "env_sha256 $ENV_SHA" ]; then
  echo "phASER already installed in $PHASER_HOME"; exit 0
fi
mkdir -p "$PHASER_HOME" || die "cannot create $PHASER_HOME"
L=$(ls -A "$PHASER_HOME") || die "cannot list $PHASER_HOME"
[ -e "$PHASER_HOME/$MARK" ] || [ -z "$L" ] || die "$PHASER_HOME is not empty and has no $MARK (it was not created by setup_phaser_env.sh): choose another, new or empty directory for phASER"
touch "$PHASER_HOME/$MARK" || die "cannot write $PHASER_HOME/$MARK"
mkdir "$PHASER_HOME/.installing" 2>/dev/null || die "another installation is running in $PHASER_HOME (if not, remove $PHASER_HOME/.installing)"
trap 'rmdir "$PHASER_HOME/.installing" 2>/dev/null' EXIT
rm -f "$PHASER_HOME/install_ok.txt"
rm -rf "$PHASER_HOME/src/phaser" "$PHASER_HOME/env" "$PHASER_HOME/pkgs"
git clone -q "$PHASER_REPO" "$PHASER_HOME/src/phaser" || die "git clone of $PHASER_REPO failed"
git -C "$PHASER_HOME/src/phaser" checkout -q "$PHASER_COMMIT" || die "commit $PHASER_COMMIT not found"
[ "$(git -C "$PHASER_HOME/src/phaser" rev-parse HEAD)" = "$PHASER_COMMIT" ] || die "the checked-out commit differs from $PHASER_COMMIT"
printf '%s\n' "$ENV_YML" > "$PHASER_HOME/environment.phaser.yml" || die "cannot write $PHASER_HOME/environment.phaser.yml"
module add miniconda3/v4 || die "cannot load the miniconda3/v4 module"
source /home/software/conda/miniconda3/bin/condainit || die "cannot initialise conda"
export CONDA_PKGS_DIRS="$PHASER_HOME/pkgs"
conda env create -p "$PHASER_HOME/env" -f "$PHASER_HOME/environment.phaser.yml" || die "conda env create failed"
rm -rf "$PHASER_HOME/pkgs"     # the package cache (about 1.9 GB) is not needed after the installation
export PYTHONNOUSERSITE=1      # the user site-packages directory (~/.local/lib/python3.14) would shadow the environment
unset PYTHONPATH PYTHONHOME    # an inherited Python search path would shadow the environment the same way
export PATH="$PHASER_HOME/env/bin:$PATH"; PY="$PHASER_HOME/env/bin/python"
"$PY" -c "import numpy, scipy, pysam, pandas, intervaltree" || die "the environment's Python packages do not import"
"$PY" "$PHASER_HOME/src/phaser/phaser/phaser.py" --help > /dev/null || die "phaser.py --help failed"
"$PY" "$PHASER_HOME/src/phaser/phaser_gene_ae/phaser_gene_ae.py" --help > /dev/null || die "phaser_gene_ae.py --help failed"
for t in python samtools bcftools bgzip tabix bedtools; do [ -x "$PHASER_HOME/env/bin/$t" ] || die "$t is missing from $PHASER_HOME/env/bin"; done
conda list -p "$PHASER_HOME/env" > "$PHASER_HOME/conda_list.txt" || die "conda list failed"
{ echo "commit $PHASER_COMMIT"; echo "env_sha256 $ENV_SHA"; "$PY" -c 'import sys; print("python", sys.version.split()[0])'; date -u; } > "$PHASER_HOME/install_ok.txt" \
  || die "cannot write $PHASER_HOME/install_ok.txt"
echo "phASER installed in $PHASER_HOME"
```

The versions in `environment.phaser.yml` are the ones conda resolved from phASER's own unpinned `environment.yml` when this installation was tested; they are pinned so that a later installation gets the tested set. A pin to another Python (for example 3.12) has not been tested.

**Bumping the pinned commit.** phASER has no releases or tags, so the skill pins a commit SHA, and a newer commit is adopted on purpose, never by accident. To adopt one:
1. Change `PHASER_COMMIT` in `setup_phaser_env.sh` and in `phaser_count.sh`, and the commit string in the phASER check of `submit_chain.sh` (Step 15). All three must agree.
2. Install into a new `{PHASER_HOME}` and re-record `ase-pipeline/tests/fixtures/phaser_help.txt` and `phaser_gene_ae_help.txt` from it (the checker compares every phASER flag in this file with them).
3. Re-run the synthetic acceptance (README, Validation status), and only then use it.

Change the version list of `environment.phaser.yml` the same way, after a test installation: a changed package list changes this hash (line 2 of `install_ok.txt`), so the new hash must also replace the old one wherever the chain compares it. A changed commit or package list makes `setup_phaser_env.sh` reinstall into `{PHASER_HOME}`: the old clone and environment are removed, so a project that still needs the old commit or versions must use another `{PHASER_HOME}`.

### `phaser_count.sh`

A SLURM array job, one task per data row of `{SAMPLES_CSV}`, submitted by `submit_chain.sh` (Step 15) after the per-sample array job (Step 11) and, when it runs, the installation job. Write it to `{RESULTS_DIR}/scripts/phaser_count.sh` and substitute the placeholders: `{PHASER_RESOURCES}` from Step 9, and `{PHASED_GT}` and `{PHASER_HOME}` from Step 7.

**Inputs.**
- The BAM is `bam/{sample}.bam`, which is WASP-filtered and duplicate-marked; phASER does no WASP of its own and removes the duplicates itself.
- The VCF is `genotypes/{individual}.het.vcf.gz`, the heterozygous VCF that STAR and ASEReadCounter use. phASER needs the sample name inside that VCF, so the script reads it from the file.

**Flags.** The flags are the tested set:
- `--paired_end 1` (the reason phASER is offered only for paired-end data);
- `--mapq 255` (STAR's unique alignments, the same reads as ASEReadCounter's `--min-mapping-quality` 10 keeps);
- `--baseq {MIN_BASEQ}`;
- `--pass_only 0` (the prepared VCF has FILTER `.`; with the default 1 phASER keeps no site);
- `--unique_ids 1` (its ID column is `.`);
- `--gw_phase_vcf` 1 for phased genotypes, 0 otherwise.
- `--threads 1`, always. **Why:** with the pinned Python 3.14 environment, `phaser.py` with more than one thread stops with `NameError: name 'args' is not defined`, because the default process start method of Python 3.14 no longer passes phASER's settings (a global variable set in `main()`) to its worker processes (seen on the synthetic test, every array task). A pin to Python 3.13 or older, which would start workers the old way, has never been tested. **Consequence:** request one core for this job (`{PHASER_RESOURCES}` holds `-n 1`); each sample runs single-threaded, so the time grows with the sample's reads and cannot be shortened by more cores, only by the array running samples side by side.

**Guards.** The script stops, and leaves no output of the sample, in these cases:
- an installation of another commit, or one whose `install_ok.txt` lacks its second line (`env_sha256 ...`, an interrupted installation); the script only refuses to run, and `setup_phaser_env.sh` decides about a re-installation;
- a contig name containing `_` (`phaser_gene_ae` splits variant IDs on `_`);
- het sites on contigs missing from the BAM header, or no gene span on a BAM contig (mismatched names give silent zero counts in phASER);
- `{PHASED_GT}` = 1 while fewer than 90 % of the het genotypes are phased;
- an empty haplotype table, or a gene table without any gene with counts.

The script loads no module and activates no conda environment: it puts `{PHASER_HOME}/env/bin` first on `PATH`, exports `PYTHONNOUSERSITE=1` and unsets `PYTHONPATH` and `PYTHONHOME` (Step 18, Environment).

```bash
#!/bin/bash
#SBATCH -J ase_phaser
#SBATCH -N 1 -p bcc
#SBATCH --array=1-{ARRAY_N}
#SBATCH {PHASER_RESOURCES}
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={USER_EMAIL}
#SBATCH -o {RESULTS_DIR}/logs/phaser_count_%A_%a.out
set -uo pipefail
# phASER for one sample (array task = data row of samples.csv): read-backed phasing on the WASP-filtered, duplicate-marked BAM
# with the heterozygous VCF that STAR used, then gene-level haplotype counts. No module and no container: the pinned environment
# of setup_phaser_env.sh is used through its bin directory (python, samtools, bcftools and tabix come from it).
PHASER_HOME="{PHASER_HOME}"; PHASER_COMMIT=aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301; PHASED_GT={PHASED_GT}
R="{RESULTS_DIR}"; BED="$R/reference/genes_span.bed"
export PYTHONNOUSERSITE=1      # a per-user site-packages directory of the same Python version would shadow the environment
unset PYTHONPATH PYTHONHOME    # an inherited Python search path would shadow the environment the same way
export PATH="$PHASER_HOME/env/bin:$PATH"; PY="$PHASER_HOME/env/bin/python"; SRC="$PHASER_HOME/src/phaser"
SAMPLE=""
die() { echo "ERROR: sample ${SAMPLE:-?}: $*" >&2; exit 1; }
# install_ok.txt (Step 18): line 1 the pinned commit, line 2 the package-list hash; this script only refuses to run without them
[ "$(sed -n 1p "$PHASER_HOME/install_ok.txt" 2>/dev/null)" = "commit $PHASER_COMMIT" ] && sed -n 2p "$PHASER_HOME/install_ok.txt" 2>/dev/null | grep -qE '^env_sha256 [0-9a-f]{64}$' && [ -x "$PY" ] \
  || die "no complete phASER installation of commit $PHASER_COMMIT in $PHASER_HOME (run setup_phaser_env.sh, Step 18)"
ROW=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1' "{SAMPLES_CSV}" | tr -d '\r')
[ -n "$ROW" ] || die "no row $SLURM_ARRAY_TASK_ID in {SAMPLES_CSV}"
IFS=, read -r SAMPLE FQ1 FQ2 CONDITION CROSS INDIVIDUAL <<< "$ROW"
[ -n "$SAMPLE" ] && [ -n "$INDIVIDUAL" ] || die "empty sample or individual in row $SLURM_ARRAY_TASK_ID of {SAMPLES_CSV}"
BAM="$R/bam/$SAMPLE.bam"; VCF="$R/genotypes/$INDIVIDUAL.het.vcf.gz"
OUT="$R/phaser/$SAMPLE"; LOG="$R/logs/phaser_$SAMPLE.log"; TMPD="$R/tmp/phaser_$SAMPLE"
mkdir -p "$R/phaser" "$TMPD" || die "cannot create output directories"
rm -f "$OUT".*                 # never keep outputs of an earlier run: a failed sample must have none
[ -s "$BAM" ] && [ -s "$BAM.bai" ] || die "missing $BAM or its index (run the per-sample array job first)"
[ -s "$VCF" ] && [ -s "$VCF.tbi" ] || die "missing $VCF or its index (run prep_genotypes.sh first)"
[ -s "$BED" ] || die "missing $BED (written by prep_genotypes.sh)"
NS=$(bcftools query -l "$VCF" | wc -l); VS=$(bcftools query -l "$VCF" | head -1)
[ "$NS" -eq 1 ] || die "$VCF has $NS sample columns; phASER needs exactly one"
# contig names: phaser_gene_ae splits variant IDs on '_' (a contig such as chrUn_xxx breaks it), and names that differ between
# the VCF, the BAM and the gene spans give silent zero counts
N_US=$(bcftools query -f '%CHROM\n' "$VCF" | awk '$1 ~ /_/' | wc -l)
[ "$N_US" -eq 0 ] || die "$N_US heterozygous sites lie on contigs whose names contain '_' ($(bcftools query -f '%CHROM\n' "$VCF" | awk '$1 ~ /_/' | uniq | head -3 | paste -sd' ')); phaser_gene_ae splits variant IDs on '_'. Supply genotype VCFs without these contigs (Step 6)"
BAM_CTG="$TMPD/bam_contigs.txt"
samtools view -H "$BAM" | awk -F'\t' '$1 == "@SQ" {sub(/^SN:/, "", $2); print $2}' > "$BAM_CTG" || die "cannot read the header of $BAM"
N_OFF=$(bcftools query -f '%CHROM\n' "$VCF" | awk 'NR == FNR {c[$1] = 1; next} !($1 in c)' "$BAM_CTG" - | wc -l)
[ "$N_OFF" -eq 0 ] || die "$N_OFF heterozygous sites of $VCF lie on contigs that are not in the BAM header"
N_BED=$(awk 'NR == FNR {c[$1] = 1; next} ($1 in c)' "$BAM_CTG" "$BED" | wc -l)
[ "$N_BED" -gt 0 ] || die "no gene of $BED lies on a contig of the BAM (gene contigs: $(cut -f1 "$BED" | uniq | head -3 | paste -sd' '); BAM contigs: $(head -3 "$BAM_CTG" | paste -sd' '))"
# phased genotypes: --gw_phase_vcf 1 only when the user said so (Step 7) and the VCF really is phased
N_HET=$(bcftools view -H "$VCF" | wc -l); N_PH=$(bcftools query -f '[%GT]\n' "$VCF" | grep -c '|')
if [ "$PHASED_GT" = 1 ]; then
  [ $((N_PH * 10)) -ge $((N_HET * 9)) ] || die "PHASED_GT=1 but only $N_PH of $N_HET heterozygous genotypes in $VCF are phased (0|1 or 1|0); answer 'not phased' in Step 7 or supply phased VCFs"
elif [ $((N_PH * 10)) -ge $((N_HET * 9)) ]; then
  echo "WARNING: $N_PH of $N_HET genotypes are phased but PHASED_GT=0: gene counts use the most-covered block of each gene only (Step 7)"
fi
"$PY" "$SRC/phaser/phaser.py" --vcf "$VCF" --bam "$BAM" --sample "$VS" --paired_end 1 --mapq 255 --baseq {MIN_BASEQ} \
     --pass_only 0 --unique_ids 1 --gw_phase_vcf "$PHASED_GT" --threads 1 --temp_dir "$TMPD" --o "$OUT" > "$LOG" 2>&1 \
  || { tail -5 "$LOG" >&2; rm -f "$OUT".*; die "phaser.py failed (log $LOG)"; }
N_BLK=$(awk 'NR > 1' "$OUT.haplotypic_counts.txt" 2>/dev/null | wc -l)
[ "$N_BLK" -gt 0 ] || { rm -f "$OUT".*; die "phASER wrote no haplotype block (check the BAM, the VCF sample $VS and the contig names; log $LOG)"; }
"$PY" "$SRC/phaser_gene_ae/phaser_gene_ae.py" --haplotypic_counts "$OUT.haplotypic_counts.txt" --features "$BED" --o "$OUT.gene_ae.txt" >> "$LOG" 2>&1 \
  || { tail -5 "$LOG" >&2; rm -f "$OUT".*; die "phaser_gene_ae.py failed (log $LOG)"; }
N_COV=$(awk -F'\t' 'NR > 1 && $7 > 0' "$OUT.gene_ae.txt" 2>/dev/null | wc -l)
[ "$N_COV" -gt 0 ] || { rm -f "$OUT".*; die "no gene has haplotype counts (totalCount > 0): the gene spans and the variants do not overlap (contig names? log $LOG)"; }
N_GW=$(awk -F'\t' 'NR > 1 && $7 > 0 && $11 == 1' "$OUT.gene_ae.txt" | wc -l)
echo "Sample $SAMPLE: $N_BLK haplotype blocks, $N_COV genes with counts ($N_GW genome-wide phased); $(grep -m1 'PHASED' "$LOG")"
rm -rf "$TMPD"
```

| Output | Content |
|---|---|
| `{RESULTS_DIR}/phaser/{sample}.haplotypic_counts.txt` | haplotype blocks with `aCount`, `bCount`, `blockGWPhase` |
| `{RESULTS_DIR}/phaser/{sample}.gene_ae.txt` | one row per gene of `genes_span.bed`: `aCount bCount totalCount n_variants gw_phased` (Rmd 05) |
| `{RESULTS_DIR}/phaser/{sample}.allele_config.txt`, `.variant_connections.txt`, `.haplotypes.txt`, `.allelic_counts.txt`, `.vcf.gz` | other phASER outputs (not used by the Rmds) |
| `{RESULTS_DIR}/logs/phaser_{sample}.log` | phASER and phaser_gene_ae output |

---

## Step 19 — Rmd 05: phASER gene-level haplotype imbalance (outbred, optional)

Only when phASER was selected in Step 7. Write `{CWD}/{TODAY}_{WD_NAME}_05_phaser.Rmd` with the same conventions as Step 12 and the same `{AUTHOR}` and `{PROJECT_TITLE}`.

**Inputs:**
- `ase_checkpoint.rds` (Rmd 01): samples, constants, bias flags;
- `ase_imbalance_checkpoint.rds` (Rmd 02): the unphased gene results;
- `phaser/{sample}.gene_ae.txt` (Step 18), read for every sample of `samples.csv` by name.

It pastes the statistics functions of Step 14 into the chunk marked below.

**Outputs, written to `{RESULTS_DIR}`:**
- `{TODAY}_{WD_NAME}_ASE_phaser.xlsx` (sheets `Gene`, `Comparison`, `Samples`, `Summary`);
- `{TODAY}_{WD_NAME}_ASE_phaser.pdf`;
- `{TODAY}_{WD_NAME}_ASE_phaser_gene.tsv.gz` and `{TODAY}_{WD_NAME}_ASE_phaser_comparison.tsv.gz`: the full `Gene` and `Comparison` tables (tab-separated, gzip). The xlsx sheets hold the full tables while they fit in a worksheet (1,048,575 data rows); beyond that they hold only the significant rows, and the Rmd prints a NOTE;
- `ase_phaser_checkpoint.rds`;
- `summary_numbers_phaser.tsv`.

- **What it adds, and what it does not replace.** The unphased per-SNP tests and the ACAT gene test of Rmd 02 stay the main per-sample result; Rmd 05 is reported next to them, never instead of them.
  - phASER counts each read pair once per haplotype block, so a gene's count pools all the heterozygous SNPs of the block without counting a fragment twice.
  - This gives one effect size per gene: the major haplotype fraction.
  - It can test genes whose SNPs each have fewer than `MIN_DEPTH` reads.
- **Test.** Per sample and gene, `hap_gene_test` (Step 14): an exact two-sided beta-binomial test of `aCount` out of `totalCount` against 0.5.
  - The dispersion comes from the trimmed central-region fit over the sample's genes (`rho_own`); the cohort value (`rho_cohort`) is the median of `rho_own`, and `rho_used = max(rho_cohort, rho_own, RHO_MIN)`.
  - Because `rho_cohort` is a median over the samples of the run, adding or removing samples can change the p-values of the samples already analysed.
  - Any warning of the dispersion fit is printed with its sample and kept in the `rho_warning` column of the `Samples` sheet and of the checkpoint.
  - BH within each sample.
  - `sig` requires `padj < FDR_SIG` and a major haplotype fraction `max(aCount, bCount) / totalCount` at least `0.5 + ABS_DEV_SIG`.
  - Genes with `totalCount` below `MIN_DEPTH` are not tested: `hap_gene_test` does not apply `MIN_DEPTH`, so the Rmd sets these rows to `NA` before the call.
  - Genes that phASER wrote with zero counts (no covered heterozygous SNP; written with `log2_aFC` `inf` and `gw_phased` 1) are dropped.
  - If no sample has 20 genes at `MIN_DEPTH`, the Rmd stops.
  - The Rmd also stops, naming the samples, when a sample's phASER table has no gene with haplotype reads (every `totalCount` 0, or no data rows), when a count field is empty, and when the columns differ from those of `phaser_gene_ae` at the pinned commit. Such a sample is not reported with 0 genes: fix it or drop it (Step 15).
- **Direction.** The test is symmetric in the two haplotypes, because the haplotype labels carry no meaning without phased genotypes.
  - **Unphased genotypes** (`{PHASED_GT}` = 0): A and B are arbitrary per gene and sample. Only the size of the imbalance is reported, with the direction "no direction (haplotype labels arbitrary)". Each gene's counts come from its most-covered haplotype block only, because `phaser_gene_ae` keeps a single block per gene when the genotypes are unphased. The other blocks are dropped, which costs coverage (on the synthetic test, one sample used 113 of 153 heterozygous SNPs).
  - **Phased genotypes** (`{PHASED_GT}` = 1, phASER run with `--gw_phase_vcf 1`): in genome-wide phased genes (`gw_phased` TRUE), haplotype A is the haplotype of the first (left) allele of the phased genotype (`0|1`: A carries the REF allele at that SNP), and the direction is "haplotype A higher" or "haplotype B higher". Genes that are not genome-wide phased keep arbitrary labels and get "no direction (haplotype labels arbitrary)".
- **What comparisons the results support.** Across samples and individuals, only the size of the imbalance can be compared.
  - Haplotype A of one individual is unrelated to haplotype A of another, even with phased genotypes.
  - With unphased genotypes, even two samples of one individual may label a gene's haplotypes differently.
  - Differential ASE on phASER counts is therefore not offered (use Rmd 04).
- **Gene spans.** phASER counts every heterozygous SNP inside the gene span (`reference/genes_span.bed`, introns included), while Rmd 02 uses exonic SNPs only. Overlapping genes share reads. The gene names of the phASER tables must be the GTF `gene_id`s that Rmd 02 uses; if none matches, the Rmd prints a WARNING and the comparison pairs nothing.
- **Comparison table.** Every sample and gene of either analysis gets one row and one category:
  - "both";
  - "phASER only";
  - "phASER only (no SNP tested unphased)": the gene has no row in the gene table of Rmd 02, because no SNP of the gene was tested in Rmd 02 (none reached `MIN_DEPTH`, or none lies in its exons after the Rmd 01 filters), for example a gene whose heterozygous SNPs are all intronic (phASER counts the whole span);
  - "unphased only";
  - "unphased only (not tested by phASER)";
  - "neither".

````rmd
---
title: "{PROJECT_TITLE} - ASE phASER haplotype counts"
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
library(openxlsx)   # no Bioconductor package is needed in this Rmd
library(tidyverse)
```

## Constants and checkpoints

```{r constants}
MODE        <- "{MODE}"
RESULTS_DIR <- "{RESULTS_DIR}"
DATE_TAG    <- "{TODAY}_{WD_NAME}"
MIN_DEPTH   <- {MIN_DEPTH}
FDR_SIG     <- {FDR_SIG}
ABS_DEV_SIG <- {ABS_DEV_SIG}
RHO_MIN     <- {RHO_MIN}             # floor for the dispersion used by the tests
BIAS_TOL    <- {BIAS_TOL}
PHASED_GT   <- {PHASED_GT}           # 1: phased genotype VCFs, phASER ran with --gw_phase_vcf 1 (Step 7); 0: unphased
if (MODE != "outbred") stop("Rmd 05 (phASER) needs MODE = outbred", call. = FALSE)
if (!PHASED_GT %in% c(0, 1)) stop("PHASED_GT must be 0 or 1 (Step 7)", call. = FALSE)
ck <- readRDS(file.path(RESULTS_DIR, "ase_checkpoint.rds"))
same <- c(MODE = identical(ck$constants$MODE, MODE),
          vapply(c("MIN_DEPTH", "FDR_SIG", "ABS_DEV_SIG", "BIAS_TOL"),
                 function(k) isTRUE(all.equal(ck$constants[[k]], get(k))), logical(1)))
if (!all(same)) stop("constants differ from the Rmd 01 checkpoint (", paste(names(same)[!same], collapse = ", "),
                     "); re-render Rmd 01 with the same values", call. = FALSE)
ck2 <- readRDS(file.path(RESULTS_DIR, "ase_imbalance_checkpoint.rds"))   # Rmd 02: the unphased per-SNP and ACAT gene results
if (!isTRUE(all.equal(ck2$constants, ck$constants)))
  stop("the Rmd 02 checkpoint does not come from the current Rmd 01 checkpoint; re-render Rmd 02", call. = FALSE)
bias <- ck$bias
knitr::kable(bias, digits = 4, caption = "Reference-bias diagnostic from Rmd 01 (a screen, not a test)")
if (any(bias$flagged)) cat("FLAGGED samples:", paste(bias$sample[bias$flagged], collapse = ", "), "- read their results with caution\n")
orient_note <- if (PHASED_GT == 1) {
  paste("Phased genotypes: in genome-wide phased genes, haplotype A is the haplotype of the first (left) allele of the phased genotype",
        "(0|1: A carries REF at that SNP); haplotype A of one individual is unrelated to haplotype A of another;",
        "genes that are not genome-wide phased have arbitrary labels.")
} else {
  paste("Unphased genotypes: the haplotype labels A and B are arbitrary per gene and sample, so only the size of the imbalance",
        "(major haplotype fraction) is reported, without a direction; each gene's counts come from its most-covered haplotype block.")
}
cat("NOTE:", orient_note, "\n")
```

## Statistics functions

```{r stats}
# <<< paste here the code of Step 14 between the two marker lines (the marker lines themselves are not pasted) >>>
```

## phASER gene counts

The gene tables are read by sample name from `samples.csv`, never by globbing. Rows with no haplotype reads (phASER writes every gene of the span file, with 0 counts, `log2_aFC` `inf` and `gw_phased` 1 when no heterozygous SNP of the gene is covered) are dropped.

```{r read}
samples <- unique(ck$samples[, c("sample", "condition", "individual")])
f <- file.path(RESULTS_DIR, "phaser", paste0(samples$sample, ".gene_ae.txt"))
miss <- samples$sample[!file.exists(f) | is.na(file.size(f)) | file.size(f) == 0]
if (length(miss) > 0)
  stop("phASER gene counts missing for sample(s) ", paste(miss, collapse = ", "), " (", file.path(RESULTS_DIR, "phaser"),
       "/<sample>.gene_ae.txt): check their phaser_count logs, re-run them, or drop the samples (Step 15)", call. = FALSE)
cls <- c(contig = "character", start = "numeric", stop = "numeric", name = "character", aCount = "numeric", bCount = "numeric",
         totalCount = "numeric", log2_aFC = "character", n_variants = "numeric", variants = "character", gw_phased = "character",
         bam = "character")
ga <- dplyr::bind_rows(lapply(seq_len(nrow(samples)), function(i) {
  d <- read.delim(f[i], colClasses = cls, quote = "", na.strings = character(0))
  if (!identical(names(d), names(cls)))
    stop("unexpected columns in ", f[i], ": found ", paste(names(d), collapse = ", "),
         "; expected (columns of phaser_gene_ae at the commit pinned in Step 18): ", paste(names(cls), collapse = ", "), call. = FALSE)
  dplyr::mutate(d, sample = samples$sample[i])
}))
n_rows <- nrow(ga)
na_cnt <- is.na(ga$aCount) | is.na(ga$bCount) | is.na(ga$totalCount)   # an empty count field is read as NA
if (any(na_cnt))
  stop("missing counts (aCount, bCount or totalCount) in the phASER gene tables: ", sum(na_cnt), " row(s), for example ",
       paste(utils::head(paste0(ga$sample[na_cnt], " / ", ga$name[na_cnt]), 5), collapse = ", "),
       "; re-run phaser_gene_ae for these samples", call. = FALSE)
genes <- ga %>% dplyr::filter(totalCount > 0) %>%   # genes without counts are written with 0 counts, log2_aFC inf and gw_phased 1
  dplyr::transmute(sample, gene_id = name, contig, aCount, bCount, totalCount, n_variants, gw_phased = gw_phased == "1", variants) %>%
  dplyr::left_join(samples, by = "sample")
no_gene <- setdiff(samples$sample, genes$sample)   # every totalCount 0, or a header without data rows
if (length(no_gene) > 0)
  stop("no gene with haplotype reads in the phASER table of sample(s) ", paste(no_gene, collapse = ", "),
       ": every totalCount is 0 or the table has no data rows (for example contig names that differ between the BAM, the VCF and ",
       "genes_span.bed). Check their phaser_count logs, re-run them, or drop the samples (Step 15)", call. = FALSE)
if (any(genes$aCount + genes$bCount != genes$totalCount)) stop("aCount + bCount differs from totalCount in the phASER gene tables", call. = FALSE)
if (anyDuplicated(genes[, c("sample", "gene_id")]))
  stop("a gene name occurs twice in one sample's phASER gene table (gene names in genes_span.bed must be unique)", call. = FALSE)
if (PHASED_GT == 1 && mean(genes$gw_phased) < 0.5)
  cat("WARNING: PHASED_GT = 1 but only", sum(genes$gw_phased), "of", nrow(genes), "genes with counts are genome-wide phased\n")
if (PHASED_GT == 0 && any(genes$gw_phased))
  cat("WARNING: PHASED_GT = 0 but", sum(genes$gw_phased), "genes are genome-wide phased; their labels are treated as arbitrary\n")
cat(n_rows, "sample-gene rows read;", nrow(genes), "with haplotype counts;", sum(genes$totalCount >= MIN_DEPTH), "at or above MIN_DEPTH\n")
```

## Haplotype-level gene test

Per sample and gene, `hap_gene_test` (Step 14) tests the haplotype counts against 0.5. The test is symmetric in the two haplotypes. `hap_gene_test` does not apply `MIN_DEPTH`, so rows below it are set to `NA` before the call and are not tested. The dispersion is the larger of the sample's own trimmed fit and the cohort median, floored at `RHO_MIN`. `rho_cohort` is the median over the samples of this run, so adding or removing samples can change the p-values of the samples already analysed. A warning of the dispersion fit is printed with its sample and kept in `rho_warning`.

```{r test}
collect_warnings <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") })
  list(value = v, warn = paste(unique(w), collapse = "; "))
}
tst <- genes$totalCount >= MIN_DEPTH   # hap_gene_test does not apply MIN_DEPTH: rows below it are NA in the call
ht_w <- character()
withCallingHandlers({
  ht <- hap_gene_test(ifelse(tst, genes$aCount, NA), ifelse(tst, genes$bCount, NA), genes$sample, RHO_MIN)
}, warning = function(cond) { ht_w <<- c(ht_w, conditionMessage(cond)); invokeRestart("muffleWarning") })
ht$samples$rho_warning <- ""
if (length(ht_w) > 0) {   # the dispersion fit warned: re-fit each sample alone (same rows as hap_gene_test) to name the sample
  for (j in seq_len(nrow(ht$samples))) {
    u <- tst & genes$sample == ht$samples$sample[j] & genes$aCount == round(genes$aCount) & genes$bCount == round(genes$bCount)
    if (sum(u) < 20) next
    w <- collect_warnings(bb_estimate_rho_trim(genes$aCount[u], genes$totalCount[u]))$warn
    ht$samples$rho_warning[j] <- w
  }
  if (all(ht$samples$rho_warning == "")) ht$samples$rho_warning <- paste("cohort fit:", paste(unique(ht_w), collapse = "; "))
  for (j in which(ht$samples$rho_warning != ""))
    cat("WARNING (rho estimation), sample", ht$samples$sample[j], ":", ht$samples$rho_warning[j], "\n")
}
if (all(is.na(ht$p)))
  stop("no sample has 20 genes with at least MIN_DEPTH haplotype reads, so no dispersion can be estimated; genes per sample at MIN_DEPTH: ",
       paste(sprintf("%s %d", ht$samples$sample, ht$samples$n_genes), collapse = ", "), call. = FALSE)
genes <- genes %>%
  dplyr::mutate(tested = tst, p = ht$p, rho_used = ht$rho_used, hap_A_frac = aCount / totalCount,
                major_frac = pmax(aCount, bCount) / totalCount, oriented = PHASED_GT == 1 & gw_phased) %>%
  dplyr::group_by(sample) %>% dplyr::mutate(padj = p.adjust(p, method = "BH")) %>% dplyr::ungroup() %>%
  dplyr::mutate(sig = !is.na(padj) & padj < FDR_SIG & major_frac - 0.5 >= ABS_DEV_SIG,
                direction = dplyr::case_when(!sig ~ "none", !oriented ~ "no direction (haplotype labels arbitrary)",
                                             hap_A_frac > 0.5 ~ "haplotype A higher", TRUE ~ "haplotype B higher")) %>%
  dplyr::arrange(sample, gene_id)
knitr::kable(ht$samples, digits = 4, caption = paste("Dispersion per sample: rho_own = trimmed central-region fit over the sample's genes",
             "(NA below 20 genes), rho_cohort = median of rho_own, rho_used = max(rho_cohort, rho_own, RHO_MIN);",
             "rho_warning = warning of the dispersion fit, empty if none"))
knitr::kable(head(dplyr::filter(genes, sig), 40), digits = 4, caption = paste("Significant genes (first 40).", orient_note))
```

## Next to the unphased result of Rmd 02

The unphased ACAT gene test of Rmd 02 stays the main per-sample result; this table puts both results side by side.

```{r compare}
unph <- ck2$gene %>% dplyr::transmute(sample, gene_id, n_snps_unphased = n_snps, acat_p, padj_unphased = padj, sig_unphased = sig)
if (nrow(unph) > 0 && !any(genes$gene_id %in% unph$gene_id))
  cat("WARNING: no phASER gene name is a gene_id of Rmd 02, so the comparison pairs nothing; the gene names of genes_span.bed must be the GTF gene_id\n")
cmp <- genes %>% dplyr::select(sample, gene_id, totalCount, n_variants, gw_phased, major_frac, p, padj, sig) %>%
  dplyr::full_join(unph, by = c("sample", "gene_id")) %>%
  dplyr::mutate(sig = !is.na(sig) & sig, sig_unphased = !is.na(sig_unphased) & sig_unphased,
                category = dplyr::case_when(sig & sig_unphased ~ "both",
                                            sig & is.na(acat_p) ~ "phASER only (no SNP tested unphased)",
                                            sig ~ "phASER only",
                                            sig_unphased & is.na(p) ~ "unphased only (not tested by phASER)",
                                            sig_unphased ~ "unphased only",
                                            TRUE ~ "neither")) %>%
  dplyr::arrange(sample, gene_id)
knitr::kable(dplyr::count(cmp, category), caption = "Sample-gene pairs by category: phASER haplotype test next to the unphased ACAT test of Rmd 02")
```

## Summary

```{r summary}
z0 <- function(x) dplyr::coalesce(as.integer(x), 0L)
summary_ph <- samples %>% dplyr::select(sample) %>%
  dplyr::left_join(dplyr::select(ht$samples, sample, rho_used, rho_own, rho_cohort), by = "sample") %>%
  dplyr::left_join(genes %>% dplyr::group_by(sample) %>%
                     dplyr::summarise(genes_with_counts = dplyr::n(), genes_tested = sum(!is.na(p)), sig_genes = sum(sig),
                                      gw_phased_genes = sum(gw_phased), .groups = "drop"), by = "sample") %>%
  dplyr::left_join(cmp %>% dplyr::group_by(sample) %>%
                     dplyr::summarise(sig_unphased = sum(sig_unphased), sig_both = sum(category == "both"),
                                      sig_phaser_only = sum(startsWith(category, "phASER only")),
                                      sig_phaser_only_untested_unphased = sum(category == "phASER only (no SNP tested unphased)"),
                                      sig_unphased_only = sum(startsWith(category, "unphased only")), .groups = "drop"), by = "sample") %>%
  dplyr::mutate(dplyr::across(c(genes_with_counts, genes_tested, sig_genes, gw_phased_genes, sig_unphased, sig_both, sig_phaser_only,
                                sig_phaser_only_untested_unphased, sig_unphased_only), z0),
                genotypes = if (PHASED_GT == 1) "phased" else "unphased") %>%
  dplyr::select(sample, genotypes, genes_with_counts, genes_tested, sig_genes, gw_phased_genes, sig_unphased, sig_both, sig_phaser_only,
                sig_phaser_only_untested_unphased, sig_unphased_only, rho_used, rho_own, rho_cohort)
knitr::kable(summary_ph, digits = 4, caption = "Summary per sample (FDR_SIG, ABS_DEV_SIG and MIN_DEPTH as in the constants block)")
```

## Figure

```{r figure}
fig <- dplyr::filter(genes, !is.na(p))
fig_n <- vapply(summary_ph$sample, function(s) sum(fig$sig[fig$sample == s]), integer(1))
stopifnot(identical(unname(fig_n), summary_ph$sig_genes))
cat("Figure sig counts equal Summary counts: TRUE\n")
lab <- fig %>% dplyr::group_by(sample) %>%
  dplyr::summarise(lab = paste0(dplyr::first(sample), ": ", sum(sig), " significant of ", dplyr::n(), " genes"), .groups = "drop")
p1 <- dplyr::left_join(fig, lab, by = "sample") %>%
  ggplot(aes(totalCount, major_frac, colour = sig)) + geom_hline(yintercept = 0.5, linetype = 2) + geom_point(alpha = 0.7) +
  scale_x_log10() + scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick")) + facet_wrap(~lab) +
  labs(x = "haplotype reads per gene (log10)", y = "major haplotype fraction (folded, no direction)", colour = "sig",
       caption = paste0("sig: BH < ", FDR_SIG, " and major fraction - 0.5 >= ", ABS_DEV_SIG,
                        if (any(bias$flagged)) "; some samples flagged for reference bias" else "")) + theme_bw()
print(p1)
ggsave(file.path(RESULTS_DIR, paste0(DATE_TAG, "_ASE_phaser.pdf")), p1, width = 10, height = 7)
```

## Export

```{r export}
xlsx_file <- file.path(RESULTS_DIR, paste0(DATE_TAG, "_ASE_phaser.xlsx"))
cut_cell <- function(x) ifelse(nchar(x) > 32000, paste0(substr(x, 1, 32000), " ... (cut: an Excel cell holds at most 32,767 characters; the full list is in ase_phaser_checkpoint.rds)"), x)
# the full tables, whatever their size: tab-separated and gzipped
write_gz <- function(d, name) {
  con <- gzfile(file.path(RESULTS_DIR, paste0(DATE_TAG, name)), "w")
  write.table(d, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)
}
write_gz(genes, "_ASE_phaser_gene.tsv.gz"); write_gz(cmp, "_ASE_phaser_comparison.tsv.gz")
XL_MAX <- 1048575   # data rows per worksheet: Excel holds 1,048,576 rows including the header
xl_rows <- function(d, keep, what) {   # the full table if it fits, otherwise its significant rows
  if (nrow(d) <= XL_MAX) return(d)
  cat("NOTE: the", what, "table has", nrow(d), "rows, more than a worksheet holds; the", what, "sheet holds its", sum(keep),
      "significant rows only; the full table is in", paste0(DATE_TAG, "_ASE_phaser_", tolower(what), ".tsv.gz"), "\n")
  utils::head(d[keep, , drop = FALSE], XL_MAX)
}
wb <- openxlsx::createWorkbook()
for (nm in c("Gene", "Comparison", "Samples", "Summary")) openxlsx::addWorksheet(wb, nm)
openxlsx::writeData(wb, "Gene", dplyr::mutate(xl_rows(genes, genes$sig, "Gene"), variants = cut_cell(variants)))
openxlsx::writeData(wb, "Comparison", xl_rows(cmp, cmp$sig | cmp$sig_unphased, "Comparison"))
openxlsx::writeData(wb, "Samples", ht$samples); openxlsx::writeData(wb, "Summary", summary_ph)
openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)
saveRDS(list(gene = genes, comparison = cmp, samples = ht$samples, summary = summary_ph,
             constants = c(ck$constants, list(RHO_MIN = RHO_MIN, PHASED_GT = PHASED_GT))),
        file.path(RESULTS_DIR, "ase_phaser_checkpoint.rds"))
write.table(summary_ph, file.path(RESULTS_DIR, "summary_numbers_phaser.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", xlsx_file, ", the gene and comparison tables (.tsv.gz) and summary_numbers_phaser.tsv\n")
sessionInfo()
```
````

Render it with `{RESULTS_DIR}/scripts/run_05_phaser.sh`. This is the same script as `run_01_import_qc.sh` (Step 12), with these differences:
- job name `ase_05_phaser`, log `run_05_phaser_%j.out`, the Rmd 05 file name;
- the first case of the `--bind` rule (`--bind {CWD}`): the Rmd reads no GTF;
- `-n 1 --mem=16G -t 4:00:00`. It runs single-threaded, one sample after another.

**Runtime and memory.** Only the test step has been measured, by the unit-test check of Step 14: 55 s and about 2 GB for one sample of 20,000 genes holding about 12 million haplotype reads (`totalCount` summed over its genes), on one core. Everything else below is a projection from that one measurement, not a genome-scale render.
- **Time.** The trimmed dispersion fit enumerates every possible count of every gene, so its time grows in proportion to a sample's haplotype reads: about 4.6 s per million haplotype reads (projected). The samples run one after another, so the test step takes about `4.6 s × (sum over all samples of their haplotype reads in millions)`. For example, 100 samples of 12 million reads take about 1.5 h, and samples of 50 million reads take about 4 minutes each, so the 4 h request fits about 55 of them. Reading, the tables and the export add time that was not measured: keep the projection below about 3 h, and otherwise raise `-t` (for example to twice the projection). Do not split the samples into several runs to save time: `rho_cohort` is the median over the samples of one run, so split runs give different p-values from one joint run.
- **Memory.** The fit of one sample needs about 2 GB per 12 million haplotype reads (projected: about 8 GB at 50 million), and the samples are fitted one after another, so the largest sample sets this part. The gene and comparison tables of all samples are held in memory together, and the xlsx is written from them; that part was not measured. The 16 GB request leaves room up to about 60 million haplotype reads in the largest sample; raise `--mem` beyond that, or when the cohort's gene table is very large.
- **Table size.** The `Gene` and `Comparison` tables have one row per sample and gene (100 samples × 20,000 genes = 2 million rows). The full tables are always written as `{TODAY}_{WD_NAME}_ASE_phaser_gene.tsv.gz` and `{TODAY}_{WD_NAME}_ASE_phaser_comparison.tsv.gz` and kept in the checkpoint. An xlsx sheet holds at most 1,048,575 data rows, so beyond that the `Gene` sheet holds only the significant genes and the `Comparison` sheet only the rows significant in either analysis, and the Rmd prints a NOTE that names the `.tsv.gz` file.

If the job stops with "phASER gene counts missing for sample(s)", the named samples have no phASER output: look at their `phaser_count` logs, re-run them, or drop them (Step 15). If it stops with "the Rmd 02 checkpoint does not come from the current Rmd 01 checkpoint", re-render Rmd 02 first. If it stops with "a gene name occurs twice", the gene span file has duplicated names; rebuild `reference/genes_span.bed` with unique names. If it stops with "no gene with haplotype reads in the phASER table of sample(s)", phASER counted nothing for these samples: compare the contig names of their BAM, the genotype VCF and `genes_span.bed`, look at their `phaser_count` logs, re-run them, or drop them (Step 15). If it stops with "missing counts" or "unexpected columns", the named tables were not written by `phaser_gene_ae` at the pinned commit, or were cut short: re-run `phaser_count` for these samples. Submission order: Step 15.

---

## Notes for the assistant

- **Containers, not modules.** Tools come from the cached Singularity images of Step 2; the only module ever loaded is `singularity/3.10.4`, always with a checked `|| exit 1`. R exists only in the `bulkrnaseq` image (`{R_SIF}`).
- **No R and no heavy work on the login node.** STAR, GATK, Picard, samtools, bcftools, R, the unit tests and the Rmd renders all run as `sbatch -p bcc` jobs; the login node only writes text, reads small files, submits and waits. The summary page is written from the summary TSVs (`summary_numbers.tsv`, and `summary_numbers_reciprocal.tsv` / `summary_numbers_differential.tsv` / `summary_numbers_phaser.tsv` when Rmd 03 / 04 / 05 ran), not from R.
- **`/tmp` is node-local.** Logs, temporary files (`--tmp-dir`, Picard `TMP_DIR`, `_STARtmp`) and every output live under `{RESULTS_DIR}` on the shared filesystem; create `{RESULTS_DIR}/logs` before submitting because `sbatch` silently drops output when the `-o` directory is missing.
- **Every job is submitted with `sbatch -p bcc`**, and chained jobs use `--parsable` and `--dependency=afterok`.
- **Unverified URLs are never embedded.** The five container URLs of Step 2 were verified; Ensembl and Mouse Genomes Project locations are resolved at run time and shown to the user first, never typed from memory.
- **Duplicates are marked, not removed** (`REMOVE_DUPLICATES=false`); ASEReadCounter skips flagged duplicates itself. The read group is written by STAR (`--outSAMattrRGline`) and the het VCF must be bgzipped, tabix-indexed and carry a genotype column, otherwise ASEReadCounter silently counts nothing (the empty-table guard catches it).
- **Orientation.** F1: REF = strain A (the reference assembly's allele) and ALT = strain B; every table and figure says which strain is "higher". Outbred: REF and ALT, gene calls "unphased, no direction".
- **Honesty.** The reference-bias flag is printed next to every ratio table and figure and shown at the top of the summary page. The flag is a screen, not a test. No claim of reduced bias is made for masking or WASP.
- **Statistics block.** The code between `# --- ase-stats-begin` and `# --- ase-stats-end` (Step 14) is pasted verbatim into Rmd 02, 03, 04 and 05; never edit it inside a generated Rmd, because the unit tests check exactly that text.
- **Analyses available.** Per-sample imbalance (Rmd 01/02), reciprocal F1 (Rmd 03) and differential ASE (Rmd 04) are implemented. phASER is a later stage: say "available in a later stage" and continue. Never add a model or package that the Step 14 block does not contain. `lme4` and `aod` are not used by the Rmds (the Step 14 models are base R; Step 2 lists the mandatory packages).
- **Never leave a `{...}` placeholder** in a generated script or Rmd, never push to a remote, and never modify raw FASTQ, BAM or VCF inputs.
