# ase-pipeline Stage 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build Stage 1 of the `ase-pipeline` skill: the wizard, the F1 (third-allele masked reference) and outbred (STAR WASP) upstream scripts, Rmd 01 (import and QC) and Rmd 02 (per-sample allelic imbalance), the static checker, the README, and a synthetic ground-truth acceptance run on the cluster.

**Architecture:** One Markdown wizard (`ase-pipeline/ase-pipeline.md`) that generates sbatch scripts (one prep job, then a per-sample array job) and self-contained Rmds. All tools run from Singularity biocontainers (STAR 2.7.10b, GATK 4.4.0.0, samtools, bcftools, Picard) through wrapper functions with `singularity exec --bind`. Correctness is enforced by (a) a static checker that validates every tool flag against recorded `--help` fixtures, (b) an R unit-test script that extracts the statistics functions from the shipped skill text and tests them, and (c) a synthetic ground-truth dataset with planted imbalance run end to end on the cluster.

**Tech Stack:** Markdown skill file, bash/awk/grep (no python on this host), R in the `bulkrnaseq` image (aod, lme4, openxlsx, tidyverse, GenomicRanges, rtracklayer), STAR 2.7.10b, GATK 4.4.0.0 ASEReadCounter, bcftools 1.20, SLURM.

**Spec:** `docs/superpowers/specs/2026-09-29-ase-pipeline-design.md` (binding). Spec amendments made by this plan are listed in Task 3.

## Global Constraints

- Skill file `ase-pipeline/ase-pipeline.md` is self-contained (no references to other skills' steps); all finite-choice questions are numbered lists; open questions only for email and paths.
- Results go to `results/{YYMMDD}_{WD_NAME}` (same convention as `bulk-rnaseq-pipeline`); raw FASTQ/BAM/VCF inputs are read-only.
- Never run R, Singularity or heavy work on the login node; every job is `sbatch -p bcc`; default requests at most 64 G and 4 h.
- Tools come from Singularity containers, never `module add htslib/bcftools/samtools` (this cluster has no such modules and a failed `module add` silently breaks later `module add` calls); the only module used is `module add singularity/3.10.4` (checked with `|| exit 1` and `command -v singularity`). `BIND` lists each directory once. /tmp is node-local: job outputs go to the shared filesystem.
- Container URLs use the pattern `https://depot.galaxyproject.org/singularity/<name>:<tag>`; every URL embedded in the skill must have been verified with a HEAD request in the session that writes it and recorded in `ase-pipeline/tests/fixtures/verified_urls.txt`.
- Generated scripts use `set -uo pipefail` with explicit exit-code checks (no blind `set -e`); `head -n1` pipes on long streams are avoided.
- STAR index and alignment both use the STAR 2.7.10b container. `sjdbOverhang = read length − 1`, a separate index per read length; `--genomeSAindexNbases` = `min(14, max(4, floor(log2(genome length)/2 − 1)))` computed in shell at run time.
- Duplicates are **marked, not removed** (ASEReadCounter's default `NotDuplicateReadFilter` skips them; verified: 38 marked duplicates lowered the count).
- **Proven by the Task 2 verification gate** (`ase-pipeline/tests/fixtures/verification.md`): (a) STAR must add a read group (`--outSAMattrRGline ID:{sample} SM:{sample} PL:ILLUMINA`), otherwise ASEReadCounter's read-group filter silently drops every read and writes an empty table with exit 0; (b) ASEReadCounter needs a **bgzipped, tabix-indexed** sites VCF with a **heterozygous genotype column** (a sites-only VCF gives 0 rows; homozygous sites are skipped); (c) ASEReadCounter needs `-R` with `.fai` and `.dict` (created by the prep job); (d) STAR WASP does not ignore homozygous genotypes, so the STAR `--varVCFfile` must be heterozygous-only; (e) the reference-fraction diagnostic is judged with a standard-error-based tolerance, not a fixed 0.02.
- Shared constants defined once in each Rmd: `MIN_DEPTH` (default 10), `FDR_SIG` (0.05), `ABS_DEV_SIG` (0.1); tables and figures use only these.
- Rmds are self-contained (no `source()`), `knitr::opts_chunk$set(cache = FALSE)`, `options(scipen = 9)`, Bioconductor packages loaded before tidyverse, `dplyr::`-prefixed verbs where masking is possible.
- Never emit an unverified flag: every `--flag` in the skill text must appear in the recorded fixtures (`star_help.txt`, `gatk_ASEReadCounter_help.txt`) or in the checker's explicit allowlist of other tools' flags.
- Git: show `git status` and `git diff --stat` before each commit; never `git push` without explicit user approval; commit messages end with `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`.

## Review Focus

1. **Third-allele masking**: the substituted base equals neither parental allele and differs from the FASTA base's role; multi-allelic and indel sites are excluded, never half-masked; masked positions equal the SNP count exactly. (Task 5)
2. **WASP filtering**: alignments with no `vW` tag are kept, `vW` 2–7 are dropped, `vW:i:1` are kept; duplicates marked not removed. (Task 5)
3. **Allele-to-strain orientation**: REF↔strain A, ALT↔strain B in every output, verified against planted truth (a swap would invert every imbalance). (Tasks 5, 6, 8)
4. **Beta-binomial edge cases**: ρ → 0 falls back to binomial; zero counts and n < `MIN_DEPTH` are excluded; `optim` non-convergence is reported, not ignored. (Task 6)
5. **Wizard constraints**: reciprocal/differential/phASER options are only offered when their preconditions hold; sample names with dashes or spaces are rejected. (Tasks 3, 4)
6. **Array-job indexing**: `--array=1-N` matches the sample list length; one failing sample fails loudly and never yields an empty count table that Rmd 01 accepts. (Task 5)

## File Structure

| File | Responsibility |
|---|---|
| `ase-pipeline/ase-pipeline.md` | The wizard prompt (Steps 0–13, Notes). Created Task 3, extended Tasks 4–6. |
| `ase-pipeline/tests/check_skill.sh` | Static checker. Created Task 3, extended Tasks 4–7. |
| `ase-pipeline/tests/fixtures/` | Recorded `--help` outputs, `verification.md`, `verified_urls.txt`. Task 2 and later. |
| `ase-pipeline/tests/synthetic/simulate_ase_data.R` | Ground-truth data generator (F1 and outbred). Task 1. |
| `ase-pipeline/tests/r/test_ase_stats.R` | Unit tests for the statistics functions extracted from the skill text. Task 6. |
| `ase-pipeline/README.md`, root `README.md` | Docs and skills-table row. Task 7. |
| `~/.claude/commands/ase-pipeline.md` | Installed copy (Task 8). |

All repo paths are relative to `/net/bmc-lab3/data/bcc/yannvrb/Github_repos/bioinformatics-claude-skills`. The synthetic base data is the nf-core rnavar test genome already on disk in `/net/bmc-lab3/data/bcc/yannvrb/rnavar_test/` (`genome.fasta`, `genome.gtf`, 40 kb of chr22 in fragment-local coordinates).

---

### Task 1: Synthetic ground-truth data generator

**Files:**
- Create: `ase-pipeline/tests/synthetic/simulate_ase_data.R`

**Interfaces:**
- Produces: `Rscript simulate_ase_data.R <genome.fasta> <genome.gtf> <outdir> <seed>` writes, under `<outdir>/f1/` and `<outdir>/outbred/`: paired-end FASTQs `{sample}_1.fastq.gz`/`_2.fastq.gz` (read length 100), `samples.csv`, `truth_genes.tsv` (columns `gene_id, n_snps, p_alt, class`), `truth_snps.tsv` (`chrom, pos, ref, alt, gene_id`), plus for F1 `parental_snps.vcf` (sites VCF, REF = strain A allele, ALT = strain B allele, INFO-free) and for outbred one single-sample VCF per individual `{individual}.het.vcf` (GT `0/1`) and `{individual}.all.vcf` (all genotypes incl. `0/0`/`1/1`). Later tasks depend on these names. **Synthetic genome mode:** if `<genome.fasta>` and `<genome.gtf>` are both the literal `SYNTHETIC`, the script first generates a deterministic synthetic genome (one contig `chr1`, 300,000 bp of random sequence with about 45 percent GC; 20 non-overlapping `+`-strand genes with 3–5 exons of 120–300 bp and introns of 300–2000 bp; GTF lines `gene`, `transcript` and `exon` carrying `gene_id` and `transcript_id`) into `<outdir>/genome/genome.fa` and `<outdir>/genome/genome.gtf` and uses it. This mode is the default for all acceptance work, because the nf-core test genome has only 6 usable non-overlapping loci (it stays supported for smoke tests with quotas scaled down; `check_simulation.sh` reads `MIN_NULL` from the environment, default 4).

- [ ] **Step 1: Write the generator test first**

Create `ase-pipeline/tests/synthetic/check_simulation.sh`:

```bash
#!/bin/bash
# Usage: check_simulation.sh <outdir>   (run after simulate_ase_data.R)
set -u; D=${1:?usage: check_simulation.sh <outdir>}; fail=0
chk() { [ -s "$1" ] || { echo "FAIL: missing/empty $1"; fail=1; }; }
for m in f1 outbred; do
  chk $D/$m/samples.csv; chk $D/$m/truth_genes.tsv; chk $D/$m/truth_snps.tsv
  for s in $(tail -n +2 $D/$m/samples.csv | cut -d, -f1); do chk $D/$m/${s}_1.fastq.gz; chk $D/$m/${s}_2.fastq.gz; done
done
chk $D/f1/parental_snps.vcf
for i in $(tail -n +2 $D/outbred/samples.csv | cut -d, -f6 | sort -u); do chk $D/outbred/$i.het.vcf; chk $D/outbred/$i.all.vcf; done
# truth: at least 2 planted imbalanced genes and 4 null genes; alt fraction values in [0,1]
awk -F'\t' 'NR>1{c[$4]++; if($3<0||$3>1) bad=1} END{ if(c["imbalanced"]<2||c["null"]<4||bad) exit 1 }' $D/f1/truth_genes.tsv \
  || { echo "FAIL: f1 truth_genes.tsv needs >=2 imbalanced and >=4 null genes with p_alt in [0,1]"; fail=1; }
# reads: 4 lines per record, equal read counts in mates, read length 100
for f in $D/f1/*_1.fastq.gz; do
  n1=$(zcat $f | awk 'END{print NR/4}'); n2=$(zcat ${f%_1.fastq.gz}_2.fastq.gz | awk 'END{print NR/4}')
  [ "$n1" = "$n2" ] && [ "$n1" -ge 500 ] || { echo "FAIL: $f read counts ($n1 vs $n2)"; fail=1; }
  zcat $f | awk 'NR%4==2{print length($0)}' | sort -u | grep -qx 100 || { echo "FAIL: $f read length not 100"; fail=1; }
done
# determinism marker
[ -s $D/seed.txt ] || { echo "FAIL: missing seed.txt"; fail=1; }
[ $fail -eq 0 ] && echo PASS || exit 1
```
Run it on an empty directory. Expected: FAIL (files missing).

- [ ] **Step 2: Write the generator**

Implement `simulate_ase_data.R` (R only; uses `Biostrings`, `rtracklayer`, base R) with this algorithm, seeded by the fourth argument (write it to `seed.txt`):
  1. Read the FASTA (`Biostrings::readDNAStringSet`) and the GTF (`rtracklayer::import`). Group exons per gene using the first transcript. Keep genes whose spliced exon length is at least 400 bp; take up to 12 genes.
  2. For each kept gene choose `k` SNP positions inside exons (k = 4 for ordinary genes; **one designated "bias" gene gets 5 SNPs within a 60 bp window** so an alt-allele read carries several mismatches and maps worse than a ref read), spaced at least 40 bp apart otherwise. The alt base is a base different from the reference base (seeded choice). Record `truth_snps.tsv`.
  3. Classes for `truth_genes.tsv`: `imbalanced` = 2 genes with `p_alt` 0.70, 1 gene 0.85, 1 gene 0.30; `bias` = the dense gene with `p_alt` 0.5; `null` = all remaining genes with `p_alt` 0.5. `p_alt` is the fraction of the gene's fragments drawn from the ALT haplotype (F1: strain B).
  4. **F1 design:** six samples `f1_s1..f1_s6` (three replicates each of two conditions `ctrl`/`treat`; `cross_direction` = `AxB` for all; no planted parent-of-origin effect in Stage 1). Per sample and gene draw the number of fragments from `Poisson(450)`, split ALT/REF by `Binomial(n, p_alt)`; the same planted p_alt in every sample. `parental_snps.vcf` lists every SNP (`REF`= strain A allele, `ALT` = strain B allele, `##INFO`-free, sorted).
  5. **Outbred design:** four individuals `ind1..ind4`, one sample each (`outbred_s1..s4`; two conditions of two individuals). At each SNP an individual's genotype is `0/1` with probability 0.6, else `0/0` or `1/1` with equal probability; only heterozygous SNPs contribute allelic imbalance (`p_alt` planted as above for genes with at least one het SNP; homozygous SNPs give all-REF or all-ALT reads). Write `{individual}.all.vcf` (single sample column, genotype per SNP) and `{individual}.het.vcf` (only `0/1`).
  6. Fragments: a fragment start is uniform over the spliced transcript, length `round(rnorm(1, 250, 30))` clipped to [150, 350] and to the transcript length; the read pair is the first and last 100 bp (mate 2 reverse-complemented); the spliced transcript sequence is built per haplotype by substituting the haplotype's SNP alleles; add substitution errors at rate 0.002; quality string `I` × 100. Names `@{sample}_{i}/1` and `/2`.
  7. Write gzip FASTQs and `samples.csv` with columns `sample,fastq_1,fastq_2,condition,cross_direction,individual` (`cross_direction` empty for outbred, `individual` empty for F1).

- [ ] **Step 3: Run the generator and the test on the cluster**

Run via `sbatch -p bcc` with the `bulkrnaseq` image (`singularity exec --bind /net/bmc-lab3 /net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif Rscript ...`; job output on the shared filesystem, e.g. `/net/bmc-lab3/data/bcc/yannvrb/ase_synthetic/`):
`Rscript ase-pipeline/tests/synthetic/simulate_ase_data.R SYNTHETIC SYNTHETIC /net/bmc-lab3/data/bcc/yannvrb/ase_synthetic 20260929`
then `bash ase-pipeline/tests/synthetic/check_simulation.sh /net/bmc-lab3/data/bcc/yannvrb/ase_synthetic`. Expected: `PASS`. Also run the generator a second time with the same seed into a different directory and confirm the FASTQ md5 sums are identical (determinism).

- [ ] **Step 4: Commit**

```bash
git status && git diff --stat
git add ase-pipeline/tests/synthetic
git commit -m "ase-pipeline: synthetic ground-truth data generator (F1 and outbred)"
```
(Do not commit generated data; only the scripts.)

---

### Task 2: Verification gate (controller-run on the cluster)

Run by the controller (needs cluster jobs and judgment about results); records facts that later tasks depend on. **Files:** create `ase-pipeline/tests/fixtures/{star_help.txt,gatk_ASEReadCounter_help.txt,r_packages.txt,star_version.txt,verification.md,verified_urls.txt}`.

- [ ] **Step 1: Copy the already-captured fixtures**

Copy `/net/bmc-lab3/data/bcc/yannvrb/rnavar_test2/ase_fixtures/{star_help.txt,gatk_ASEReadCounter_help.txt,r_packages.txt,star_version.txt}` into `ase-pipeline/tests/fixtures/`.

- [ ] **Step 2: Verify container URLs**

`curl -sI` each container the skill will use and append `<url> <HTTP status> <content-length> <date>` to `verified_urls.txt`: `star:2.7.10b--h9ee0642_0`, `gatk4:4.4.0.0--py36hdfd78af_0`, `bcftools:1.20--h8b25389_0`, `samtools:1.21--h50ea8bc_0`, `picard:3.1.1--hdfd78af_0` under `https://depot.galaxyproject.org/singularity/`. Expected: HTTP 200 for each; any other status means that tool needs another tag (find it in `~/.singularity/cache` names) before continuing.

- [ ] **Step 3: STAR WASP behaviour test**

On the synthetic outbred data: build a STAR index for `/net/bmc-lab3/data/bcc/yannvrb/ase_synthetic/genome/genome.fa` with its GTF (container 2.7.10b, `--genomeSAindexNbases 6`, `--sjdbOverhang 99`), align one outbred sample with `--varVCFfile <individual>.all.vcf --waspOutputMode SAMtag --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM vA vG vW`, then with `<individual>.het.vcf`. Record in `verification.md`: (a) whether reads overlapping homozygous SNPs receive a `vW` tag when the all-genotypes VCF is used (i.e. whether homozygous sites are ignored); (b) counts of `vW:i:1`, `vW:i:2..7` and untagged reads; (c) the mean ALT fraction at the bias gene's SNPs before and after removing `vW`≠1 reads; (d) that a multi-sample VCF uses only its first sample (test with a two-sample VCF built by `bcftools merge`).

- [ ] **Step 4: ASEReadCounter behaviour and output columns**

On the BAM from Step 3 run `gatk ASEReadCounter` (container 4.4.0.0) with `--variant` = the het VCF, then with the all-genotypes VCF, then with a sites-only VCF (no genotype column). Record: which variants are counted in each case (does the tool require a het genotype?), the exact output column names, and the effect of `--min-mapping-quality 10` and `--min-base-quality 10`. This decides whether the F1 prep must write a het-sites VCF (`GT 0/1`) beside `parental_snps.vcf`.

- [ ] **Step 5: Mouse Genomes Project region query**

From a compute node: `bcftools view -r 19:1000000-1100000 -s A_J,C57BL_6NJ https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz` (bcftools 1.20 container) with a timer. Record whether region queries work through the `.csi` index remotely, the elapsed time, and the estimated time for the whole genome (sizes the mouse helper's sbatch request and decides remote streaming versus download).

- [ ] **Step 6: Write `verification.md` and commit**

`verification.md` lists each finding with the command and its output snippet, and states the resulting design decisions (any change to the spec is made in Task 3). Commit `git add ase-pipeline/tests/fixtures && git commit -m "ase-pipeline: verification gate results and recorded tool help fixtures"`.

---

### Task 3: Checker scaffold and wizard Steps 0–3 (plus spec amendments)

**Files:**
- Create: `ase-pipeline/tests/check_skill.sh`, `ase-pipeline/ase-pipeline.md`
- Modify: `docs/superpowers/specs/2026-09-29-ase-pipeline-design.md`

**Interfaces:**
- Produces: `check_skill.sh <skill.md>` with helpers `need "<text>"` / `forbid "<text>"` and a flag rule; skill Steps 0–3 defining `{CWD}`, `{WD_NAME}`, `{TODAY_YYMMDD}`, `{RESULTS_DIR}`, `{USER_EMAIL}`, `{MODE}` ∈ {`f1`, `outbred`} and the container variables `{SIF_DIR}`, `{STAR_SIF}`, `{GATK_SIF}`, `{BCFTOOLS_SIF}`, `{SAMTOOLS_SIF}`, `{PICARD_SIF}`.

- [ ] **Step 1: Write the checker first**

`check_skill.sh` mirrors `nfcore-rnavar-setup/tests/check_skill.sh` (`need`, `forbid`, exit codes) with this flag rule instead of the nf-core schema:

```bash
FIX=$(dirname "$0")/fixtures
tool_flags=$( { grep -ohE '^ *--[A-Za-z][A-Za-z0-9-]*|^--[A-Za-z][A-Za-z0-9-]*|^[A-Za-z]+ +.*' "$FIX/gatk_ASEReadCounter_help.txt" | grep -oE '^ *--[A-Za-z][A-Za-z0-9-]*' ; \
                grep -oE '^[A-Za-z][A-Za-z0-9]+' "$FIX/star_help.txt"; } | tr -d ' ' | sed 's/^--//' | sort -u )
allow="runMode genomeDir genomeFastaFiles sjdbGTFfile sjdbOverhang runThreadN genomeSAindexNbases readFilesIn readFilesCommand outFileNamePrefix outSAMtype outSAMattributes outSAMunmapped twopassMode varVCFfile waspOutputMode bind mail-type mail-user mem array output-fmt regions targets exclude-uncalled genotype types min-alleles max-alleles samples apply-filters include ref-allele-fasta reference fasta-ref file"
for flag in $(grep -oE '(^|[ `(=])--[A-Za-z][A-Za-z0-9-]*' "$SKILL" | sed -E 's/^[^-]*--//' | sort -u); do
  echo "$tool_flags $allow" | tr ' ' '\n' | grep -qx -- "$flag" || { echo "FAIL: flag not in fixtures or allowlist: --$flag"; fail=1; }
done
```
(Adjust the extraction so `tool_flags` contains every long option name that appears in the two help fixtures; STAR parameters are listed one per line at line start in `star_help.txt`, GATK options as `--name,-x <Type>` at line start. The implementer verifies the extraction by printing `tool_flags` and confirming it contains `waspOutputMode`, `varVCFfile`, `min-base-quality`, `count-overlap-reads-handling`.) Also add: `need "# ase-pipeline — allele-specific expression"`, `need "## Step 0"`, `need "## Step 3"`, `forbid "module add htslib"`, `forbid "module add bcftools"`, `need "singularity exec --bind"`, `need "|| { echo"` (checked module add). Run on the missing file: expected FAIL.

- [ ] **Step 2: Write skill Steps 0–3**

Create `ase-pipeline/ase-pipeline.md` with a title line `# ase-pipeline — allele-specific expression (F1 cross and outbred/human)`, an introduction (purpose, the two modes, the honesty rule that ratios are never reported without the reference-bias diagnostic), and:
  - **Step 0** working directory, `{WD_NAME}`, `{TODAY_YYMMDD}`, `{RESULTS_DIR}` = `{CWD}/results/{TODAY_YYMMDD}_{WD_NAME}`.
  - **Step 1** email (no pre-fill).
  - **Step 2** environment check, no R on the login node: (a) the Singularity module load `module add singularity/3.10.4 || { echo "ERROR: cannot load singularity" >&2; exit 1; }` with `command -v singularity`; (b) locate the cached containers in `${NXF_SINGULARITY_CACHEDIR:-$HOME/.singularity/cache}` by the file names `depot.galaxyproject.org-singularity-<name>-<tag>.img`, and for any missing one the helper downloads it (`wget -c -O X.part URL && mv`), URLs from `verified_urls.txt`; (c) the R package check of the `bulkrnaseq` image (`aod`, `lme4`, `openxlsx`, `tidyverse`, `GenomicRanges`, `rtracklayer` mandatory) run as a tiny `sbatch` job, never on the login node; stop and ask for an image rebuild if a mandatory package is missing.
  - **Step 3** mode question (numbered: 1. F1 cross · 2. Outbred/human) with a two-line explanation of each and what it means for genotypes and mapping; store `{MODE}`.

- [ ] **Step 3: Amend the spec (same commit)**

In `docs/superpowers/specs/2026-09-29-ase-pipeline-design.md` make these edits, each as a sentence replacement: (1) the sample sheet is `{WD_NAME}_samples.csv` (not xlsx) because no R runs on the login node; (2) in the F1 upstream section add: "the reference-prep job also writes `f1_het_sites.vcf.gz` (single sample `F1`, genotype `0/1` at every parental SNP) because ASEReadCounter counts at heterozygous sites; the parental-difference VCF stays the source of truth" — (Task 2 confirmed: a sites-only VCF gives 0 rows, so this file is required); (3) in the architecture section state that STAR runs from the cached 2.7.10b container for both index and alignment.

- [ ] **Step 4: Run the checker; expected PASS; commit**

```bash
bash ase-pipeline/tests/check_skill.sh ase-pipeline/ase-pipeline.md
git status && git diff --stat
git add ase-pipeline docs/superpowers/specs/2026-09-29-ase-pipeline-design.md
git commit -m "ase-pipeline: checker scaffold, wizard Steps 0-3, spec amendments"
```

---

### Task 4: Wizard Steps 4–9 (reference, samples, genotypes, analyses, constants, resources)

**Files:** Modify `ase-pipeline/tests/check_skill.sh`, `ase-pipeline/ase-pipeline.md`.

**Interfaces:**
- Consumes: `{MODE}`, `{RESULTS_DIR}`, container variables (Task 3).
- Produces: `{FASTA_PATH}`, `{GTF_PATH}`, `{GENOME_DIR}`, `{READ_LENGTH}`, `{SJDB_OVERHANG}`, `{STAR_INDEX}`, `{SAMPLES_CSV}`, `{PARENTAL_VCF}` or `{GENOTYPE_VCFS}`, `{STRAIN_A}`, `{STRAIN_B}`, `{ANALYSES}`, `MIN_DEPTH`, `FDR_SIG`, `ABS_DEV_SIG`, `{ARRAY_N}`, `{MEM}`, `{TIME}`.

- [ ] **Step 1: Failing checker lines**

Add `need` lines for the distinctive phrases you will write: `need "## Step 4"` … `need "## Step 9"`, `need "star_ase_sjdb"` (index folder name per read length), `need "genomeSAindexNbases"`, `need "{SAMPLES_CSV}"`, `need "cross_direction"`, `need "individual"`, `need "MIN_DEPTH"`, `need "FDR_SIG"`, `need "ABS_DEV_SIG"`, `need "reciprocal"`, `need "phASER"`, `need "--array"`, `need "-t 4:00:00"`, and `forbid "strandedness"`. Run: expected FAIL.

- [ ] **Step 2: Write Steps 4–9**

  - **Step 4 Reference and read length:** FASTA and GTF (custom paths, or the folder convention `{genome_base}/{organism}/{assembly}_ens{version}/` used by the other skills; ask which); gzipped inputs are decompressed by the prep job with `gunzip -c file.gz > {GENOME_DIR}/<name>` (originals untouched); read-length detection with the same `awk` one-liner as the other skills, mode of lengths, `{SJDB_OVERHANG}`; F1 mode index folder `{GENOME_DIR}/index/star_ase_masked_sjdb{N-1}/`, outbred `{GENOME_DIR}/index/star_ase_sjdb{N-1}/` (F1 uses the masked genome, so its index never collides with an unmasked one).
  - **Step 5 Sample sheet:** scaffold `{WD_NAME}_samples.csv` from the FASTQs (paired-end detection, name sanitisation exactly as in the other skills: `-` and spaces to `_`, collisions warn); columns `sample,fastq_1,fastq_2,condition,cross_direction,individual` (F1: `cross_direction` and the two strain names asked once; outbred: `individual` per row). Validation rules: no dashes or spaces; each condition needs at least one replicate (warn); reciprocal analysis needs both cross directions.
  - **Step 6 Genotype source:** F1: `{PARENTAL_VCF}` (user-supplied biallelic parental-difference VCF where REF is strain A and ALT is strain B) or the mouse helper (Task 5); outbred: one VCF per `individual` from external WGS/array data or the rnavar `variant_calling/` filtered VCFs; state the RNA-derived-genotype caveat verbatim from the spec whenever the rnavar VCF is chosen.
  - **Step 7 Analysis menu:** numbered multi-select; per-sample imbalance always; reciprocal only when `{MODE}` = `f1` and both `cross_direction` values are present; differential only with at least two conditions with replicates; phASER only when `{MODE}` = `outbred`; Stage 1 implements only the always-on analysis, and choosing the others must say "available in a later stage" and continue.
  - **Step 8 Constants:** `MIN_DEPTH` 10, `FDR_SIG` 0.05, `ABS_DEV_SIG` 0.1, plus ASEReadCounter `min-mapping-quality` 10 and `min-base-quality` 10 defaults; shown and editable.
  - **Step 9 Resources:** per-sample array job `#SBATCH --array=1-{ARRAY_N}` (N = number of samples), STAR alignment human default `-n 8 --mem=48G -t 4:00:00`, small genomes `-n 4 --mem=8G -t 1:00:00`, prep jobs scaled by genome size as in the other skills; never above 64 G / 4 h without the user asking.

- [ ] **Step 3: Run the checker; expected PASS; commit**

```bash
bash ase-pipeline/tests/check_skill.sh ase-pipeline/ase-pipeline.md
git status && git diff --stat
git add ase-pipeline && git commit -m "ase-pipeline: wizard Steps 4-9 (reference, samples, genotypes, analyses, constants, resources)"
```

---

### Task 5: Upstream script templates (Steps 10–11)

**Files:** Modify `ase-pipeline/tests/check_skill.sh`, `ase-pipeline/ase-pipeline.md`.

**Interfaces:**
- Consumes: Task 4 variables.
- Produces (generated files, named `{RESULTS_DIR}/scripts/`): `prep_f1_reference.sh` (F1), `prep_genotypes.sh` (outbred), `align_count_f1.sh`, `align_wasp_count.sh` (array jobs), `extract_mgp_parental_vcf.sh` (optional). Outputs per sample: `{RESULTS_DIR}/bam/{sample}.bam` (+`.bai`), `{RESULTS_DIR}/ase_counts/{sample}.table`. Rmd 01 (Task 6) reads `ase_counts/*.table` and `samples.csv`.

- [ ] **Step 1: Failing checker lines**

`need` for: the container wrapper block (`singularity exec --bind "$BIND" "$SIF"`), `bcftools consensus`, `third allele`, `--waspOutputMode SAMtag`, `--varVCFfile`, `vW`, `MarkDuplicates`, `ASEReadCounter`, `--min-mapping-quality`, `--min-base-quality`, `--count-overlap-reads-handling`, `SLURM_ARRAY_TASK_ID`, `set -uo pipefail`, an existence check of every `ase_counts/{sample}.table` (`needs non-empty`), and `forbid` for `--outSAMattributes All` (WASP needs explicit attributes) and `module add bcftools`. Run: expected FAIL.

- [ ] **Step 2: F1 prep (`prep_f1_reference.sh`) text and snippets**

The skill instructs the wizard to write a prep script whose logic is:
  1. Filter `{PARENTAL_VCF}` to biallelic SNPs: `bcftools view -m2 -M2 -v snps`; count sites.
  2. For each site pick the third base with this function (same rule everywhere):
     ```bash
     third() { for c in A C G T; do if [ "$c" != "$1" ] && [ "$c" != "$2" ]; then echo "$c"; return; fi; done; }
     ```
     write `masked_sites.vcf` (`REF` = the base in the FASTA, read with `samtools faidx`; `ALT` = `third $REF_STRAIN_A $ALT_STRAIN_B`) and record `masked_positions.bed`.
  3. `bcftools consensus -f {FASTA_PATH} masked_sites.vcf.gz > masked.fa` (bgzip and `tabix -f -p vcf` first), then verify: number of positions where `masked.fa` differs from the FASTA equals the number of sites, and at every site the masked base equals the chosen third allele; the job exits 1 otherwise.
  4. Write `f1_het_sites.vcf.gz` (single sample `F1`, genotype `0/1` at every parental SNP, bgzipped and `tabix -f -p vcf` indexed): required, because ASEReadCounter counts only heterozygous genotype sites and a sites-only VCF gives 0 rows. Also create the reference files ASEReadCounter needs: `samtools faidx` (`.fai`) and `gatk CreateSequenceDictionary` (`.dict`) next to the FASTA used for counting (the original FASTA, not the masked one).
  5. Build the STAR index from `masked.fa` (`--genomeSAindexNbases` computed in shell, `--sjdbOverhang {SJDB_OVERHANG}`, container 2.7.10b).

- [ ] **Step 3: Outbred genotype prep (`prep_genotypes.sh`)**

Per individual: `bcftools view -s <sample> -g het -v snps -m2 -M2` (and, for rnavar VCFs, `-f PASS` plus `-i 'INFO/DP>=10'` style filters that the wizard shows and lets the user edit; for external genotypes a `GQ` filter) into `{individual}.het.vcf.gz` (single sample, `tabix -f -p vcf`); print the number of sites per individual and warn when fewer than 1000 (human) or fewer than 20 (test data) heterozygous sites remain.

- [ ] **Step 4: Per-sample array scripts**

Both scripts start `set -uo pipefail`, take the sample row from `{SAMPLES_CSV}` by `SLURM_ARRAY_TASK_ID`, and check every step's exit code. F1 (`align_count_f1.sh`): STAR (container) to the masked index with `--outSAMattrRGline ID:{sample} SM:{sample} PL:ILLUMINA` (mandatory: without a read group ASEReadCounter silently counts nothing) → `samtools sort`/`index` → Picard `MarkDuplicates` (mark, not remove; the Picard step must keep the read group) → `gatk ASEReadCounter -R {FASTA_PATH} -I bam -V f1_het_sites.vcf.gz --min-mapping-quality M --min-base-quality Q --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE -O ase_counts/{sample}.table`. Outbred (`align_wasp_count.sh`): STAR with `--varVCFfile {individual}.het.vcf` (the heterozygous-only VCF; STAR does not ignore homozygous genotypes) `--waspOutputMode SAMtag --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM vA vG vW --outSAMattrRGline ID:{sample} SM:{sample} PL:ILLUMINA`, keep alignments with `vW:i:1` and alignments with no `vW` tag (drop `vW` 2–7) with `samtools view -h` plus `awk`, index, `MarkDuplicates`, ASEReadCounter as above with `-V {individual}.het.vcf.gz` (bgzipped and `tabix -f -p vcf` indexed: ASEReadCounter fails on a plain VCF). Each script also writes `ase_counts/{sample}.wasp_stats.tsv` (outbred) recording, from the STAR BAM before filtering, the number of alignments by `vW` code (`none`, 1–7) and by `vA` allele class (1 = ref, 2 = alt, 3 = no match), so Rmd 01 can show how many reads WASP removed and from which allele. After ASEReadCounter, each script fails (exit 1, message naming the sample) if `ase_counts/{sample}.table` is missing or has no data rows, so an empty table can never reach Rmd 01.

- [ ] **Step 5: Mouse helper (`extract_mgp_parental_vcf.sh`)**

Text and snippet only (the helper was sized in Task 2 Step 5): for the two named strains (`{STRAIN_A}` = the reference strain, e.g. `C57BL_6NJ`, and `{STRAIN_B}`, e.g. `A_J`; names exactly as in the VCF header) run **one array task per chromosome** (a whole-genome serial run would take about 6 hours: remote access has about 80 s fixed overhead and about 80 s per 10 Mb of region) with `bcftools view -r <chr> -s A,B -f PASS -m2 -M2 -v snps`, keep sites where the two genotypes differ and both have `FMT/FI=1`, and concatenate the per-chromosome results with `bcftools concat` in a final dependent job. MGP contigs are `1`, `2` ... (no `chr` prefix), so the FASTA must use the same names. The wizard confirms the URL with a HEAD request and states that the release is GRCm39, so the FASTA must be GRCm39 too (contig names checked with the same first-contig guard used in `nfcore-rnavar-setup`).

- [ ] **Step 6: Run the checker; expected PASS; commit**

```bash
bash ase-pipeline/tests/check_skill.sh ase-pipeline/ase-pipeline.md
git status && git diff --stat
git add ase-pipeline && git commit -m "ase-pipeline: upstream script templates (F1 masked reference, WASP, array jobs, mouse helper)"
```

---

### Task 6: Rmd 01 and Rmd 02 templates, statistics functions and their unit tests

**Files:** Create `ase-pipeline/tests/r/test_ase_stats.R`; modify `ase-pipeline/tests/check_skill.sh`, `ase-pipeline/ase-pipeline.md`.

**Interfaces:**
- Consumes: `ase_counts/{sample}.table` (columns recorded in Task 2 Step 4), `samples.csv`.
- Produces: the statistics function block in the skill between the marker lines `# --- ase-stats-begin` and `# --- ase-stats-end` (functions `bb_negloglik`, `bb_estimate_rho`, `bb_pvalue`, `acat`), Rmd 01 (`{DATE}_{WD_NAME}_01_import_qc.Rmd`) writing `ase_checkpoint.rds`, Rmd 02 (`..._02_imbalance.Rmd`) writing `*_ASE_imbalance.xlsx`.

- [ ] **Step 1: Write the unit-test script first**

`ase-pipeline/tests/r/test_ase_stats.R` extracts the functions from the shipped skill text so the tests exercise the shipped code:

```r
args <- commandArgs(trailingOnly = TRUE); skill <- args[1]
txt <- readLines(skill)
b <- grep("^# --- ase-stats-begin", txt); e <- grep("^# --- ase-stats-end", txt)
stopifnot(length(b) == 1, length(e) == 1, e > b)
eval(parse(text = txt[(b + 1):(e - 1)]))
set.seed(1); fail <- 0
ok <- function(cond, msg) { if (!isTRUE(cond)) { cat("FAIL:", msg, "\n"); fail <<- 1 } else cat("ok  ", msg, "\n") }

# 1. null calibration: p-values approximately uniform under H0 (rho = 0.02, n = 50)
rho0 <- 0.02; n <- rep(50, 3000); a <- 0.5 * (1 - rho0) / rho0
x <- rbinom(3000, n, rbeta(3000, a, a))
p <- mapply(bb_pvalue, x, n, MoreArgs = list(rho = rho0))
ok(abs(mean(p < 0.05) - 0.05) < 0.02, sprintf("null size at 0.05 = %.3f", mean(p < 0.05)))
ok(ks.test(p, "punif")$p.value > 0.001 || abs(mean(p) - 0.5) < 0.05, "null p-values roughly uniform")

# 2. power: alt fraction 0.7, n = 60, low overdispersion
x1 <- rbinom(500, 60, 0.7); p1 <- sapply(x1, bb_pvalue, n = 60, rho = 0.02)
ok(mean(p1 < 0.05) > 0.5, sprintf("power at p=0.7 = %.2f", mean(p1 < 0.05)))

# 3. rho estimation recovers roughly the truth and handles rho -> 0
xx <- rbinom(2000, 50, rbeta(2000, a, a)); est <- bb_estimate_rho(xx, rep(50, 2000))
ok(est > 0.01 && est < 0.04, sprintf("rho estimate %.4f near 0.02", est))
xb <- rbinom(2000, 50, 0.5); estb <- bb_estimate_rho(xb, rep(50, 2000))
ok(estb < 0.005, sprintf("binomial data gives rho %.5f near 0", estb))
ok(is.finite(bb_pvalue(25, 50, rho = 0)), "rho = 0 falls back to the binomial")

# 4. edge cases: zero total returns NA; x = 0 and x = n give small p at high depth
ok(is.na(bb_pvalue(0, 0, rho = 0.02)), "n = 0 gives NA")
ok(bb_pvalue(0, 80, rho = 0.02) < 1e-6, "all-ref at depth 80 is significant")

# 5. ACAT combination
ok(abs(acat(c(1, 1, 1)) - 1) < 1e-6, "acat of all ones is 1")
ok(acat(c(1e-10, 0.5, 0.5)) < 1e-6, "acat dominated by a tiny p")
ok(acat(c(0.2, 0.3, 0.4)) > 0.1 && acat(c(0.2, 0.3, 0.4)) < 0.5, "acat of moderate p-values is moderate")
quit(status = fail)
```
Run it on the not-yet-written skill: expected FAIL (markers missing).

- [ ] **Step 2: Write the statistics block in the skill**

Add a `## Step 14 — Statistics functions` section with a ```r fence containing exactly:

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
  rho <- plogis(fit$minimum)
  if (rho < 1e-6) 0 else rho
}
bb_pvalue <- function(x, n, rho, p0 = 0.5) {
  if (is.na(n) || n <= 0) return(NA_real_)
  k <- 0:n
  lp <- if (rho < 1e-8) dbinom(k, n, p0, log = TRUE) else {
    a <- p0 * (1 - rho) / rho; b <- (1 - p0) * (1 - rho) / rho
    lchoose(n, k) + lbeta(k + a, n - k + b) - lbeta(a, b)
  }
  min(1, sum(exp(lp[lp <= lp[x + 1] + 1e-9])))
}
acat <- function(p) {   # Cauchy combination with equal weights
  p <- pmin(pmax(p[!is.na(p)], 1e-15), 1); if (length(p) == 0) return(NA_real_)
  0.5 - atan(mean(tan((0.5 - p) * pi))) / pi
}
# --- ase-stats-end
```
Explain in prose that the overdispersion is estimated per sample under H0 (p = 0.5) from all filtered sites, which is conservative when true imbalance exists (inflated ρ lowers power, never inflates false positives).

- [ ] **Step 3: Write Rmd 01 and Rmd 02 templates**

Step 12 (Rmd 01): libraries (`tidyverse`, `openxlsx`, `GenomicRanges`, `rtracklayer`, Bioconductor before tidyverse), constants block (`MIN_DEPTH`, `FDR_SIG`, `ABS_DEV_SIG`), read `samples.csv`, read `ase_counts/{sample}.table` **for every sample listed in `samples.csv`, not by globbing `*.table`**, and stop with a clear error naming each sample whose table is missing, header-only or has no data rows (a stray or stale table from an earlier run must never be accepted; the Task 5 scripts delete failed tables, this is the second line of defence), site filters (total depth at least `MIN_DEPTH`; drop sites where `lowMAPQDepth + lowBaseQDepth + otherBases` exceeds 10 percent of `rawDepth`; drop `improperPairs > 0`), allele mapping (F1: REF count = strain A, ALT count = strain B; outbred: REF/ALT), the reference-bias diagnostic (per-sample mean of REF/(REF+ALT) across sites, a violin/histogram, flagged in a table when the deviation exceeds `max(0.03, 3 × SE)` where SE = `sqrt(0.25 / total reads)`; the fixed 0.02 of the spec is replaced because the verification gate showed a benign +0.02 to −0.03 spread at n of about 2000), for outbred samples the WASP removal table from `ase_counts/{sample}.wasp_stats.tsv` (alignments by `vW` code and by allele class), coverage and site counts, `saveRDS(list(sites, samples, constants), ase_checkpoint.rds)`. Step 13 (Rmd 02): load the checkpoint, paste the statistics block from Step 14, per-sample `rho <- bb_estimate_rho(...)`, per-SNP `bb_pvalue`, BH per sample, `sig` column from `FDR_SIG` and `ABS_DEV_SIG` (used identically by the summary table and every plot), F1 gene-level (sum strain counts per gene from the GTF overlap with `GenomicRanges::findOverlaps`), outbred gene-level (`acat` over SNP p-values, labelled "unphased, no direction"), plots (allelic ratio against depth per sample with significant sites coloured by `sig`), xlsx (`SNP`, `Gene`, `Summary` sheets) and a summary table of counts per sample. Rmds contain no `source()` and use `dplyr::` prefixes.

- [ ] **Step 4: Run the R unit tests (on the cluster) and the checker**

Run `singularity exec --bind /net/bmc-lab3 <bulkrnaseq sif> Rscript ase-pipeline/tests/r/test_ase_stats.R ase-pipeline/ase-pipeline.md` inside an `sbatch -p bcc` job (job output on the shared filesystem); expected: every line `ok` and exit 0. Then `bash ase-pipeline/tests/check_skill.sh ase-pipeline/ase-pipeline.md` (with new `need` lines for `# --- ase-stats-begin`, `bb_pvalue`, `acat`, `ase_checkpoint.rds`, `_ASE_imbalance.xlsx`, `MIN_DEPTH`): expected PASS. If a statistical check fails, fix the function in the skill (never weaken the test).

- [ ] **Step 5: Commit**

```bash
git status && git diff --stat
git add ase-pipeline && git commit -m "ase-pipeline: Rmd 01/02 templates and unit-tested statistics functions"
```

---

### Task 7: Summary report, README, registration

**Files:** Modify `ase-pipeline/ase-pipeline.md`, `ase-pipeline/tests/check_skill.sh`; create `ase-pipeline/README.md`; modify root `README.md`.

- [ ] **Step 1: Failing checker lines** — `need` for `## Step 15`, `_summary_report.html`, `relative links`, `Validation status`, and in the root README a row containing `/ase-pipeline` (checked by a separate `grep -q "/ase-pipeline" README.md`). Run: expected FAIL.
- [ ] **Step 2: Step 15 summary report and Notes** — the standalone `{WD_NAME}_summary_report.html` (inline CSS, relative links to each stage HTML and xlsx, headline numbers per sample, the reference-bias table), the submission order (`sbatch --dependency=afterok:` chains prep → array → Rmd 01 → Rmd 02, each `sbatch -p bcc`), and a "Notes for the assistant" section listing the rules from Global Constraints (containers not modules, no R on the login node, /tmp is node-local, unverified URLs never embedded, duplicates marked not removed, orientation REF=strain A).
- [ ] **Step 3: README** — in the shape of `nfcore-rnavar-setup/README.md`: what it does, prerequisites (containers, `bulkrnaseq` image), usage (`/ase-pipeline`), the wizard steps table, outputs, key design points (third-allele masking, WASP, honest bias diagnostic, RNA-derived genotype caveat), **Validation status** stating exactly what Task 8 exercised and what remains untested (stages 2 and 3, real data), and known limitations. Add one skills-table row to the root `README.md`.
- [ ] **Step 4: Run the checker; expected PASS; commit**

```bash
bash ase-pipeline/tests/check_skill.sh ase-pipeline/ase-pipeline.md
git status && git diff --stat
git add ase-pipeline README.md && git commit -m "ase-pipeline: summary report, README, registration"
```

---

### Task 8: Synthetic acceptance on the cluster (controller-run, not a subagent)

Run by the controller because it needs cluster jobs and judgment; the wizard is followed by hand using the answers below.

- [ ] **Step 1:** Install the branch skill to `~/.claude/commands/ase-pipeline.md` (no file exists there, so nothing is overwritten) and generate the synthetic datasets with Task 1's script into `/net/bmc-lab3/data/bcc/yannvrb/ase_synthetic/`.
- [ ] **Step 2 (F1):** In a fresh directory follow the wizard for `f1/` (custom reference `ase_synthetic/genome/genome.fa` and `genome.gtf`, read length 100, parental VCF from the generator, constants at defaults); submit prep → array → Rmd 01 → Rmd 02 with `--dependency=afterok`.
- [ ] **Step 3 (outbred):** Same for `outbred/` with the individual VCFs; submit the chain.
- [ ] **Step 4:** Evaluate against `truth_genes.tsv`:
  - **F1:** the four planted imbalanced genes are called significant at `FDR_SIG` 0.05 and `ABS_DEV_SIG` 0.1 with the correct **direction** (strain-B fraction above 0.5 for `p_alt` 0.70/0.85, below 0.5 for 0.30); at most one of the null genes is called (nominal rate); on the **null genes only** (truth 0.5) the mean REF fraction after masking is within `max(0.03, 3 × SE)` of 0.5.
  - **Outbred:** same recovery for genes with at least one het SNP; on the null genes the mean REF fraction is within `max(0.03, 3 × SE)` of 0.5 both without and with the WASP filter (the verification gate measured +0.02 unfiltered and −0.03 filtered, so both are reported and no improvement is promised), and the WASP removal table shows the alignments removed by `vW` code and allele class; run the same BAM without the WASP filter and report the difference honestly.
  - Both Rmds build without hand edits; `execution`-style checks: every `ase_counts/*.table` non-empty.
- [ ] **Step 5:** Fix defects through the process (one fix subagent, scoped re-review); record the results and the Validation status wording; do not push without the user's approval.

---

## Self-review notes

- Spec coverage: wizard Steps 0–13 (Tasks 3, 4, 5, 7), F1 and outbred upstream and the mouse helper (Task 5), Rmd 01/02 and statistics (Task 6), verification gate (Task 2), synthetic ground truth and acceptance (Tasks 1, 8), checker and README (Tasks 3, 7). Stages 2 and 3 (Rmd 03/04, phASER) are deliberately outside this plan, as the spec's staging says.
- Consistency: placeholder names are defined once and reused (`{MODE}`, `{SAMPLES_CSV}`, `{PARENTAL_VCF}`, `{STAR_INDEX}`, `{READ_LENGTH}`, `{SJDB_OVERHANG}`); orientation is fixed as REF = strain A / ALT = strain B in Tasks 1, 5, 6 and 8.
- The unit tests test the shipped code (extracted from the skill text between markers), so a passing test cannot drift from what the skill emits.
