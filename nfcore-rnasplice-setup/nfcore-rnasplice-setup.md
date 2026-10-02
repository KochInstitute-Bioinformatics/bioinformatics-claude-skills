# nf-core/rnasplice Pipeline Setup Skill

You are helping the user set up and submit nf-core/rnasplice, a pipeline for differential alternative splicing of bulk RNA-seq, on an HPC cluster using SLURM and Singularity. rnasplice aligns the reads with STAR, quantifies transcripts with Salmon, and runs the analyses the user chooses: rMATS (differential splicing events: skipped exons, alternative 5' and 3' splice sites, mutually exclusive exons, retained introns), SUPPA2 (event and isoform PSI from transcript abundance), DEXSeq and edgeR differential exon usage, and DEXSeq differential transcript usage (DRIMSeq filter, DEXSeq, stageR). Walk through each step below in order, ask only what you need, and perform automated steps silently. This skill configures and submits the run; it does not analyse the results.

---

## Step 0 — Working directory

Before anything else, run `pwd` to record the current working directory (`{CWD}`). The samplesheet, the contrasts sheet, the params file, `nextflow.config`, the scripts and the output directory are written here. Tell the user: "Working directory: `{CWD}`"

Also derive:
- `{WD_NAME}` = basename of `{CWD}`
- `{TODAY_YYMMDD}` = today's date as `YYMMDD`
- `{TODAY_ISO}` = today's date as `YYYY-MM-DD` (used for the output directory `results/{TODAY_ISO}_{WD_NAME}`)

---

## Step 1 — Email address

Ask the user: "What email address should SLURM use for pipeline notifications (completion/failure)?"

Do **not** pre-suggest or pre-fill any email address. Wait for the answer before Step 2. Store `{USER_EMAIL}`.

---

## Step 2 — Conda environment

Ask the user: "Which conda environment should be activated for Nextflow?"

Do **not** pre-suggest any environment name. Ask this as a separate question after Step 1 — never combine Steps 1 and 2 in one message. Store `{CONDA_ENV}`.

Tell the user that this pipeline revision needs Nextflow 26.04.0 or newer (verified with Nextflow 26.04.6).

Do not run `nextflow` here to check the version: Nextflow runs only on compute nodes. The submission script (Step 11) checks the version when the job starts and stops with a clear message if it is outside the verified range.

---

## Step 3 — nf-core/rnasplice version

This skill was verified against nf-core/rnasplice 1b447239488097651d8eac44bca2c1556865eb0f with Nextflow 26.04.6; that is a commit of the development branch (`dev-1b44723`), pinned because release 1.0.4 does not run under that Nextflow (see the README).

Fetch the latest release with the `WebFetch` tool against:
```
https://api.github.com/repos/nf-core/rnasplice/releases/latest
```
Extract `tag_name` (strip a leading `v`) as `{TAG}`.
- Set `{VERSION}` = 1b447239488097651d8eac44bca2c1556865eb0f (`nextflow run -r` accepts a commit).
- If the fetch fails (no network, rate limit, or no `tag_name` in the answer): tell the user that the latest release could not be checked and keep the pinned commit.
- If `{TAG}` is 1.0.4: tell the user that the pinned commit is used (release 1.0.4 does not run under Nextflow 26.04.6).
- If `{TAG}` is newer than 1.0.4: tell the user, fetch `https://raw.githubusercontent.com/nf-core/rnasplice/{TAG}/nextflow_schema.json` with `WebFetch`, list the keys this skill writes (the Step 8 module block and the Step 11 template) that are missing from it, and ask (numbered): 1. Use the pinned commit (verified; default) · 2. Use release {TAG} (not verified by this skill; offer it only when no key is missing). If the user picks the release, `{VERSION}` = {TAG}, and warn the user that the Nextflow requirement (Step 2) and the Nextflow version range the submission script checks (Step 11) were verified for the pinned commit only, not for release {TAG}.

Set `{VERSION_TAG}` = `{VERSION}` when it is a release tag, or `dev-` followed by the first 7 characters when it is a commit; it is used in file names.

**Important:** Do NOT use the `gh` CLI — it is not installed on this HPC cluster. Always use `WebFetch`.

---

## Step 4 — Input files and samplesheet rows

**Scan the working directory first** (directories with their file counts; Nextflow `work/` directories are skipped):
```bash
find "{CWD}" -path "*/work" -prune -o \( -name "*.fastq.gz" -o -name "*.fq.gz" \) -print | sed 's|/[^/]*$||' | sort | uniq -c
find "{CWD}" -path "*/work" -prune -o -name "*.bam" -print | sed 's|/[^/]*$||' | sort | uniq -c
```
Ask (numbered), mentioning what the scan found: "Which input do you want to start from?" 1. FASTQ files (default; the pipeline aligns them with STAR and quantifies them with Salmon) · 2. Genome BAM files from a splice-aware aligner (for example the STAR BAMs of an nf-core/rnaseq run). Store `{SOURCE}` = `fastq` or `genome_bam`.

### Step 4a — FASTQ input

- FASTQ files found: list the unique directories and ask "I found FASTQ files in: `{FOUND_DIRS}`. Use this directory, or specify another?" Nothing found: ask for the full path. Store `{FASTQ_DIR}`.
- **Paired-end detection:** a file whose name has an R1 token — `_R1_`/`_R2_`, `_R1.`/`_R2.`, `_1.fastq.gz`/`_2.fastq.gz`, `_1.fq.gz`/`_2.fq.gz`, or `_1_sequence`/`_2_sequence` — and whose R2 file exists → paired-end; otherwise single-end. Use `fastq_r2_name` (block below) on every file. Store `{LAYOUT}` = `paired` or `single`. rMATS needs one strandedness and one read type for all samples (the pipeline stops with "Cannot run rMats with mixed single and paired end samples"), so if the directory mixes paired-end and single-end files, stop and ask which files to use.
- **Sequencing date:** a leading 6-digit `YYMMDD` filename prefix → `{SEQ_DATE}` (fall back to `{TODAY_YYMMDD}`).
- **Path style:** if `{FASTQ_DIR}` is inside `{CWD}`, use paths relative to `{CWD}`; otherwise absolute paths.
- **Rows:** for each R1 file, `fastq_r2_name` gives its R2 file; if an R2 is missing, warn and stop (one read type for all samples). Single-end rows leave `fastq_2` empty. Sample name before sanitisation = `fastq_sample_name`: the file name up to `_S\d+`, `_R1`/`_R2`, `_1.f`/`_2.f` or `_1_sequence`/`_2_sequence`; a single-end file without these tokens (for example `ctrl.fastq.gz`) only loses its `.fastq.gz`/`.fq.gz` extension. FASTQ files must end in `.fastq.gz` or `.fq.gz`; the pipeline accepts nothing else.
- **Paths:** run `input_path_ok` on every path. A space, tab or comma stops the wizard with a clear message (the samplesheet schemas reject whitespace in file names, and a comma breaks the CSV): the user renames the files or makes symlinks without them.

**Sample names and read pairs.** The pinned samplesheet schemas (FASTQ and genome BAM, recorded in this skill's verification) accept a sample name only if it matches `^(?!\.\.\d+)(?!\.$)[a-zA-Z.]([a-zA-Z0-9._]*)?$` and is not an R reserved word (`if`, `else`, `repeat`, `while`, `function`, `for`, `in`, `next`, `break`, `TRUE`, `FALSE`, `NULL`, `Inf`, `NaN`, `NA`, `NA_integer_`, `NA_real_`, `NA_complex_`, `NA_character_`): "Sample name must be a valid R identifier". Otherwise the pipeline rejects the samplesheet when the job starts. Sanitisation (always, before showing the user) therefore replaces `-`, spaces, `/`, `(`, `)`, `.` and every other character outside `A-Za-z0-9_` with `_`; if a name then does not start with a letter (a digit or `_`), prefix `S`; if it is an R reserved word, append `_S`. Define these functions in the Bash tool and use them for every file (light work on the login node); note every substitution in the preview.
```bash
fastq_r2_name() {
  # usage: fastq_r2_name <R1 file>; prints the R2 file (same directory), or prints nothing and returns 1 when the name has no R1 token
  local dir="" b=$1 r
  case "$1" in */*) dir=${1%/*}/; b=${1##*/} ;; esac
  r=$(printf '%s\n' "$b" | sed -E -n \
    -e 's/^(.*)_R1_/\1_R2_/p;t' \
    -e 's/^(.*)_R1\./\1_R2./p;t' \
    -e 's/^(.*)_1(\.f(ast)?q\.gz)$/\1_2\2/p;t' \
    -e 's/^(.*)_1_sequence/\1_2_sequence/p')
  [ -n "$r" ] || return 1
  printf '%s%s\n' "$dir" "$r"
}
fastq_sample_name() {
  # usage: fastq_sample_name <FASTQ file>; prints the sample name before sanitisation
  printf '%s\n' "${1##*/}" | sed -E 's/(_S[0-9]+[._]|_R[12][._]|_[12]\.f(ast)?q\.gz$|_[12]_sequence|\.f(ast)?q\.gz$).*$//'
}
sanitize_sample_name() {
  # usage: sanitize_sample_name <name>; prints a name that the pinned samplesheet schemas accept
  local s
  s=$(printf '%s' "$1" | LC_ALL=C sed 's/[^A-Za-z0-9_]/_/g')
  case "${s:0:1}" in [ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz]) ;; *) s="S$s" ;; esac
  case "$s" in
    if|else|repeat|while|function|for|in|next|break|TRUE|FALSE|NULL|Inf|NaN|NA|NA_integer_|NA_real_|NA_complex_|NA_character_) s="${s}_S" ;;
  esac
  printf '%s\n' "$s"
}
sample_name_collisions() {
  # usage: one name before sanitisation per row on stdin; prints one line per sanitised name used by more than one row
  local raw
  while IFS= read -r raw; do printf '%s\t%s\n' "$(sanitize_sample_name "$raw")" "$raw"; done | awk -F'\t' '
    !($1 in n) { order[++m] = $1 }
    { n[$1]++ }
    !(($1, $2) in seen) { seen[$1, $2] = 1; k[$1]++; raws[$1] = (k[$1] == 1) ? $2 : raws[$1] ", " $2 }
    END { for (i = 1; i <= m; i++) { s = order[i]; if (n[s] < 2) continue
      if (k[s] == 1) print s ": same name before sanitisation (lanes or technical replicates?)"
      else print s ": created by sanitisation from " raws[s] " (different samples?)" } }'
}
input_path_ok() {
  # usage: input_path_ok <path>; prints STOP and returns 1 when the path has whitespace (the samplesheet schemas reject it) or a comma (it breaks the CSV)
  case "$1" in
    *[[:space:]]*|*,*) echo "STOP: '$1' contains a space, tab or comma, which the samplesheet cannot hold. Rename the file or make a symlink without them, and use that path."; return 1 ;;
  esac
  return 0
}
```
Then run `sample_name_collisions` on the names before sanitisation (one line per row). For every line it prints, warn: "⚠️ Name collision '{NAME}': in nf-core/rnasplice, rows with the same sample name are one sample, and their reads are merged before alignment." Show the original file names of those rows and ask (numbered): 1. the same sample (lanes or technical replicates; the pipeline merges their reads) · 2. different samples — rename them (go to custom naming).
- `same name before sanitisation`: lane files (`X_S1_L001_R1_001`, `X_S1_L002_R1_001`) collide by design; option 1 is expected for them.
- `created by sanitisation` (for example `WT-1` and `WT_1`): the collision comes from the sanitisation, not from lanes. The original names differ, so option 2 is expected unless the user confirms they are one sample.

**Review and naming — order is mandatory:**
1. Show all rows as a table (sample, fastq_1, fastq_2).
2. Ask (numbered): 1. Use auto-generated names · 2. Provide custom names. For custom names: show numbered auto names next to the filenames, ask for a plain-language description, build the mapping, show it as an auto → new table and ask "Does this mapping look correct?" Custom names follow the same rule: `sanitize_sample_name` must return them unchanged (otherwise propose its result); duplicates only where the user chose option 1 above.

The rows stay in memory: Step 5 adds the condition column (the samplesheet columns are `sample,fastq_1,fastq_2,strandedness,condition`) and writes the file.

**Strandedness is asked, never guessed.** rnasplice has no `auto` strandedness, and rMATS needs one strandedness and one read type for all samples. Ask (numbered): "Which strandedness does the library have?" 1. unstranded · 2. forward · 3. reverse · 4. I don't know. Explain: *reverse* = read 1 comes from the strand opposite to the transcript (dUTP libraries such as Illumina TruSeq Stranded and NEBNext Ultra II Directional — the most common stranded libraries); *forward* = read 1 comes from the transcript strand; *unstranded* = no strand information (for example TruSeq non-stranded). A wrong value makes rMATS discard or misassign junction reads without any error. If the user answers 4, or has an nf-core/rnaseq output directory for the same samples, use the block below. Never continue with "I don't know". Store `{STRANDEDNESS}` (one value for every row).

**Strandedness from an existing nf-core/rnaseq run.** Ask for that run's output directory (`{RNASEQ_OUTDIR}`) and count its RSeQC files with `find "{RNASEQ_OUTDIR}" -path "*/work" -prune -o -name "*.infer_experiment.txt" -print | wc -l`. Define the function below in the Bash tool and run it on every one of those files, not a subset (a loop over the same `find`; small text files, light work on the login node). Show a table file → result; with many files, show the count per result and list every file whose result differs from the others. The rule is nf-core/rnaseq's: forward if the forward fraction is at least 0.8, reverse if the reverse fraction is at least 0.8, unstranded if the two fractions differ by at most 0.1, otherwise unclear. A missing or unreadable file, a missing line, or a line that appears twice (reports concatenated into one file) gives unclear. No `*.infer_experiment.txt` file (for example a run with RSeQC skipped): tell the user, and ask the strandedness question again without a proposal; continue only with answer 1, 2 or 3.
```bash
infer_strandedness() {
  # usage: infer_strandedness <file.infer_experiment.txt>; prints forward, reverse, unstranded or unclear
  [ -f "$1" ] && [ -r "$1" ] || { echo "unclear"; return 0; }
  awk -F': ' '
    /explained by "1\+\+,1--,2\+-,2-\+"|explained by "\+\+,--"/ { f = $2 + 0; nf++ }
    /explained by "1\+-,1-\+,2\+\+,2--"|explained by "\+-,-\+"/ { r = $2 + 0; nr++ }
    END {
      if (nf != 1 || nr != 1) { print "unclear"; exit }
      d = f - r; if (d < 0) d = -d
      if (f >= 0.8) print "forward"
      else if (r >= 0.8) print "reverse"
      else if (d <= 0.1 + 1e-9) print "unstranded"
      else print "unclear"
    }' "$1"
}
```
- All files give the same `forward`, `reverse` or `unstranded`: propose it and ask the user to confirm.
- Any file gives `unclear`, or the files disagree: show the fractions (`grep "explained by" FILE`), explain that the library type cannot be read reliably, and ask the user to check the kit; do not continue until they answer 1, 2 or 3.
- The user answered 1-3 and a result disagrees: show both and ask which is right before continuing.

### Step 4b — Genome BAM input

Ask for the directory of the BAM files (for example the `star_salmon/` directory of an nf-core/rnaseq run) and store `{BAM_DIR}`. Genome BAMs must come from a splice-aware aligner (STAR, as in nf-core/rnaseq: `star_salmon/{sample}.markdup.sorted.bam`) and be aligned to the same FASTA and GTF that Step 7 gives the pipeline (same contig names). They need not be sorted or indexed: the pipeline re-sorts and indexes them itself (`samtools sort`, then `samtools index`), so no `.bai` file is needed; a splice-aware aligner is still required. Find them with `find "{BAM_DIR}" -path "*/work" -prune -o -name "*.bam" ! -name "*toTranscriptome*" -print`, derive sample names by removing `.markdup.sorted.bam`, `.umi_dedup.sorted.bam`, `.sorted.bam` or `.bam`, and use one BAM per sample (BAM rows are never merged). When several BAMs give the same name (an nf-core/rnaseq run can keep `X.sorted.bam` next to `X.markdup.sorted.bam`), prefer `.markdup.sorted.bam` (or `.umi_dedup.sorted.bam` in a UMI run), drop the others and tell the user which file was kept. Then apply `input_path_ok`, `sanitize_sample_name` and `sample_name_collisions` (block in Step 4a) as for FASTQ rows. Ask the strandedness question of Step 4a (with the nf-core/rnaseq block when the BAMs come from such a run) and the read type (numbered): 1. paired-end · 2. single-end; store `{STRANDEDNESS}` and `{LAYOUT}`. Single-end and `forward` BAM input were verified from the pipeline code only, not run. Then apply the BAM input rule.

**BAM input rule.** Define this function in the Bash tool and run `bam_input_allowed {STRANDEDNESS} {LAYOUT}`:
```bash
bam_input_allowed() {
  # usage: bam_input_allowed <unstranded|forward|reverse> <paired|single>
  # The four values below are what nf-core/rnasplice dev-1b44723 does with genome-BAM input (recorded by this skill's verification run).
  local RMATS_LIBTYPE="from_sheet" RMATS_READTYPE="from_sheet" DEXSEQ_STRAND="from_sheet" FC_STRAND="from_sheet"
  local lt ds fc
  case "$1" in
    unstranded) lt=fr-unstranded; ds=no; fc=0 ;;
    forward) lt=fr-secondstrand; ds=yes; fc=1 ;;
    reverse) lt=fr-firststrand; ds=reverse; fc=2 ;;
    *) echo "REFUSED: strandedness must be unstranded, forward or reverse (got '$1')"; return 1 ;;
  esac
  case "$2" in paired|single) ;; *) echo "REFUSED: read type must be paired or single (got '$2')"; return 1 ;; esac
  case "$RMATS_LIBTYPE" in
    none) echo "REFUSED: genome-BAM input does not work with this pipeline revision. Start from FASTQ."; return 1 ;;
    from_sheet) echo "ALLOWED: strandedness and read type are written to the BAM samplesheet"; return 0 ;;
  esac
  if [ "$lt" = "$RMATS_LIBTYPE" ] && [ "$2" = "$RMATS_READTYPE" ] \
     && { [ "$DEXSEQ_STRAND" = n/a ] || [ "$ds" = "$DEXSEQ_STRAND" ]; } \
     && { [ "$FC_STRAND" = n/a ] || [ "$fc" = "$FC_STRAND" ]; }; then
    echo "ALLOWED: this revision runs rMATS on BAM input as $RMATS_LIBTYPE, $RMATS_READTYPE-end, which matches the library"; return 0
  fi
  echo "REFUSED: with BAM input this revision runs rMATS as $RMATS_LIBTYPE, $RMATS_READTYPE-end (DEXSeq -s $DEXSEQ_STRAND, featureCounts -s $FC_STRAND), but the library is $1, $2-end. Start from FASTQ instead."
  return 1
}
```
- `ALLOWED`: continue. The samplesheet header is `sample,condition,genome_bam,strandedness,single_end`; write `{STRANDEDNESS}` and `true` (single-end) or `false` in every row (without these two columns the pipeline runs every BAM as unstranded paired-end).
- `REFUSED`: with this revision the rule refuses only a strandedness or read type outside the values above, never a library type (the BAM samplesheet carries both). Show the message and ask the question again with the allowed values; if the user cannot answer, ask (numbered): 1. Start from the FASTQ files instead (go to Step 4a) · 2. Stop here.

With BAM input only rMATS, DEXSeq exon usage and edgeR exon usage can run: DTU and SUPPA2 need Salmon quantification from reads, so Step 8 switches them off. The read length is asked in Step 6.
