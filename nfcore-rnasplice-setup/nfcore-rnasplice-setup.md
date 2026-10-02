# nf-core/rnasplice Pipeline Setup Skill

You are helping the user set up and submit nf-core/rnasplice, a pipeline for differential alternative splicing of bulk RNA-seq, on an HPC cluster using SLURM and Singularity. rnasplice aligns the reads with STAR, quantifies transcripts with Salmon, and runs the analyses the user chooses: rMATS (differential splicing events: skipped exons, alternative 5' and 3' splice sites, mutually exclusive exons, retained introns), SUPPA2 (event and isoform PSI from transcript abundance), DEXSeq and edgeR differential exon usage, and DEXSeq differential transcript usage (DRIMSeq filter, DEXSeq, stageR). Walk through each step below in order, ask only what you need, and perform automated steps silently. This skill configures and submits the run; it does not analyse the results.

---

## Step 0 — Working directory

Before anything else, run `pwd` to record the current working directory (`{CWD}`). The samplesheet, the contrasts sheet, the params file, `nextflow.config`, the scripts and the output directory are written here. Tell the user: "Working directory: `{CWD}`"

Also derive:
- `{WD_NAME}` = basename of `{CWD}`
- `{TODAY_YYMMDD}` = today's date as `YYMMDD`
- `{TODAY_ISO}` = today's date as `YYYY-MM-DD` (used for the output directory `results/{TODAY_ISO}_{WD_NAME}`)

**Bash tool calls.** Shell variables and functions do not persist between Bash tool calls. When a step defines functions in a code block, paste the whole block at the start of every Bash call that uses its functions, in that same call. A value needed in a later call (a path, a result) is written out literally there (or kept as a placeholder such as `{CWD}`), never taken from a shell variable of an earlier call.

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

**Sample names and read pairs.** The pinned samplesheet schemas (FASTQ and genome BAM, recorded in this skill's verification) accept a sample name only if it matches `^(?!\.\.\d+)(?!\.$)[a-zA-Z.]([a-zA-Z0-9._]*)?$` and is not an R reserved word (`if`, `else`, `repeat`, `while`, `function`, `for`, `in`, `next`, `break`, `TRUE`, `FALSE`, `NULL`, `Inf`, `NaN`, `NA`, `NA_integer_`, `NA_real_`, `NA_complex_`, `NA_character_`): "Sample name must be a valid R identifier". Otherwise the pipeline rejects the samplesheet when the job starts. Sanitisation (always, before showing the user) therefore replaces `-`, spaces, `/`, `(`, `)`, `.` and every other character outside `A-Za-z0-9_` with `_`; if a name then does not start with a letter (a digit or `_`), prefix `S`; if it is an R reserved word, append `_S`. Paste this block at the start of every Bash call that uses these functions and use them for every file (light work on the login node); note every substitution in the preview.
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

**Strandedness from an existing nf-core/rnaseq run.** Ask for that run's output directory (`{RNASEQ_OUTDIR}`) and count its RSeQC files with `find "{RNASEQ_OUTDIR}" -path "*/work" -prune -o -name "*.infer_experiment.txt" -print | wc -l`. In one Bash call, define the function below and run it on every one of those files, not a subset (a loop over the same `find`; small text files, light work on the login node). Show a table file → result; with many files, show the count per result and list every file whose result differs from the others. The rule is nf-core/rnaseq's: forward if the forward fraction is at least 0.8, reverse if the reverse fraction is at least 0.8, unstranded if the two fractions differ by at most 0.1, otherwise unclear. A missing or unreadable file, a missing line, or a line that appears twice (reports concatenated into one file) gives unclear. No `*.infer_experiment.txt` file (for example a run with RSeQC skipped): tell the user, and ask the strandedness question again without a proposal; continue only with answer 1, 2 or 3.
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

Ask for the directory of the BAM files (for example the `star_salmon/` directory of an nf-core/rnaseq run) and store `{BAM_DIR}`. For BAM input, `{SEQ_DATE}` = a leading 6-digit `YYMMDD` prefix shared by the BAM file names, else `{TODAY_YYMMDD}` (used for the file names in Step 5). Genome BAMs must come from a splice-aware aligner (STAR, as in nf-core/rnaseq: `star_salmon/{sample}.markdup.sorted.bam`) and be aligned to the same FASTA and GTF that Step 7 gives the pipeline (same contig names). They need not be sorted or indexed: the pipeline re-sorts and indexes them itself (`samtools sort`, then `samtools index`), so no `.bai` file is needed; a splice-aware aligner is still required. Find them with `find "{BAM_DIR}" -path "*/work" -prune -o -name "*.bam" ! -name "*toTranscriptome*" -print`, derive sample names by removing `.markdup.sorted.bam`, `.umi_dedup.sorted.bam`, `.sorted.bam` or `.bam`, and use one BAM per sample (BAM rows are never merged). When several BAMs give the same name (an nf-core/rnaseq run can keep `X.sorted.bam` next to `X.markdup.sorted.bam`), prefer `.markdup.sorted.bam` (or `.umi_dedup.sorted.bam` in a UMI run), drop the others and tell the user which file was kept. Then apply `input_path_ok`, `sanitize_sample_name` and `sample_name_collisions` (block in Step 4a) as for FASTQ rows. Ask the strandedness question of Step 4a (with the nf-core/rnaseq block when the BAMs come from such a run) and the read type (numbered): 1. paired-end · 2. single-end; store `{STRANDEDNESS}` and `{LAYOUT}`. Single-end and `forward` BAM input were verified from the pipeline code only, not run. Then apply the BAM input rule.

**BAM input rule.** In one Bash call, define this function and run `bam_input_allowed {STRANDEDNESS} {LAYOUT}`:
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

---

## Step 5 — Conditions, contrasts and the two sheets

**Conditions.** Ask the user to describe the groups in plain language (for example "samples 1-3 are wild type, 4-6 are Rbpms2 knockout"), interpret it, and show a table sample → condition; ask "Is this correct?" until confirmed. Condition labels use letters, digits and `_`, start with a letter and are not R reserved words such as `NA`, `TRUE`, `FALSE`, `in` (the pipeline's R scripts read the sheet with `read.csv`, and edgeR builds its contrasts from the labels; suggest short labels such as `WT`, `KO`); rows of one sample (technical replicates) share its condition. Every condition used in a contrast needs at least 2 samples (distinct sample names), because rMATS, DEXSeq and edgeR estimate variability from replicates; if one has fewer, say so and stop until the design is changed. With only one condition there is nothing to compare: stop.

Then ask for the reference: with two conditions, "Which condition is the reference (control)?" (numbered); with more, ask for an order with the reference first. Store `{CONDITIONS}` = the labels in that order.

**Contrasts.** Ask (numbered): 1. All pairwise comparisons · 2. Only the comparisons I list.
- Option 1: for every pair of conditions, treatment = the later and control = the earlier one in `{CONDITIONS}`; contrast name `{treatment}_vs_{control}`. With two conditions this is one contrast.
- Option 2: ask for lines "treatment vs control", one per contrast; names as above.

Show the contrasts as a table (contrast, treatment, control). rMATS runs one prep/post pair per contrast; with more than 6 contrasts tell the user that the run takes correspondingly longer. DTU and SUPPA2 name their outputs `<treatment>-<control>` (not the contrast name), so the same treatment and control may appear in one contrast only.
When showing the table, tell the user the direction: treatment = the later condition (in option 2, the condition named before "vs"), and a positive rMATS IncLevelDifference means more inclusion in the treatment.

**Paired design.** Ask only when there are exactly two conditions with the same number of samples; otherwise set `{PAIRED_DESIGN}` = `false` without asking. Ask (numbered): "Are the samples paired (each sample of one condition has a partner from the same individual, litter or batch in the other)?" 1. No, the samples are independent (default) · 2. Yes, every sample has a partner. Explain: the pipeline default `diffsplice_paired: true` assumes a pairing (SUPPA2; `rmats_paired_stats` defaults to `false` at this revision); a paired test on independent samples gives wrong statistics, so this skill writes both `false` unless the user confirms pairing. Answer 1: `{PAIRED_DESIGN}` = `false`.
Answer 2: rMATS pairs the i-th sample of one condition with the i-th sample of the other, in samplesheet order (verified for this revision). Ask for the partner of each sample as a subject label (letters, digits and `_`, for example `mouse1`), check that each subject has exactly one sample in each condition, and write the rows ordered by subject within each condition: sort the rows by condition (in `{CONDITIONS}` order), then by subject, so that the i-th sample of each condition belongs to the same subject; the rows of one sample (technical replicates) stay next to each other. Show the pairs as a table (subject, sample of each condition) and ask "Is this pairing correct?" until confirmed. The labels are written to a pairs file (header `sample,subject`, one row per sample) that the validator below checks position by position. `{PAIRED_DESIGN}` = `true`.
`{PAIRED_DESIGN}` sets both `rmats_paired_stats` and SUPPA2's `diffsplice_paired` in Step 8.

**File names** (numbered; for BAM input `{SEQ_DATE}` comes from Step 4b): 1. `{SEQ_DATE}_{WD_NAME}` · 2. `{TODAY_YYMMDD}_{WD_NAME}` · 3. Custom prefix. Store `{SHEET_PREFIX}`; `{SAMPLESHEET_CSV}` = `{SHEET_PREFIX}_samplesheet.csv` and `{CONTRASTS_CSV}` = `{SHEET_PREFIX}_contrasts.csv` (with a paired design also `{SHEET_PREFIX}_pairs.csv`, a record of the pairing that the pipeline does not read).

**Sheet validation.** Both sheets are written to a scratch directory first, validated there, and only then moved into `{CWD}`. The block below defines `validate_rnasplice_sheets` and `install_rnasplice_sheets` (awk and coreutils only; light work on the login node); it is pasted into the one Bash call of **Writing the sheets** below.
```bash
validate_rnasplice_sheets() {
  # usage: validate_rnasplice_sheets <samplesheet.csv> <contrasts.csv> <fastq|genome_bam> <paired design: 0|1> [pairs.csv]
  # pairs.csv (header sample,subject) is required for a paired design: the i-th sample of each condition must share a subject
  # prints "SHEETS OK", or one "ERROR: ..." line per problem and returns 1
  [ -s "$1" ] && [ -s "$2" ] || { echo "ERROR: samplesheet or contrasts file missing or empty"; return 1; }
  case "$3" in fastq|genome_bam) ;; *) echo "ERROR: source must be fastq or genome_bam (got '$3')"; return 1 ;; esac
  case "$4" in 0|1) ;; *) echo "ERROR: paired design must be 0 or 1 (got '$4': 1 = paired design, 0 = not paired)"; return 1 ;; esac
  [ -z "${5:-}" ] || [ -s "$5" ] || { echo "ERROR: pairs file $5 missing or empty"; return 1; }
  awk -F',' -v src="$3" -v paired="$4" \
      -v rw="^(if|else|repeat|while|function|for|in|next|break|TRUE|FALSE|NULL|Inf|NaN|NA|NA_integer_|NA_real_|NA_complex_|NA_character_)$" \
      -v fqhdr="sample,fastq_1,fastq_2,strandedness,condition" -v bamhdr="sample,condition,genome_bam,strandedness,single_end" '
    function err(m) { print "ERROR: " m; bad = 1 }
    FNR == 1 { file++ }
    { sub(/\r$/, "") }
    /"/ { err("file " file " line " FNR " contains a double quote; write plain comma-separated values"); next }
    file == 1 && FNR == 1 {
      want = (src == "fastq") ? fqhdr : bamhdr
      if ($0 != want) { err("samplesheet header is \"" $0 "\", expected \"" want "\""); hdrbad = 1 }
      ncol = split(want, h, ","); for (i = 1; i <= ncol; i++) col[h[i]] = i
      next
    }
    file == 1 && $0 == "" { next }
    file == 1 {
      if (NF != ncol) { err("samplesheet line " FNR " has " NF " fields, expected " ncol); next }
      s = $col["sample"]; c = $col["condition"]
      if (s !~ /^[A-Za-z][A-Za-z0-9_]*$/) err("sample name \"" s "\" (line " FNR "): letters, digits and _ only, starting with a letter")
      else if (s ~ rw) err("sample name \"" s "\" is an R reserved word; the pipeline rejects it")
      if (c !~ /^[A-Za-z][A-Za-z0-9_]*$/) err("condition \"" c "\" of " s ": letters, digits and _ only, starting with a letter")
      else if (c ~ rw) err("condition \"" c "\" of " s " is an R reserved word; the R scripts of the pipeline would read it as a value")
      if (s in cond) {
        if (src != "fastq") err("sample " s " appears twice; a BAM samplesheet has one row per sample")
        else if (cond[s] != c) err("rows of sample " s " have different conditions (" cond[s] ", " c "); rows with one sample name are merged into one sample")
        else if (paired == 1 && s != last) err("rows of sample " s " are not next to each other; in a paired design keep the rows of one sample together")
      } else { nsamp[c]++; ord[c, nsamp[c]] = s }
      cond[s] = c; conds[c] = 1; last = s
      if (src == "fastq") {
        if ($col["fastq_1"] !~ /\.f(ast)?q\.gz$/) err("fastq_1 of " s " does not end in .fastq.gz or .fq.gz")
        if ($col["fastq_2"] != "" && $col["fastq_2"] !~ /\.f(ast)?q\.gz$/) err("fastq_2 of " s " does not end in .fastq.gz or .fq.gz")
        lay = ($col["fastq_2"] == "") ? "single-end" : "paired-end"
      } else {
        if ($col["genome_bam"] !~ /\.bam$/) err("genome_bam of " s " does not end in .bam")
        lay = "n/a"
        if ("single_end" in col) {
          if ($col["single_end"] !~ /^(true|false)$/) err("single_end of " s " must be true or false")
          lay = ($col["single_end"] == "true") ? "single-end" : "paired-end"
        }
      }
      if (lay != "n/a") { if (layout == "") layout = lay; else if (lay != layout) err("mixed single-end and paired-end samples (" s "); rMATS needs one read type for all samples") }
      if ("strandedness" in col) {
        st = $col["strandedness"]
        if (st !~ /^(unstranded|forward|reverse)$/) err("strandedness \"" st "\" of " s ": use unstranded, forward or reverse")
        else if (strand == "") strand = st; else if (st != strand) err("mixed strandedness (" strand ", " st "); rMATS needs one strandedness for all samples")
      }
      next
    }
    file == 2 && FNR == 1 { if ($0 != "contrast,treatment,control") err("contrasts header is \"" $0 "\", expected \"contrast,treatment,control\""); next }
    file == 2 && $0 == "" { next }
    file == 2 {
      if (NF != 3) { err("contrasts line " FNR " has " NF " fields, expected 3"); next }
      ncon++
      if ($1 !~ /^[A-Za-z][A-Za-z0-9_-]*$/) err("contrast name \"" $1 "\": letters, digits, _ and - only, starting with a letter")
      if ($1 in seen) err("contrast name " $1 " is used twice"); seen[$1] = 1
      if (!($2 in conds)) err("treatment \"" $2 "\" of contrast " $1 " is not a value of the condition column")
      if (!($3 in conds)) err("control \"" $3 "\" of contrast " $1 " is not a value of the condition column")
      if ($2 == $3) err("contrast " $1 " compares " $2 " with itself")
      if (($2, $3) in tc) err("contrasts " tc[$2, $3] " and " $1 " both compare " $2 " with " $3 "; DTU and SUPPA2 name their outputs <treatment>-<control>, so one would overwrite the other")
      else tc[$2, $3] = $1
      used[$2] = 1; used[$3] = 1
      next
    }
    file == 3 && paired != 1 { next }
    file == 3 && FNR == 1 { if ($0 != "sample,subject") err("pairs header is \"" $0 "\", expected \"sample,subject\""); next }
    file == 3 && $0 == "" { next }
    file == 3 {
      if (NF != 2) { err("pairs line " FNR " has " NF " fields, expected 2"); next }
      if (!($1 in cond)) { err("pairs file names sample " $1 ", which is not in the samplesheet"); next }
      if ($1 in subj) { err("sample " $1 " appears twice in the pairs file"); next }
      if ($2 !~ /^[A-Za-z0-9][A-Za-z0-9_]*$/) err("subject \"" $2 "\" of " $1 ": letters, digits and _ only")
      subj[$1] = $2
      if (++per[cond[$1], $2] == 2) err("subject " $2 " has 2 samples in condition " cond[$1] "; each subject needs exactly one sample per condition")
      next
    }
    END {
      if (file < 2) err("contrasts file was not read")
      if (ncon == 0) err("contrasts file has no contrast rows")
      if (!hdrbad) for (c in used) if ((c in nsamp) && nsamp[c] < 2) err("condition " c " has " nsamp[c] " sample; each compared condition needs at least 2")
      if (paired == 1) {
        k = 0; for (c in conds) { k++; cn[k] = c; sz[k] = nsamp[c] }
        if (k != 2) err("a paired design needs exactly two conditions, found " k)
        else if (sz[1] != sz[2]) err("a paired design needs the same number of samples in both conditions (" sz[1] ", " sz[2] ")")
        else if (file < 3) err("a paired design needs the pairs file (header sample,subject) as fifth argument")
        else {
          for (s in cond) if (!(s in subj)) err("sample " s " has no subject in the pairs file")
          for (s in subj) { has[cond[s], subj[s]] = 1; subjects[subj[s]] = 1 }
          for (u in subjects) for (j = 1; j <= 2; j++) if (!((cn[j], u) in has)) { err("subject " u " has no sample in condition " cn[j]); lone = 1 }
          if (!lone) for (i = 1; i <= sz[1]; i++) {
            s1 = ord[cn[1], i]; s2 = ord[cn[2], i]
            if ((s1 in subj) && (s2 in subj) && subj[s1] != subj[s2]) err("paired design: sample " i " of condition " cn[1] " is " s1 " (subject " subj[s1] ") but sample " i " of condition " cn[2] " is " s2 " (subject " subj[s2] "); sort the rows of each condition by subject")
          }
        }
      }
      if (!bad) print "SHEETS OK"
      exit bad
    }' "$1" "$2" ${5:+"$5"}
}
install_rnasplice_sheets() {
  # usage: install_rnasplice_sheets <scratch dir> <fastq|genome_bam> <paired design: 0|1> <samplesheet name> <contrasts name> [pairs name]
  # Run in the directory {CWD}, in the same Bash call that wrote <scratch dir>/samplesheet.csv, contrasts.csv (and pairs.csv).
  # Validates the sheets, checks that every input file exists (paths relative to the current directory; URLs are not checked),
  # then moves the sheets into the current directory under the given names.
  # The scratch files and directory are removed on every path (success or error), so nothing is left behind
  # when the user stops after an error.
  local T=$1 src=$2 paired=$3 rc=0 p
  [ -d "$T" ] || { echo "ERROR: scratch directory '$T' not found"; return 1; }
  case "$paired" in 0|1) ;; *) echo "ERROR: paired design must be 0 or 1 (got '$paired')"; rc=1 ;; esac
  for p in "$4" "$5" ${6:+"$6"}; do case "$p" in ''|*/*) echo "ERROR: target name '$p' must be a plain file name in the current directory"; rc=1 ;; esac; done
  [ "$paired" != 1 ] || [ -n "${6:-}" ] || { echo "ERROR: a paired design needs the pairs file name"; rc=1; }
  if [ $rc -eq 0 ]; then
    if [ "$paired" = 1 ]; then
      validate_rnasplice_sheets "$T/samplesheet.csv" "$T/contrasts.csv" "$src" 1 "$T/pairs.csv" || rc=1
    else
      validate_rnasplice_sheets "$T/samplesheet.csv" "$T/contrasts.csv" "$src" 0 || rc=1
    fi
  fi
  if [ $rc -eq 0 ]; then
    while IFS= read -r p; do
      [ -s "$p" ] || { echo "ERROR: input file $p is missing or empty: fix the path with the user (Step 4), then run the whole call again"; rc=1; }
    done < <(awk -F',' '{ sub(/\r$/, "") } NR == 1 { for (i = 1; i <= NF; i++) if ($i ~ /^(fastq_1|fastq_2|genome_bam)$/) k[i] = 1; next }
                        { for (i in k) if ($i != "" && $i !~ /:\/\//) print $i }' "$T/samplesheet.csv")
  fi
  if [ $rc -eq 0 ]; then
    mv -f "$T/samplesheet.csv" "./$4" && mv -f "$T/contrasts.csv" "./$5" || rc=1
    [ "$paired" != 1 ] || mv -f "$T/pairs.csv" "./$6" || rc=1
  fi
  rm -f "$T/samplesheet.csv" "$T/contrasts.csv" "$T/pairs.csv"
  rmdir "$T" || echo "WARNING: scratch directory $T is not empty: remove its files by name, then rmdir it"
  [ $rc -eq 0 ] && echo "WROTE $4 $5${6:+ $6}"
  return $rc
}
```
**Writing the sheets.** Shell variables and functions do not persist between Bash tool calls, so the whole procedure is one Bash call:
1. Before that call, check which of `{SAMPLESHEET_CSV}`, `{CONTRASTS_CSV}` (and `{SHEET_PREFIX}_pairs.csv` with a paired design) already exist in `{CWD}` (`ls`). For each, ask (numbered): 1. overwrite · 2. choose another filename (store the new name in `{SAMPLESHEET_CSV}` or `{CONTRASTS_CSV}`, which later steps read; for the pairs file, replace `{SHEET_PREFIX}_pairs.csv` in the call).
2. In one Bash call: paste the whole sheet-validation block above (both functions), then the block below with the rows filled in (samplesheet header `sample,fastq_1,fastq_2,strandedness,condition` for FASTQ with `{STRANDEDNESS}` in every row, or `sample,condition,genome_bam,strandedness,single_end` for BAM; the sanitised names; with a paired design the rows sorted as described under **Paired design**). Use the last line for an independent design, or replace it by the commented paired line (with the pairs file) when `{PAIRED_DESIGN}` is `true`.
```bash
cd "{CWD}" || exit 1
T=$(mktemp -d) || exit 1
cat > "$T/samplesheet.csv" <<'END_OF_SHEET'
{SAMPLESHEET HEADER AND ROWS}
END_OF_SHEET
cat > "$T/contrasts.csv" <<'END_OF_SHEET'
contrast,treatment,control
{CONTRAST ROWS}
END_OF_SHEET
# paired design only:
# cat > "$T/pairs.csv" <<'END_OF_SHEET'
# sample,subject
# {SAMPLE,SUBJECT ROWS}
# END_OF_SHEET
# install_rnasplice_sheets "$T" {SOURCE} 1 {SAMPLESHEET_CSV} {CONTRASTS_CSV} {SHEET_PREFIX}_pairs.csv
install_rnasplice_sheets "$T" {SOURCE} 0 {SAMPLESHEET_CSV} {CONTRASTS_CSV}
```
3. On `ERROR` lines nothing was moved and the scratch directory is already removed. Fix the cause with the user — never edit the validator — and run the whole call again. A missing or empty input file stops the call the same way: correct the path in the rows (Step 4) or stop.
4. On `WROTE`, show the files in full.

Keep for Step 8: the number of distinct samples, the size of the smallest compared condition and, per contrast, the number of treatment and control samples.

---

## Step 6 — Read length for rMATS

rMATS needs the read length (`rmats_read_len`), and the pipeline default of 40 is wrong for almost all data, so the length is always set here and written to the params file. rMATS uses this value to normalise the inclusion and skipping counts of every sample, and always runs with `--variable-read-length` and `--allow-clipping`, so trimmed reads of other lengths are still counted. No other read-length setting is used by the analyses this skill runs (the STAR index the pipeline builds does not depend on it; `miso_read_len` belongs to the MISO sashimi plots, which this skill leaves off).

**Read-length detection.** FASTQ input. Every row of the samplesheet is checked: its R1 file and, for paired-end data, its R2 file, each with its first 1000 reads only (`head -n 4000`), which is light work on the login node. In one Bash call, paste this block, then the lines under **Running the read-length check** below.
```bash
read_lengths_of() {
  # usage: read_lengths_of <fastq.gz>; prints the length of each of the first 1000 reads (a CR at the line end is ignored),
  # or one "ERROR: ..." line and returns 1 when the file is missing, empty, unreadable, truncated, corrupt or not FASTQ
  local f=$1 e out msg
  [ -f "$f" ] && [ -r "$f" ] && [ -s "$f" ] || { echo "ERROR: $f is missing, empty or unreadable"; return 1; }
  e=$(mktemp) || { echo "ERROR: mktemp failed"; return 1; }
  out=$(zcat -f -- "$f" 2> "$e" | head -n 4000 | awk '
    { sub(/\r$/, "") }
    NR % 4 == 1 && !/^@/ { bad = 1 }
    NR % 4 == 2 { s = length($0) }
    NR % 4 == 3 && !/^\+/ { bad = 1 }
    NR % 4 == 0 { if (length($0) != s) bad = 1; print s }
    END { if (bad || NR == 0 || NR % 4 != 0) print "BAD" }')
  msg=$(grep -v 'Broken pipe' "$e" | head -n 1)
  rm -f "$e"
  [ -z "$msg" ] || { echo "ERROR: $f cannot be read: $msg"; return 1; }
  case "$out" in *BAD*) echo "ERROR: $f is not a complete FASTQ file (truncated, corrupt or not FASTQ)"; return 1 ;; esac
  printf '%s\n' "$out"
}
detect_read_length() {
  # usage: detect_read_length <fastq.gz> ...
  # prints "<length> <reads>": the most common read length among the first 1000 reads of each file (a tie goes to the longer
  # length), or the "ERROR: ..." line of the first file that cannot be read, and returns 1
  local f l all=""
  [ $# -gt 0 ] || { echo "ERROR: no FASTQ file given"; return 1; }
  for f in "$@"; do l=$(read_lengths_of "$f") || { echo "$l"; return 1; }; all+="$l"$'\n'; done
  printf '%s' "$all" | sort -n | uniq -c | sort -k1,1nr -k2,2nr | head -1 | awk '{print $2, $1}'
}
read_length_report() {
  # usage: read_length_report < rows (one per samplesheet row: sample<TAB>condition<TAB>fastq_1<TAB>fastq_2, fastq_2 empty if single-end)
  # prints "SAMPLE <sample> <condition> R1 <length> [R2 <length>]" per row, then "READ_LENGTH <n>" (the most common length over
  # the reads of every R1 file) and one verdict: "OK: ..." (one length everywhere), "WARNING: ..." (lengths differ, but every
  # condition has every length) or "CONFOUNDED: ..." (lengths differ between conditions).
  # A file that cannot be read gives its ERROR line, a final "ERROR: read length not set ..." line and return 1 (no READ_LENGTH).
  local s c r1 r2 l1 l2 rc=0 rows="" r1s=()
  while IFS=$'\t' read -r s c r1 r2; do
    [ -n "$s" ] || continue
    l1=$(detect_read_length "$r1") || { echo "$l1"; rc=1; continue; }
    l2=""
    if [ -n "$r2" ]; then l2=$(detect_read_length "$r2") || { echo "$l2"; rc=1; continue; }; fi
    echo "SAMPLE $s $c R1 ${l1%% *}${l2:+ R2 ${l2%% *}}"
    rows+="$c"$'\t'"${l1%% *}"$'\t'"${l2%% *}"$'\n'; r1s+=("$r1")
  done
  [ $rc -eq 0 ] || { echo "ERROR: read length not set: fix or replace the files above (Step 4), then run the whole call again"; return 1; }
  [ ${#r1s[@]} -gt 0 ] || { echo "ERROR: read length not set: no samplesheet rows were given"; return 1; }
  l1=$(detect_read_length "${r1s[@]}") || { echo "$l1"; return 1; }
  echo "READ_LENGTH ${l1%% *}"
  printf '%s' "$rows" | awk -F'\t' '
    function add(c, l) {
      if (!((c, l) in has)) { has[c, l] = 1; per[c] = per[c] (per[c] == "" ? "" : "/") l }
      if (!(l in len)) { len[l] = 1; nl++; all = all (all == "" ? "" : ", ") l }
    }
    !($1 in seen) { seen[$1] = 1; order[++nc] = $1 }
    { add($1, $2); if ($3 != "") add($1, $3) }
    END {
      if (nl == 1) { print "OK: every file has " all " bp reads"; exit }
      for (i = 1; i <= nc; i++) for (l in len) if (!((order[i], l) in has)) conf = 1
      for (i = 1; i <= nc; i++) sum = sum (i > 1 ? "; " : "") order[i] ": " per[order[i]] " bp"
      if (conf) print "CONFOUNDED: read lengths differ between conditions (" sum ")"
      else print "WARNING: read lengths differ (" all " bp), but every condition has every length (" sum ")"
    }'
}
```
**Running the read-length check.** After the block above, in the same Bash call (the rows come from the samplesheet written in Step 5; paths are relative to `{CWD}`):
```bash
cd "{CWD}" || exit 1
awk -F',' '{ sub(/\r$/, "") } NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next } $0 != "" {
  print $c["sample"] "\t" $c["condition"] "\t" $c["fastq_1"] "\t" $c["fastq_2"] }' "{SAMPLESHEET_CSV}" | read_length_report
```
Show the `SAMPLE` lines as a table grouped by condition (condition, sample, R1 length, R2 length). `{READ_LENGTH}` = the number on the `READ_LENGTH` line (the most common R1 length over all samples), never a single sample's value. A `gzip: stdout: Broken pipe` message is harmless (`head` stops reading early) and is already filtered out.
- `ERROR` lines: there is no `READ_LENGTH` line, and `{READ_LENGTH}` is not set. Show the errors and ask the user to fix or replace the files (back to Step 4) or stop; never continue with an empty or guessed read length.
- `OK`: tell the user.
- `WARNING`: warn loudly: "⚠️ Read lengths differ between files ({the lengths}). rMATS normalises every sample with one read length ({READ_LENGTH} bp); the lengths are spread over all conditions, so the comparison is not biased by them, but tell me if this is unexpected." Continue.
- `CONFOUNDED`: warn loudly: "⚠️ Read lengths differ between conditions ({the per-condition lengths}). rMATS normalises the inclusion and skipping counts of every sample with one read length ({READ_LENGTH} bp), so the samples of the condition with other lengths are normalised with the wrong length: their PSI values, and the differences between the conditions, are biased, and rMATS gives no error." Ask (numbered): 1. Stop here (default) — trim all reads to one common length (for example the shortest) outside this wizard, then run it again · 2. Continue with {READ_LENGTH} bp anyway (the bias stays).

Tell the user: "Read length for rMATS: {READ_LENGTH} bp (`rmats_read_len: {READ_LENGTH}`)."

BAM input: ask "What is the read length of the sequencing (for example 100 or 150)?" — an integer between 20 and 1000; store it as `{READ_LENGTH}`.

---

## Step 7 — Organism and genome files

Ask, in this order:
1. "What organism is this data from? (e.g. mouse, human)"
2. (numbered): 1. Ensembl release in the standard folder (default) · 2. Custom reference — I already have a FASTA and GTF.
3. The option-specific questions below. Ask the base directory only for option 1.

- **Option 1 (Ensembl):** ask "What is the base directory where genome files and indexes are stored?" (`{genome_base}`). Folder convention, shared with `/nfcore-rnaseq-setup` so that FASTA, GTF and the STAR index are reused, never downloaded or built twice:
  ```
  {genome_base}/{organism}/{assembly}_ens{version}/
  ├── {FASTA}.fa            ← primary assembly FASTA
  ├── {GTF}.gtf             ← annotation GTF
  └── index/
      ├── star/             ← STAR index (as built by /nfcore-rnaseq-setup); reused only if compatible
      └── salmon/           ← not used: rnasplice always builds its own Salmon index
  ```
  Mouse: GRCm39, FASTA `Mus_musculus.GRCm39.dna.primary_assembly.fa`, GTF `Mus_musculus.GRCm39.{version}.gtf`, directory `{genome_base}/mouse/mm39_ens{version}/`. Human: GRCh38, FASTA `Homo_sapiens.GRCh38.dna.primary_assembly.fa`, GTF `Homo_sapiens.GRCh38.{version}.gtf`, directory `{genome_base}/human/hg38_ens{version}/`. Other organisms: option 2. Version: the highest existing `{assembly}_ens{N}` directory unless the user asks otherwise; if none exists, the latest Ensembl release from `https://ftp.ensembl.org/pub/current/README` (`WebFetch`; its line "Ensembl Release N Databases." gives N). If that fetch fails, use the highest `release-N/` directory in the listing of `https://ftp.ensembl.org/pub/` (`WebFetch`); if both fail, ask the user for the release. Store `{GENOME_DIR}`, `{FASTA_PATH}`, `{GTF_PATH}`, `{ORGANISM}`, `{ASSEMBLY}`, `{ENS_VERSION}` and `{REF_TAG}` = `{ASSEMBLY}_ens{ENS_VERSION}`.
  If the FASTA or the GTF is missing: resolve both download URLs at run time with `WebFetch` on the Ensembl FTP listing of that release (`https://ftp.ensembl.org/pub/release-{version}/fasta/{species}/dna/` and `https://ftp.ensembl.org/pub/release-{version}/gtf/{species}/`), check each with a HEAD request (`curl -sI URL | head -1` must show 200), show them to the user, store `{ENSEMBL_FASTA_URL}` and `{ENSEMBL_GTF_URL}`, and generate `download_genome_{REF_TAG}.sh` (Step 11). Never type a URL from memory.
- **Option 2 (Custom reference):** ask for the FASTA path and the GTF path; `{GENOME_DIR}` = the directory of the FASTA, `{REF_TAG}` = `custom_{WD_NAME}`; no download. The pipeline accepts `.fa`, `.fasta`, `.fna` and the same with `.gz`, and `.gtf` or `.gtf.gz`. Check that both files exist (`test -s`). Ask (numbered) whether a STAR index built from exactly this FASTA and GTF exists: 1. No (default) · 2. Yes — ask for its directory and apply the compatibility checks below.

**GTF source.** `zcat -f {GTF_PATH} | grep -v "^#" | head -3`: gene IDs with a version suffix (`ENSG00000000003.15`) indicate GENCODE, without one Ensembl. Report it. `gencode: false` is written either way: the pipeline builds the transcript FASTA from the GTF itself, so it never receives a GENCODE-format transcript FASTA, the only case `gencode: true` is for.

**Existing STAR index — reused only when compatible; otherwise the pipeline builds its own.** (BAM input uses no index: skip this part; `{STAR_INDEX}` is empty.) The wizard never builds an index itself.
- STAR (`{STAR_DIR}` = `{GENOME_DIR}/index/star` for option 1, or the directory the user gave): in one Bash call, check the files and read the two header lines (small text files, light work on the login node):
  ```bash
  for f in SA SAindex Genome sjdbList.out.tab genomeParameters.txt; do [ -s "{STAR_DIR}/$f" ] || echo "MISSING: $f"; done
  grep -E '^(versionGenome|sjdbOverhang)[[:space:]]' "{STAR_DIR}/genomeParameters.txt"
  ```
  Any `MISSING` line: there is no usable index; let the pipeline build it. Otherwise offer reuse — (numbered) 1. Reuse (default) · 2. Let the pipeline build its own — only when `versionGenome` is `2.7.4a` (the index format of the STAR inside this pipeline revision); with any other value, say why and let the pipeline build it. If its `sjdbOverhang` is not 100 (the value the pipeline uses when it builds the index itself), add a note (information only; it is not a reason to rebuild, and indexes from /nfcore-rnaseq-setup usually have read length − 1): "This STAR index was built with sjdbOverhang {N}, not the pipeline's 100. STAR aligns with it; the junction database is tuned for reads of {N}+1 bp (your reads: {READ_LENGTH} bp)." Reuse: `{STAR_INDEX}` = that directory; otherwise `{STAR_INDEX}` is empty. Tell the user, when reuse is chosen: at this revision the pipeline copies a given STAR index into its `work/` directory before aligning (about 30 GB for a human index, for every run), and this skill's verification run did not exercise the reuse path (unverified).
- Salmon: `{SALMON_INDEX}` is always empty: the pipeline always builds its own Salmon index from the transcripts it extracts from the GTF. A Salmon index from /nfcore-rnaseq-setup is built from Ensembl cDNA and lacks the GTF's non-coding transcripts, which would get no quantification (DTU, SUPPA2) without any error; its index version cannot show this, so it is never reused.

When the pipeline builds the STAR index of a human or mouse genome, STAR_GENOMEGENERATE needs about 32 GB of memory and an hour or more; the `nextflow.config` of Step 10 gives it 64 GB and 8 h. The pipeline also builds a decoy-aware Salmon index (SALMON_INDEX: 6 CPUs, 36 GB and 8 h from its process label); for a human or mouse genome this likely takes more than an hour (not measured by this skill's verification). Indexes built by the pipeline are not kept (`save_reference: false`), so a later run builds them again.

---

## Step 8 — Analyses (modules) and their settings

In the pipeline's own configuration every analysis module is switched on, so a module that is not mentioned in the params file runs anyway (the one exception at this revision is LeafCutter, which is off by default). This skill therefore writes every switch explicitly.

Ask (numbered; several can be chosen, for example "1, 3"): "Which analyses should run?"
1. **rMATS** — differential splicing events (skipped exon, alternative 5' and 3' splice sites, mutually exclusive exons, retained intron) from junction reads, for each contrast. Recommended. As an example (not a standard), Akerberg et al. 2022 kept rMATS events with 0 uncalled replicates, FDR < 0.1 (zebrafish) or FDR < 0.05 (human), and |IncLevelDifference| > 0.1 (from the parts of the paper's Methods that were accessible); choose your own cut-offs when filtering the results.
2. **SUPPA2** — PSI of local events (SE, SS, MX, RI, FL) and of isoforms from Salmon transcript abundance, with differential splicing between conditions. FASTQ input only.
3. **DEXSeq exon usage** — differential usage of exon bins within a gene (exon counts relative to the gene).
4. **edgeR exon usage** — the same question with edgeR on featureCounts exon counts (edgeR function in this revision: diffSpliceDGE glmQLFit glmQLFTest).
5. **DEXSeq transcript usage (DTU)** — changes in the share of each transcript within its gene (DRIMSeq filter, DEXSeq, stageR), from Salmon. FASTQ input only. ⚠️ `/bulk-rnaseq-pipeline` has the same DTU workflow (DRIMSeq → DEXSeq → stageR) on the Salmon output of an nf-core/rnaseq run; if it ran (or will run) there for these samples, do not choose it here: it is the same analysis twice. Choose DTU here only to get it in the same pipeline run from the FASTQ files; for event-level splicing questions (which exons or events change) use rMATS, and differential expression plus DTU on the Salmon output of an existing nf-core/rnaseq run stays in `/bulk-rnaseq-pipeline`.

Default (empty answer): 1 only. At least one analysis must be chosen. With BAM input, options 2 and 5 are not offered (they need Salmon quantification from reads), and `{RUN_SUPPA}` and `{RUN_DEXSEQ_DTU}` are `false`. Set `{RUN_RMATS}`, `{RUN_SUPPA}`, `{RUN_DEXSEQ_EXON}`, `{RUN_EDGER_EXON}` and `{RUN_DEXSEQ_DTU}` to `true` for the chosen analyses and `false` for all others. Show the choice as a table (analysis, on/off).

**MISO is not used**: in this pipeline it only draws sashimi plots for a short gene list (by default three human Ensembl IDs), it is no genome-wide splicing test, and MISO itself is unmaintained Python 2 software. This skill always writes `sashimi_plot: false`. For sashimi plots, use rmats2sashimiplot or ggsashimi on the BAMs afterwards.

**rMATS settings** (asked only when rMATS is chosen; otherwise the values below are written unchanged):
- Novel splice sites (numbered): 1. Annotated splice sites only (default) · 2. Also detect unannotated splice sites. `{RMATS_NOVEL}` = `false` (answer 1, and whenever rMATS is not chosen) or `true` (answer 2). With 2, rMATS also uses `rmats_min_intron_len` (50) and `rmats_max_exon_len` (500), written at their defaults.
- `rmats_splice_diff_cutoff` stays at 0.0001: it is the threshold of rMATS's null hypothesis (an inclusion difference larger than this counts as differential), not a reporting cut-off; filter the results by FDR and |IncLevelDifference| afterwards.
- `rmats_read_len` = `{READ_LENGTH}` (Step 6); `rmats_paired_stats` = `{PAIRED_DESIGN}` (Step 5).

**DTU filter values** (always written; used when DTU runs). The pipeline defaults (`min_samps_gene_expr` 4, `min_samps_feature_expr` 2, `min_samps_feature_prop` 2) fit one design size only. This skill uses the rule of Love et al. 2018 (as `/bulk-rnaseq-pipeline` does), with the counts kept in Step 5:
- `{MIN_SAMPS_GENE_EXPR}` = the number of samples in the samplesheet (distinct names), because the filter is applied once to all samples.
- `{MIN_SAMPS_FEATURE}` = the size of the smallest condition used in a contrast (for both `min_samps_feature_expr` and `min_samps_feature_prop`).
- `min_gene_expr` 10, `min_feature_expr` 10, `min_feature_prop` 0.1, `dtu_txi` `dtuScaledTPM` (the pipeline defaults).

Tell the user the six filter values and how the first three were derived (for example: "6 samples, smallest compared condition 3: `min_samps_gene_expr` 6, `min_samps_feature_expr` 3, `min_samps_feature_prop` 3").

**SUPPA2, DEXSeq and edgeR settings:** pipeline defaults, written explicitly; SUPPA2's paired test `diffsplice_paired` = `{PAIRED_DESIGN}` (its pipeline default `true` assumes paired samples).

**Salmon route.** `aligner: "star"` and `pseudo_aligner: "salmon"`: the STAR alignments feed rMATS, DEXSeq and edgeR; Salmon, run on the reads, feeds DTU and SUPPA2 once. (The pipeline's own default, `star_salmon`, quantifies with Salmon twice and runs DTU and SUPPA2 on both.) With FASTQ input, Salmon (index build and quantification) runs even when neither DTU nor SUPPA2 is chosen, because `pseudo_aligner` has no off value at this revision; the cost is mainly the Salmon index build (Step 7). With BAM input no Salmon step runs (verified for this revision).

**Module and option keys (written to the params file).** Step 11 inserts this block, with its placeholders filled, into the params file:
```yaml
# alignment and quantification
aligner: "star"
pseudo_aligner: "salmon"
# analysis modules: every switch explicit (pipeline defaults: on for all of them except leafcutter)
rmats: {RUN_RMATS}
dexseq_exon: {RUN_DEXSEQ_EXON}
edger_exon: {RUN_EDGER_EXON}
dexseq_dtu: {RUN_DEXSEQ_DTU}
suppa: {RUN_SUPPA}
sashimi_plot: false
isoformswitchanalyzer: false
leafcutter: false
# rMATS
rmats_read_len: {READ_LENGTH}
rmats_paired_stats: {PAIRED_DESIGN}
rmats_splice_diff_cutoff: 0.0001
rmats_novel_splice_site: {RMATS_NOVEL}
rmats_min_intron_len: 50
rmats_max_exon_len: 500
# DEXSeq exon usage
alignment_quality: 10
aggregation: true
save_dexseq_annotation: false
save_dexseq_plot: true
n_dexseq_plot: 10
# edgeR exon usage
save_edger_plot: true
n_edger_plot: 10
# DEXSeq transcript usage (DRIMSeq filter, DEXSeq, stageR)
dtu_txi: "dtuScaledTPM"
min_samps_gene_expr: {MIN_SAMPS_GENE_EXPR}
min_samps_feature_expr: {MIN_SAMPS_FEATURE}
min_samps_feature_prop: {MIN_SAMPS_FEATURE}
min_gene_expr: 10
min_feature_expr: 10
min_feature_prop: 0.1
ignore_tx_version: true
# SUPPA2
suppa_per_local_event: true
suppa_per_isoform: true
generateevents_pool_genes: true
generateevents_event_type: "SE SS MX RI FL"
generateevents_boundary: "S"
generateevents_threshold: 10
generateevents_exon_length: 100
psiperevent_total_filter: 0
diffsplice_local_event: true
diffsplice_isoform: true
diffsplice_method: "empirical"
diffsplice_area: 1000
diffsplice_lower_bound: 0
diffsplice_gene_correction: true
diffsplice_paired: {PAIRED_DESIGN}
diffsplice_alpha: 0.05
diffsplice_median: false
diffsplice_tpm_threshold: 0
diffsplice_nan_threshold: 0
clusterevents_local_event: true
clusterevents_isoform: true
clusterevents_dpsithreshold: 0.05
clusterevents_eps: 0.05
clusterevents_metric: "euclidean"
clusterevents_min_pts: 20
clusterevents_method: "DBSCAN"
# MISO sashimi plots (off: sashimi_plot is false), options at the pipeline defaults
miso_genes: "ENSG00000004961, ENSG00000005302, ENSG00000147403"
miso_read_len: 75
fig_height: 7
fig_width: 7
# IsoformSwitchAnalyzeR (off), options at the pipeline defaults
isoformswitchanalyzer_alpha: 0.05
isoformswitchanalyzer_dIF: 0.1
```
Every option of every module is in this block, also for modules that are switched off, at the pipeline default except the switches, `aligner`, `rmats_read_len`, `rmats_novel_splice_site`, the paired tests and the DTU sample counts above. Options without a default value (`gff_dexseq`, `suppa_tpm`, `clusterevents_sigthreshold`, `clusterevents_separation`, `miso_genes_file`) are not written: the pipeline then derives what it needs itself (for example the DEXSeq annotation from the GTF).

Never written: `rmats_variable_read_len` and `local_events` (not parameters; rMATS always uses variable read lengths here, and the SUPPA2 event types are `generateevents_event_type`); `max_cpus`, `max_memory` and `max_time` (not parameters of this revision, and Nextflow 26.04 stops at launch on an undeclared key: resource caps go into `nextflow.config`, Step 10); `salmon_index` (Step 7: the pipeline always builds its own Salmon index); no iGenomes `genome` key (FASTA and GTF are always given). `isoformswitchanalyzer` and `leafcutter` are parameters of this development revision and are written `false` (this skill does not offer them).

---

## Step 9 — Trimming and QC

Trimming (Trim Galore) and QC use the pipeline defaults; this skill does not change them and writes no trimming key. If the user wants other trimming settings, they edit the params file after Step 11 and should know that the run is then not the configuration this skill verified.

---

## Step 10 — MultiQC title, output directory, nextflow.config

**MultiQC title** (numbered): 1. `{SEQ_DATE}_{WD_NAME}` · 2. `{TODAY_YYMMDD}_{WD_NAME}` · 3. Custom. Store `{MULTIQC_TITLE}`.
**Output directory** (numbered): 1. `results/{TODAY_ISO}_{WD_NAME}` (default) · 2. `results/{TODAY_ISO}_{SHEET_PREFIX}` · 3. Custom. Store `{OUTDIR}`. Check it with `cd "{CWD}" && ls -A "{OUTDIR}"` (a relative `{OUTDIR}` is relative to `{CWD}`, where the job starts); if it exists and is not empty, ask (numbered): 1. use it anyway (the pipeline adds to it and may overwrite files) · 2. choose another.
A custom title or output directory may contain only letters, digits, `.`, `_`, `-` and (in the directory) `/`: a double quote or a backslash would break its double-quoted value in the params file, and a space breaks the shell lines. For any other character, tell the user "⚠️ `{VALUE}` contains a character the params file cannot hold; use only letters, digits, `.`, `_`, `-` (and `/` in the directory)." and ask again.

Check for an existing config: `ls "{CWD}/nextflow.config"`. If it exists, do not overwrite it: tell the user which selectors matter for rnasplice (the `withName` lines and the `resourceLimits` line below) so they can compare. If it does not exist, write it from the template.

**nextflow.config template.**
```nextflow
// nextflow.config — nf-core/rnasplice on SLURM + Singularity (generated by /nfcore-rnasplice-setup)
profiles {
    slurm {
        process {
            executor = 'slurm'
            queue = 'bcc'
            cpus = 2
            memory = '8 GB'
            time = '4h'

            withName: '.*:STAR_GENOMEGENERATE' {
                cpus = 8
                memory = '64 GB'
                time = '8h'
            }
            withName: '.*:STAR_ALIGN' {
                cpus = 8
                memory = '48 GB'
                time = '8h'
            }
            withName: '.*:RMATS_PREP' {
                cpus = 4
                memory = '16 GB'
                time = '8h'
            }
            withName: '.*:RMATS_POST' {
                cpus = 8
                memory = '32 GB'
                time = '16h'
            }
            withName: '.*:DEXSEQ_COUNT' {
                cpus = 2
                memory = '8 GB'
                time = '8h'
            }
            withName: '.*:DEXSEQ_EXON' {
                cpus = 8
                memory = '32 GB'
                time = '8h'
            }
            withName: '.*:DEXSEQ_DTU' {
                cpus = 8
                memory = '32 GB'
                time = '8h'
            }
            withName: '.*:SALMON_QUANT.*' {
                cpus = 8
                memory = '36 GB'
                time = '8h'
            }
        }
        executor {
            queueSize = 10
            submitRateLimit = '10/1min'
            pollInterval = '30s'
        }
    }
    singularity {
        singularity {
            enabled = true
            autoMounts = true
        }
    }
}

process {
    resourceLimits = [ cpus: 16, memory: '64 GB', time: '24h' ]
}

timeline { enabled = true; overwrite = true; file = "${params.outdir}/pipeline_info/execution_timeline.html" }
report   { enabled = true; overwrite = true; file = "${params.outdir}/pipeline_info/execution_report.html"   }
trace    { enabled = true; overwrite = true; file = "${params.outdir}/pipeline_info/execution_trace.txt"; fields = 'task_id,hash,native_id,name,status,exit,cpus,memory,time,realtime,peak_rss' }
dag      { enabled = true; overwrite = true; file = "${params.outdir}/pipeline_info/pipeline_dag.svg"        }
```

The selectors (`'.*:NAME'`, matching the process whatever its workflow prefix) use the process names of this skill's verification runs of nf-core/rnasplice dev-1b44723; every other process keeps the resources of the pipeline's own process labels. The Salmon selector ends in `.*` because the pipeline names its Salmon quantification processes `SALMON_QUANT_SALMON` and `SALMON_QUANT_STAR`; `'.*:SALMON_QUANT'` matched neither in the verification run. It requests 8 CPUs, 36 GB and 8 h: the memory and time of the pipeline's own tested process label for Salmon quantification (the label gives 6 CPUs; the decoy-aware Salmon index of a human or mouse genome is large). The resources of every selector are judgement, checked on the nf-core test data of this skill's verification only: they are not measured on real data. Fixed `withName` resources do not grow when the pipeline retries a task: if a task of a process with a selector fails with exit status 137 or 140 (memory or time limit), raise that selector's `memory` or `time` in nextflow.config (up to the `resourceLimits` caps) and resubmit as described under **Resuming a run** (Step 11). After the first real run, compare the selectors with `{OUTDIR}/pipeline_info/execution_trace.txt`, whose `cpus`, `memory` and `time` columns show what each task requested. `resourceLimits` caps every task at 16 CPUs, 64 GB and 24 h; this revision has no `max_cpus`, `max_memory` or `max_time` parameters (Step 8), so the caps live only here. `overwrite = true` on the four report files lets a run that failed early be started again.

---

## Step 11 — Params file, submission script and helper script

**Params file.** `{PARAMS_YAML}` = `{SHEET_PREFIX}_params.yaml` in `{CWD}`; if it exists, ask (numbered): 1. overwrite · 2. choose another filename. Build it from the template below: replace the `{MODULE_PARAMS}` line by the Step 8 block, fill every placeholder, and delete the `star_index` line when `{STAR_INDEX}` is empty (BAM input, or no reusable STAR index in Step 7: the pipeline then builds the index itself). No other line is deleted, and no placeholder may remain. There is no `salmon_index` line: the pipeline always builds its own Salmon index (Step 7).

**Params file template.**
```yaml
# nf-core/rnasplice {VERSION} parameters — generated by /nfcore-rnasplice-setup
input: "{SAMPLESHEET_CSV}"
contrasts: "{CONTRASTS_CSV}"
source: "{SOURCE}"
outdir: "{OUTDIR}"
multiqc_title: "{MULTIQC_TITLE}"
fasta: "{FASTA_PATH}"
gtf: "{GTF_PATH}"
gencode: false
star_index: "{STAR_INDEX}"
save_reference: false
{MODULE_PARAMS}
```

Typing rules: paths and strings are double-quoted; numbers and booleans are bare (`rmats_read_len: 101`, `rmats: true`); `multiqc_title` is always double-quoted, so that a numeric-looking title (for example `261001`) stays a string: Nextflow 26.04.6 stops at launch on a bare number for a string parameter (seen in this skill's verification). Why a params file: under Nextflow 26.04, values given on the command line (a double-dash option followed by a value) reach the parameter validator as strings, and numeric parameters were rejected ("Value is [string] but should be [number]"); a params file keeps the YAML types. So every pipeline parameter is in this file, and the launch line carries only `-params-file`.

**Submission script.** Write `nf-core_rnasplice_{VERSION_TAG}.sh` in `{CWD}` (if it exists, ask: 1. overwrite · 2. choose another filename):
```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=8G
#SBATCH -t 48:00:00
#SBATCH -p bcc
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={USER_EMAIL}
#SBATCH -o nf-core_rnasplice_{VERSION_TAG}.%j.log

cd "{CWD}" || { echo "ERROR: cannot change to {CWD}" >&2; exit 1; }
for f in nextflow.config {PARAMS_YAML} {SAMPLESHEET_CSV} {CONTRASTS_CSV}; do [ -s "$f" ] || { echo "ERROR: $f is missing or empty in {CWD}" >&2; exit 1; }; done

module add miniconda3/v4 || { echo "ERROR: cannot load module miniconda3/v4" >&2; exit 1; }
source /home/software/conda/miniconda3/bin/condainit || { echo "ERROR: cannot source condainit" >&2; exit 1; }
conda activate {CONDA_ENV} || { echo "ERROR: cannot activate conda environment {CONDA_ENV}" >&2; exit 1; }
module add singularity/3.10.4 || { echo "ERROR: cannot load module singularity/3.10.4" >&2; exit 1; }
command -v singularity >/dev/null || { echo "ERROR: singularity is not on PATH" >&2; exit 1; }
command -v nextflow >/dev/null || { echo "ERROR: nextflow is not in conda environment {CONDA_ENV}" >&2; exit 1; }

# Singularity images are kept here and reused by later runs (a preset value, for example from ~/.bashrc, is kept)
export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-$HOME/.singularity/cache}"
mkdir -p "$NXF_SINGULARITY_CACHEDIR" || { echo "ERROR: cannot create $NXF_SINGULARITY_CACHEDIR" >&2; exit 1; }

# Nextflow version range verified for nf-core/rnasplice dev-1b44723 by this skill (checked here, on the compute node)
NF_MIN="26.04.0"; NF_MAX_EXCL="none"
NF_VER=$(nextflow -version 2>/dev/null | awk '/version/ {for (i = 1; i <= NF; i++) if ($i ~ /^[0-9]+\.[0-9]+\.[0-9]+/) {print $i; exit}}')
ver_ge() { [ "$(printf '%s\n%s\n' "$1" "$2" | sort -V | head -n1)" = "$2" ]; }
[ -n "$NF_VER" ] || { echo "ERROR: could not read the Nextflow version" >&2; exit 1; }
ver_ge "$NF_VER" "$NF_MIN" || { echo "ERROR: Nextflow $NF_VER is older than $NF_MIN, the minimum for this pipeline revision" >&2; exit 1; }
if [ "$NF_MAX_EXCL" != none ] && ver_ge "$NF_VER" "$NF_MAX_EXCL"; then
  echo "ERROR: Nextflow $NF_VER is $NF_MAX_EXCL or newer; this pipeline revision was verified only with older versions (see the skill README)" >&2; exit 1
fi
echo "Nextflow $NF_VER (verified with 26.04.6)"

nextflow run nf-core/rnasplice -r {VERSION} -c nextflow.config -profile slurm,singularity -params-file {PARAMS_YAML}
```
The head job only coordinates the pipeline (2 CPUs, 8 GB) but must outlive every task, hence 48 h. Every path the job uses (nextflow.config, the params file, the sheets, the output directory) is relative to `{CWD}`: it changes to `{CWD}` first, so it works wherever it is submitted from, and stops with a clear message when one of its files is missing or empty. The module and conda lines are the ones that worked in every job of this skill's verification, each stopping the job when it fails (under Lmod a failed `module add` silently breaks the later ones). Nextflow 26.04.0 or newer is required and there is no upper bound; the job reads the version from `nextflow -version` on the compute node and stops with a clear message before the pipeline starts when it is too old. The script fetches nothing itself: Nextflow downloads the pinned revision from GitHub the first time a job runs it, also when `~/.bashrc` sets `NXF_OFFLINE=TRUE` (verified for this revision; the compute nodes have internet access), and the Singularity images go into `NXF_SINGULARITY_CACHEDIR`.

**Resuming a run.** If the head job stops at its 48 h limit (SLURM state TIMEOUT), or after you raised the resources of a selector in nextflow.config, add ` -resume` by hand at the end of the `nextflow run` line of `nf-core_rnasplice_{VERSION_TAG}.sh` and submit it again: Nextflow reuses the finished tasks from `work/` in `{CWD}`. `-resume` is a Nextflow option, not a pipeline parameter; the wizard never writes it into the generated script.

**Genome download helper.** Only when Step 7 found a missing Ensembl FASTA or GTF: write `download_genome_{REF_TAG}.sh` in `{CWD}` (URLs verified in this session, Step 7). It runs on a compute node, keeps the `.gz` files next to the decompressed ones, can be re-run safely (an interrupted download resumes), and removes no file; a download that is not a valid gzip file is moved aside as `.gz.invalid`:
```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=8G
#SBATCH -t 4:00:00
#SBATCH -p bcc
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={USER_EMAIL}
#SBATCH -o download_genome_{REF_TAG}.%j.log
set -u
GENOME_DIR="{GENOME_DIR}"
case "$GENOME_DIR" in /?*) ;; *) echo "ERROR: GENOME_DIR must be an absolute path: '$GENOME_DIR'" >&2; exit 1 ;; esac
mkdir -p "$GENOME_DIR" || { echo "ERROR: cannot create $GENOME_DIR" >&2; exit 1; }
fetch() {  # fetch <url of a .gz file> <decompressed target path>
  local url=$1 out=$2
  case "$out" in /?*) ;; *) echo "ERROR: the target must be an absolute path: '$out'" >&2; exit 1 ;; esac
  if [ -s "$out" ]; then echo "present: $out"; return 0; fi
  if ! { [ -s "$out.gz" ] && gzip -t "$out.gz" 2>/dev/null; }; then
    wget -c -O "$out.gz.part" "$url" || { echo "ERROR: download failed for $url (partial file kept for resuming: $out.gz.part); run the helper again" >&2; exit 1; }
    gzip -t "$out.gz.part" 2>/dev/null || { if mv -f "$out.gz.part" "$out.gz.invalid"; then echo "ERROR: $url is not a valid gzip file (moved to $out.gz.invalid); run the helper again" >&2; else echo "ERROR: $url is not a valid gzip file, and $out.gz.part could not be moved aside: rename it by hand, then run the helper again" >&2; fi; exit 1; }
    mv -f "$out.gz.part" "$out.gz" || { echo "ERROR: cannot move $out.gz.part to $out.gz" >&2; exit 1; }
  fi
  gunzip -c "$out.gz" > "$out.part" && [ -s "$out.part" ] && mv -f "$out.part" "$out" \
    || { echo "ERROR: decompression failed for $out.gz ($out.part is incomplete); run the helper again" >&2; exit 1; }
  echo "ready: $out"
}
fetch "{ENSEMBL_FASTA_URL}" "{FASTA_PATH}"
fetch "{ENSEMBL_GTF_URL}" "{GTF_PATH}"
echo "genome files ready in $GENOME_DIR"
```

Show every written file in full. Never run these scripts on the login node, and never run `nextflow`, `conda`, Java, R or Python there to test them: they are submitted with `sbatch` and run on compute nodes. Tell the user:
```
Params file written: {PARAMS_YAML}
Script written: nf-core_rnasplice_{VERSION_TAG}.sh
To submit:  sbatch nf-core_rnasplice_{VERSION_TAG}.sh
```
If the genome download helper was written, the pipeline must wait for it; submit both in one shell (or one Bash call), from `{CWD}`:
```bash
jid=$(sbatch --parsable download_genome_{REF_TAG}.sh) && sbatch --dependency=afterok:$jid nf-core_rnasplice_{VERSION_TAG}.sh
```
