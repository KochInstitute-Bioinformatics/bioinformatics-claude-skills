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

**Scan the working directory first:**
```bash
find {CWD} \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | head -50
find {CWD} -name "*.bam" | head -20
```
Ask (numbered), mentioning what the scan found: "Which input do you want to start from?" 1. FASTQ files (default; the pipeline aligns them with STAR and quantifies them with Salmon) · 2. Genome BAM files from a splice-aware aligner (for example the STAR BAMs of an nf-core/rnaseq run). Store `{SOURCE}` = `fastq` or `genome_bam`.

### Step 4a — FASTQ input

- FASTQ files found: list the unique directories and ask "I found FASTQ files in: `{FOUND_DIRS}`. Use this directory, or specify another?" Nothing found: ask for the full path. Store `{FASTQ_DIR}`.
- **Paired-end detection:** filenames containing `_R1_`/`_R2_`, `_1.fastq.gz`/`_2.fastq.gz`, `_1.fq.gz`/`_2.fq.gz`, or `_1_sequence`/`_2_sequence` → paired-end; otherwise single-end. Store `{LAYOUT}` = `paired` or `single`. rMATS needs one strandedness and one read type for all samples (the pipeline stops with "Cannot run rMats with mixed single and paired end samples"), so if the directory mixes paired-end and single-end files, stop and ask which files to use.
- **Sequencing date:** a leading 6-digit `YYMMDD` filename prefix → `{SEQ_DATE}` (fall back to `{TODAY_YYMMDD}`).
- **Path style:** if `{FASTQ_DIR}` is inside `{CWD}`, use paths relative to `{CWD}`; otherwise absolute paths.
- **Rows:** find each R1 and pair it with its R2 by substituting `_R1_`→`_R2_`, `_1.`→`_2.` or `_1_sequence`→`_2_sequence`; if an R2 is missing, warn and stop (one read type for all samples). Single-end rows leave `fastq_2` empty. Sample name = filename up to `_S\d+`, `_R1`, `_1.f` or `_1_sequence`. FASTQ files must end in `.fastq.gz` or `.fq.gz`; the pipeline accepts nothing else.

**Sample-name sanitisation (always, before showing the user):** replace `-`, spaces, `/`, `(`, `)` and every other character outside `A-Za-z0-9_` with `_`; if a name then starts with a digit, prefix `S` (R renames columns that start with a digit, which breaks the matching of sample names in the DEXSeq and edgeR steps). Note every substitution in the preview. Then check uniqueness. On a collision warn: "⚠️ Name collision '{NAME}': in nf-core/rnasplice, rows with the same sample name are one sample, and their reads are merged before alignment." and ask (numbered): 1. the same sample (lanes or technical replicates; the pipeline merges their reads) · 2. different samples — rename them (go to custom naming). Lane files (`X_S1_L001_R1_001`, `X_S1_L002_R1_001`) collide by design; option 1 is expected for them.

**Review and naming — order is mandatory:**
1. Show all rows as a table (sample, fastq_1, fastq_2).
2. Ask (numbered): 1. Use auto-generated names · 2. Provide custom names. For custom names: show numbered auto names next to the filenames, ask for a plain-language description, build the mapping, show it as an auto → new table and ask "Does this mapping look correct?" Custom names follow the same rule (letters, digits and `_`, starting with a letter); duplicates only where the user chose option 1 above.

The rows stay in memory: Step 5 adds the condition column (the samplesheet columns are `sample,fastq_1,fastq_2,strandedness,condition`) and writes the file.

**Strandedness is asked, never guessed.** rnasplice has no `auto` strandedness, and rMATS needs one strandedness and one read type for all samples. Ask (numbered): "Which strandedness does the library have?" 1. unstranded · 2. forward · 3. reverse · 4. I don't know. Explain: *reverse* = read 1 comes from the strand opposite to the transcript (dUTP libraries such as Illumina TruSeq Stranded and NEBNext Ultra II Directional — the most common stranded libraries); *forward* = read 1 comes from the transcript strand; *unstranded* = no strand information (for example TruSeq non-stranded). A wrong value makes rMATS discard or misassign junction reads without any error. If the user answers 4, or has an nf-core/rnaseq output directory for the same samples, use the block below. Never continue with "I don't know". Store `{STRANDEDNESS}` (one value for every row).

**Strandedness from an existing nf-core/rnaseq run.** Ask for that run's output directory (`{RNASEQ_OUTDIR}`), list its RSeQC files with `find {RNASEQ_OUTDIR} -name "*.infer_experiment.txt" | head -50`, define the function below in the Bash tool and run it on each file (small text files; light work on the login node). Show a table file → result.
```bash
infer_strandedness() {
  # usage: infer_strandedness <file.infer_experiment.txt>; prints forward, reverse, unstranded or unclear
  awk -F': ' '
    /explained by "1\+\+,1--,2\+-,2-\+"|explained by "\+\+,--"/ { f = $2 + 0; nf = 1 }
    /explained by "1\+-,1-\+,2\+\+,2--"|explained by "\+-,-\+"/ { r = $2 + 0; nr = 1 }
    END {
      if (!nf || !nr) { print "unclear"; exit }
      if (f >= 0.8) print "forward"
      else if (r >= 0.8) print "reverse"
      else if (f >= 0.3 && f <= 0.7 && r >= 0.3 && r <= 0.7) print "unstranded"
      else print "unclear"
    }' "$1"
}
```
- All files give the same `forward`, `reverse` or `unstranded`: propose it and ask the user to confirm.
- Any file gives `unclear`, or the files disagree: show the fractions (`grep "explained by" FILE`), explain that the library type cannot be read reliably, and ask the user to check the kit; do not continue until they answer 1, 2 or 3.
- The user answered 1-3 and a result disagrees: show both and ask which is right before continuing.

### Step 4b — Genome BAM input

Genome BAMs must come from a splice-aware aligner (STAR, as in nf-core/rnaseq: `{RNASEQ_OUTDIR}/star_salmon/{sample}.markdup.sorted.bam`), be coordinate-sorted, and be aligned to the same FASTA and GTF that Step 7 gives the pipeline (same contig names). No `.bai` index is needed: the pipeline sorts and indexes the BAMs itself. Find them with `find {DIR} -name "*.bam" ! -name "*toTranscriptome*"`, derive sample names by removing `.markdup.sorted.bam`, `.sorted.bam` or `.bam`, apply the sanitisation rule of Step 4a, and use one BAM per sample (BAM rows are never merged). Ask the strandedness question of Step 4a (with the nf-core/rnaseq block when the BAMs come from such a run) and the read type (numbered): 1. paired-end · 2. single-end; store `{STRANDEDNESS}` and `{LAYOUT}`. Then apply the BAM input rule.

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
