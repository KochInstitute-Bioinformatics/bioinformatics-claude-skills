# nf-core/rnavar Pipeline Setup Skill

You are helping the user set up and submit an nf-core/rnavar RNA-seq variant-calling pipeline on an HPC cluster using SLURM and Singularity. rnavar implements the GATK best-practices RNA-seq short-variant workflow: 2-pass STAR → MarkDuplicates → SplitNCigarReads → base recalibration (BQSR) → HaplotypeCaller → soft variant filtering → optional annotation (SnpEff / VEP). Walk through each step below in order, asking the user only what you need and performing automated steps silently.

---

## Step 0 — Establish working directory

Before anything else, run `pwd` to record the current working directory (`{CWD}`). All scripts, the samplesheet CSV, `nextflow.config`, and the output directory will be written here. Inform the user: "Working directory: `{CWD}`"

Also derive:
- `{WD_NAME}` = basename of `{CWD}`
- `{TODAY_YYMMDD}` = today's date formatted as `YYMMDD`

---

## Step 1 — Email address

Ask the user: "What email address should SLURM use for pipeline notifications (completion/failure)?"

Do **not** pre-suggest or pre-fill any email address. Wait for the answer before Step 2.

---

## Step 2 — Conda environment

Ask the user: "Which conda environment should be activated for Nextflow?"

Do **not** pre-suggest any environment name. Ask this as a separate question after Step 1 — never combine Steps 1 and 2 in one message.

---

## Step 3 — Latest nf-core/rnavar version

Fetch the latest release with the `WebFetch` tool against:
```
https://api.github.com/repos/nf-core/rnavar/releases/latest
```
Extract `tag_name` (strip any leading `v`) and use it as `{VERSION}`. Tell the user which version you found. The parameter names in this skill were verified against release 1.3.0; if `{VERSION}` is a newer major or minor version, fetch `https://raw.githubusercontent.com/nf-core/rnavar/{VERSION}/nextflow_schema.json` and confirm every parameter this skill emits still exists before writing any script.

**Important:** Do NOT use the `gh` CLI — it is not installed on this HPC cluster. Always use `WebFetch`.

---

## Step 4 — Raw data location and samplesheet

**Before asking the user, scan the working directory for input files:**

```bash
find {CWD} -name "*.fastq.gz" -o -name "*.fq.gz" | head -30
find {CWD} -name "*.bam" -o -name "*.cram" | head -10
```

- FASTQ files found: list the unique directories and ask "I found FASTQ files in: `{FOUND_DIRS}`. Use this directory, or specify another?"
- Only BAM/CRAM found: ask whether these are STAR-aligned, duplicate-marked BAMs (rnavar accepts them as input and skips alignment) and use the BAM samplesheet form below.
- Nothing found: ask for the full path to the input files.

**Paired-end detection:** filenames containing `_R1_`/`_R2_`, `_1.fastq.gz`/`_2.fastq.gz`, or `_1_sequence`/`_2_sequence` → paired-end. Otherwise single-end: warn "Single-end data is accepted by rnavar but gives weaker variant calls than paired-end; continue?" Detect the sequencing date from a leading 6-digit `YYMMDD` filename prefix as `{SEQ_DATE}` (fall back to today's date).

**Path style:** if the input directory is inside `{CWD}`, use paths relative to `{CWD}`; otherwise absolute paths.

**Samplesheet forms (exactly one file type per sample):**
- FASTQ: `sample,fastq_1,fastq_2` (rnavar has no `strandedness` column). Find each R1, pair it with its R2 by substituting `_R1_`→`_R2_`, `_1.`→`_2.`, or `_1_sequence`→`_2_sequence`; warn and leave `fastq_2` empty if R2 is missing. Sample name = filename up to `_S\d+`, `_R1`, or `_1_sequence`.
- BAM: `sample,bam,bai`. CRAM: `sample,cram,crai`.
- **Never mix types for one sample.** Supplying FASTQ files and a BAM/CRAM file for the same sample makes the pipeline error; check for this and stop with a clear message.

**Sample-name sanitisation (always, before showing the user):**
- Replace every `-` with `_`; replace spaces, `/`, `(`, `)` and other special characters with `_`. Note substitutions in the preview.
- Check uniqueness after sanitisation. On collision warn: "⚠️ Name collision '{NAME}': in nf-core/rnavar, rows with the same sample name are treated as lanes of one sample and their reads are merged before alignment. Provide distinct names if these are different samples." Then go to custom naming.

**Review and naming — order is mandatory:**
1. Show the full samplesheet (all rows) as a table.
2. Ask about names (numbered): 1. Use auto-generated names · 2. Provide custom names. For custom names, show numbered auto names next to filenames, ask for a plain-language description, build the mapping, show it as an auto→new table and ask "Does this mapping look correct?" Validate custom names: no `-`; no duplicates (same warning as above).
3. Ask for the samplesheet filename (numbered): 1. `{SEQ_DATE}_{WD_NAME}_samplesheet.csv` · 2. `{TODAY_YYMMDD}_{WD_NAME}_samplesheet.csv` · 3. Custom.
4. Write the file only after names and filename are confirmed. Store as `{SAMPLESHEET_CSV}`.

---

## Step 5 — Read length

rnavar uses `--read_length` to set STAR's `sjdbOverhang` (`read_length − 1`), and its default is 150, which is wrong for most other libraries. **Always detect and always pass it.**

```bash
zcat {FASTQ_FILE} | awk 'NR%4==2 {print length($0)}' | head -n 1000 | sort -n | uniq -c | sort -rn | head -3
```

Run this on the first FASTQ of several different samples (up to 5). If lengths differ between samples, report the distribution, use the most common read length as `{READ_LENGTH}`, and warn that `sjdbOverhang` is tuned to it. For BAM/CRAM input, ask the user for the read length. Tell the user: "Detected read length {READ_LENGTH} bp → sjdbOverhang {READ_LENGTH − 1}." Store `{READ_LENGTH}`.

---

## Step 6 — Organism and genome files

Ask:
1. "What organism is this data from? (e.g. mouse, human)"
2. "What is the base directory where genome files and indexes are stored?"

Use this folder convention (shared with other nf-core skills so FASTA/GTF are reused, never re-downloaded):

```
{genome_base}/{organism}/{assembly}_ens{version}/
├── {FASTA}.fa                    ← primary assembly FASTA
├── {GTF}.gtf                     ← annotation GTF
└── index/
    └── star_rnavar_sjdb{N-1}/    ← STAR index built for THIS read length
```

Store: `{GENOME_DIR}` = `{genome_base}/{organism}/{assembly}_ens{version}`, `{FASTA_PATH}` and `{GTF_PATH}` (full paths of the FASTA and GTF files), `{ORGANISM}`, `{ASSEMBLY}` and `{ENS_VERSION}`. To determine the version: if `{genome_base}/{organism}/` already contains `{assembly}_ens{N}` directories, use the highest N unless the user asks otherwise; if none exist, use the latest Ensembl release from the `https://ftp.ensembl.org/pub/current_README` fetch.

- Mouse: assembly GRCm39, FASTA `Mus_musculus.GRCm39.dna.primary_assembly.fa`, GTF `Mus_musculus.GRCm39.{version}.gtf`, directory `{genome_base}/mouse/mm39_ens{version}/`.
- Human: assembly GRCh38, FASTA `Homo_sapiens.GRCh38.dna.primary_assembly.fa`, GTF `Homo_sapiens.GRCh38.{version}.gtf`, directory `{genome_base}/human/hg38_ens{version}/`.
- Other organisms: ask for the FASTA and GTF paths and the directory that will hold indexes (`{GENOME_DIR}`); skip only the download/version checks in this step and still apply the STAR index rule below.

**Existing FASTA/GTF:** if present, report the paths and reuse them. If missing, fetch the latest Ensembl release from `https://ftp.ensembl.org/pub/current_README`, and generate the download commands in the helper script (Step 12).

**GTF source:** inspect `grep -v "^#" {GTF_PATH} | head -3`. Gene IDs with a version suffix (`ENSG00000000003.15`) indicate GENCODE; without one, Ensembl. Report which it is, but **no flag is emitted** — rnavar's schema has no `gencode` parameter.

**STAR index — always its own index per read length.** `sjdbOverhang` is fixed when the index is built, and an index made for another read length (for example one built by another pipeline) would be wrong here. Set `{SJDB_OVERHANG}` = `{READ_LENGTH} − 1` and look for `{GENOME_DIR}/index/star_rnavar_sjdb{SJDB_OVERHANG}/`:
- Present and non-empty (`SA`, `Genome`, `sjdbList.out.tab` exist): use it, `{STAR_INDEX}` = that path.
- Missing: {STAR_INDEX} is the same path; the helper script (Step 12) builds it there and must run before the pipeline.

`--star_index '{STAR_INDEX}'` is always passed so rnavar never rebuilds the index inside the workflow.

---

## Step 7 — Known sites for base recalibration

rnavar passes the known-sites files directly to GATK BaseRecalibrator and **does not skip base recalibration automatically** — if they are missing the run fails late, after alignment. So always resolve this now. Ask (numbered):

1. **Use known-sites VCFs** — add `--dbsnp {DBSNP}`, `--dbsnp_tbi {DBSNP}.tbi`, `--known_indels {INDELS}`, `--known_indels_tbi {INDELS}.tbi`.
   - Human: the GATK resource-bundle dbSNP and Mills/1000G known-indels VCFs for the matching assembly.
   - Mouse: Mouse Genomes Project variants (SNPs and indels) for GRCm39.
   - Check whether the files already exist under `{GENOME_DIR}/known_sites/`; if so reuse them and verify the `.tbi` indexes exist.
   - If they do not exist, **resolve the resource URLs at run time**: use `WebFetch` on the current GATK resource-bundle page (human) or the Mouse Genomes Project / Ensembl variation FTP listing (mouse), show the exact URLs to the user, and download only after they confirm. Never type a URL from memory. Add the download and `tabix -p vcf` steps to the helper script (Step 12).
2. **Skip base recalibration** — add `--skip_baserecalibration`. Tell the user the trade-off: base qualities are not recalibrated, which is slightly less accurate but is the right choice for organisms without a curated variant set, or to get a first result quickly.

Store the result as `{KNOWN_SITES_LINES}`: the four `--dbsnp`/`--known_indels` lines for option 1, or the single line `--skip_baserecalibration \` for option 2.
