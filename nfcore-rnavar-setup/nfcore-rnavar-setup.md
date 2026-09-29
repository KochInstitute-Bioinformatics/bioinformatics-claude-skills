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
- Check uniqueness after sanitisation. On collision warn: "⚠️ Name collision '{NAME}': in nf-core/rnavar, rows with the same sample name are treated as lanes of one sample and their reads are merged before alignment. Provide distinct names if these are different samples." Then ask (numbered): 1. these are lanes of the same sample — keep the duplicate name (rows are merged before alignment) · 2. these are different samples — rename them (go to custom naming). Normal Illumina lane files (`X_S1_L001_R1_001`, `X_S1_L002_R1_001`) collide by design, so option 1 is expected for them.

**Review and naming — order is mandatory:**
1. Show the full samplesheet (all rows) as a table.
2. Ask about names (numbered): 1. Use auto-generated names · 2. Provide custom names. For custom names, show numbered auto names next to filenames, ask for a plain-language description, build the mapping, show it as an auto→new table and ask "Does this mapping look correct?" Validate custom names: no `-`; duplicates only where the user chose option 1 above (validation must allow the deliberate duplicates); any other duplicate triggers the same warning and choice.
3. Ask for the samplesheet filename (numbered): 1. `{SEQ_DATE}_{WD_NAME}_samplesheet.csv` · 2. `{TODAY_YYMMDD}_{WD_NAME}_samplesheet.csv` · 3. Custom.
4. Write the file only after names and filename are confirmed; first check whether `{SAMPLESHEET_CSV}` already exists and, if so, ask (numbered): 1. overwrite · 2. choose another filename — before writing. Store as `{SAMPLESHEET_CSV}`.

---

## Step 5 — Read length

rnavar uses the `read_length` parameter to set STAR's `sjdbOverhang` (`read_length − 1`), and its default is 150, which is wrong for most other libraries. **Always detect it and always write it to the params file.**

```bash
zcat {FASTQ_FILE} | awk 'NR%4==2 {print length($0)}' | head -n 1000 | sort -n | uniq -c | sort -rn | head -3
```

Run this on the first FASTQ of several different samples (up to 5). If lengths differ between samples, report the distribution, use the most common read length as `{READ_LENGTH}`, and warn that `sjdbOverhang` is tuned to it. For BAM/CRAM input, ask the user for the read length. Tell the user: "Detected read length {READ_LENGTH} bp → sjdbOverhang {READ_LENGTH − 1}." Store `{READ_LENGTH}` (written as `read_length: {READ_LENGTH}` in Step 11).

---

## Step 6 — Organism and genome files

Ask, in this order:
1. "What organism is this data from? (e.g. mouse, human)"
2. (numbered): "1. Ensembl release in the standard folder (default) · 2. Custom reference — I already have a FASTA and GTF".
3. The option-specific questions below. Ask the base directory only for option 1; option 2 asks for its own paths instead.

- **Option 1 (Ensembl):** ask "What is the base directory where genome files and indexes are stored?" (`{genome_base}`); follow the folder convention and version logic below; set `{REF_TAG}` = `{ASSEMBLY}_ens{ENS_VERSION}` and `{FASTA_SOURCE}` = `{FASTA_PATH}` (Step 6 option 1 binds `{FASTA_SOURCE}` = `{FASTA_PATH}`, since the Ensembl FASTA is never gzipped).
- **Option 2 (Custom reference):** ask for the FASTA path, the GTF path and the directory that will hold indexes (`{GENOME_DIR}`); skip the Ensembl download and version logic entirely; set `{FASTA_PATH}`, `{GTF_PATH}`, `{GENOME_DIR}` from the answers and `{REF_TAG}` = `custom_{WD_NAME}`. Also set `{FASTA_SOURCE}` = the original FASTA path the user gave, gzipped or not (the Step 7 contig check on the login node reads this file, because the decompressed `{FASTA_PATH}` does not exist yet when the wizard runs). Keep the per-read-length rule: still apply the STAR index rule below, and report the GTF source as usual. If the FASTA or GTF ends in `.gz`, the STAR helper decompresses it with `gunzip -c file.gz > {GENOME_DIR}/<name>` and `{FASTA_PATH}`/`{GTF_PATH}` point at the decompressed copies. The user's originals are never modified.

`{REF_TAG}` is the tag used in helper script names (Step 12).

For option 1, use this folder convention (shared with other nf-core skills so FASTA/GTF are reused, never re-downloaded):

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
- Other organisms: use option 2 (Custom reference).

**Existing FASTA/GTF:** if present, report the paths and reuse them. If missing, fetch the latest Ensembl release from `https://ftp.ensembl.org/pub/current_README`, and generate the download commands in the helper script (Step 12).

**GTF source:** inspect `grep -v "^#" {GTF_PATH} | head -3`. Gene IDs with a version suffix (`ENSG00000000003.15`) indicate GENCODE; without one, Ensembl. Report which it is, but **no flag is emitted** — rnavar's schema has no `gencode` parameter.

**STAR index — always its own index per read length.** `sjdbOverhang` is fixed when the index is built, and an index made for another read length (for example one built by another pipeline) would be wrong here. Set `{SJDB_OVERHANG}` = `{READ_LENGTH} − 1` and look for `{GENOME_DIR}/index/star_rnavar_sjdb{SJDB_OVERHANG}/`:
- Present and non-empty (`SA`, `Genome`, `sjdbList.out.tab` exist): use it, `{STAR_INDEX}` = that path.
- Missing: {STAR_INDEX} is the same path; the helper script (Step 12) builds it there and must run before the pipeline.

The `star_index` key is always written to the params file (Step 11) so rnavar never rebuilds the index inside the workflow.

---

## Step 7 — Known sites for base recalibration

rnavar passes the known-sites files directly to GATK BaseRecalibrator and **does not skip base recalibration automatically** — if they are missing the run fails late, after alignment. So always resolve this now. Ask (numbered):

1. **Use known-sites VCFs** — write the keys `dbsnp: "{DBSNP}"`, `dbsnp_tbi: "{DBSNP}.tbi"`, `known_indels: "{INDELS}"`, `known_indels_tbi: "{INDELS}.tbi"` (double-quoted paths). `{DBSNP}` and `{INDELS}` are the FINAL post-guard paths (the `.renamed.vcf.gz` names when a rename is expected, see Step 12).
   - Human: the GATK resource-bundle dbSNP and Mills/1000G known-indels VCFs for the matching assembly.
   - Mouse: Mouse Genomes Project variants (SNPs and indels) for GRCm39.
   - Check whether the files already exist under `{GENOME_DIR}/known_sites/`; if so reuse them and verify the `.tbi` indexes exist.
   - The pipeline schema requires bgzipped `.vcf.gz` files with `.tbi` indexes (`dbsnp` must match `.vcf.gz`, `dbsnp_tbi` must match `.vcf.gz.tbi`). Some sources ship plain `.vcf` (with `.vcf.idx`); those must be bgzipped and re-indexed with `tabix`, never passed as-is.
   - If they do not exist, **resolve the resource URLs at run time**: use `WebFetch` on the current GATK resource-bundle page (human) or the Mouse Genomes Project / Ensembl variation FTP listing (mouse), show the exact URLs to the user, and download only after they confirm. Never type a URL from memory. Add the download, `bgzip` and `tabix -f -p vcf` steps to the helper script (Step 12).
   - **Mandatory contig-name check, before the known-sites choice is finalised.** The FASTA from Step 6 is Ensembl-named (`1` ... `MT`) for option 1; a custom FASTA may use either style. GATK resource-bundle hg38 VCFs use `chr1` ... `chrM`. A mismatch makes GATK BaseRecalibrator stop with "incompatible contigs" only after alignment, MarkDuplicates and SplitNCigarReads have already run. On the login node (no `tabix` there), compare the first VCF contig with the first FASTA header (of `{FASTA_SOURCE}`, which may be gzipped; Step 12's guard runs later on the compute node and uses `{FASTA_PATH}`):
     ```bash
     VCF_CONTIG=$(zcat -f FILE | awk '!/^#/{print $1; exit}')
     FASTA_CONTIG=$(zcat -f {FASTA_SOURCE} | awk '/^>/{sub(/^>/,""); print $1; exit}')
     [ -n "$VCF_CONTIG" ] && [ -n "$FASTA_CONTIG" ] || echo "EMPTY contig — stop and resolve before continuing"
     ```
     Treat an empty value as a mismatch. If the FASTA does not exist yet (option 1, Ensembl download pending), skip the wizard-side check and rely on the helper's contig guard. If the VCFs already exist, this check runs in the wizard now. If they are NOT downloaded yet, tell the user that the helper's contig guard performs the check after download (see `prepare_known_sites_{REF_TAG}.sh` in Step 12), renames the contigs when a fix is possible, and stops with a non-zero exit on an unresolvable mismatch, and that the pipeline must therefore be submitted with `--dependency=afterok` on the helper job(s) (Step 11) so a failed guard prevents the pipeline from starting. Mouse Genomes Project VCFs use Ensembl-style contig names, so a rename is normally not needed for mouse, but the guard still runs. On a mismatch found in the wizard, either (i) prefer Ensembl-named variation VCFs that match the Ensembl FASTA, or (ii) fix them with the helper's rename step (`bcftools annotate --rename-chrs`, see Step 12) — the wizard also generates `prepare_known_sites_{REF_TAG}.sh` in this case, even though the VCFs already exist, so the rename and guard run before the pipeline. The skill must never proceed with mismatched contigs; if neither fix is possible, offer option 2 (skip base recalibration) instead.
2. **Skip base recalibration** — write `skip_baserecalibration: true`. Tell the user the trade-off: base qualities are not recalibrated, which is slightly less accurate but is the right choice for organisms without a curated variant set, or to get a first result quickly.

Store the result as `{KNOWN_SITES_PARAMS}`: the four key lines for option 1, or the single line `skip_baserecalibration: true` for option 2.

---

## Step 8 — Variant calling options

Ask each as a numbered choice. **Omit any key whose value equals the pipeline default** — only write a key to the params file when the user changes it. Every option below is written as a key-value line.

**a) Duplicates.** 1. Keep duplicates marked (default — omit) · 2. Remove duplicates — write `remove_duplicates: true`.

**b) STAR two-pass.** Two-pass mapping (`star_twopass`) is on by default and recommended for calling; keep it. Write `star_twopass: false` only if the user explicitly asks to disable it.

**c) Calling and filtering thresholds.** Ask: 1. Pipeline defaults (recommended; writes nothing) · 2. Customise. If customising, ask for each and write only values that differ from the default: `gatk_hc_call_conf: <int>` (default 20), `gatk_vf_qd_filter: <number>` (2), `gatk_vf_fs_filter: <number>` (30), `gatk_vf_window_size: <int>` (35), `gatk_vf_cluster_size: <int>` (3). Offer `skip_variantfiltration: true` if the user wants unfiltered calls.

**d) gVCFs.** Ask: "Will you jointly call variants across samples later?" 1. No (omit) · 2. Yes — write `generate_gvcf: true`.

**e) Large chromosomes.** Only if the genome has chromosomes longer than 512 Mb (not human or mouse): write `bam_csi_index: true` and tell the user it disables variant filtration.

**f) Save intermediates.** 1. No (omit) · 2. Yes — write `save_align_intermeds: true` (recommended if the BAMs will feed allele-specific expression analysis).

Collect the emitted key lines as `{VARIANT_PARAMS}`.

---

## Step 9 — Optional variant annotation

Ask (numbered): 1. No annotation (default — omit `tools`) · 2. SnpEff · 3. VEP · 4. Both merged. Write the key `tools: "snpeff"`, `tools: "vep"`, or `tools: "merge"` respectively. Set `{ANNOTATION_TOOL}`.

If annotation is chosen, the caches must be on disk **before** submission:
- Ask for existing cache directories → the keys `snpeff_cache: "{DIR}"` and/or `vep_cache: "{DIR}"`, plus the matching identifiers `snpeff_db: "{SNPEFF_DB}"`, `vep_genome: "{VEP_GENOME}"`, `vep_species: "{VEP_SPECIES}"`, `vep_cache_version: "{VEP_CACHE_VERSION}"` (all quoted strings; ask; do not guess versions).
- If a cache is missing, do **not** silently add `download_cache`: that option needs internet from compute nodes, which may not be available, and the job then stalls without a clear error. Instead tell the user this, and offer to generate a pre-download helper script (Step 12) that they run where internet is available.

Note for the assistant: `annotation_cache` appears on the rnavar usage page but is not a parameter in the rnavar schema — never emit it. Write `download_cache: true` only if the user explicitly chooses it after being warned.

Collect the emitted key lines as `{ANNOTATION_PARAMS}`.

---

## Step 10 — MultiQC title, output directory, nextflow.config

**MultiQC title** (numbered): 1. `{SEQ_DATE}_{WD_NAME}` · 2. `{TODAY_YYMMDD}_{WD_NAME}` · 3. Custom. Store as `{MULTIQC_TITLE}`.
**Output directory** (numbered): same three options. Store as `{OUTDIR}`.

Check for an existing config: `ls nextflow.config`. **If it exists, do not overwrite it** — instead tell the user the selectors that matter for rnavar (below) and that they can compare them. If it does not exist, write:

```nextflow
// nextflow.config — nf-core/rnavar on SLURM + Singularity
profiles {
    slurm {
        process {
            executor = 'slurm'
            queue = 'bcc'
            cpus = 2
            memory = '8 GB'
            time = '4h'

            withName: '.*:STAR_ALIGN' {
                cpus = 8
                memory = '64 GB'
                time = '8h'
            }
            withName: '.*:GATK4_SPLITNCIGARREADS' {
                cpus = 2
                memory = '16 GB'
                time = '8h'
            }
            withName: '.*:GATK4_BASERECALIBRATOR' {
                cpus = 2
                memory = '16 GB'
                time = '8h'
            }
            withName: '.*:GATK4_HAPLOTYPECALLER' {
                cpus = 2
                memory = '16 GB'
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
trace    { enabled = true; overwrite = true; file = "${params.outdir}/pipeline_info/execution_trace.txt"     }
dag      { enabled = true; overwrite = true; file = "${params.outdir}/pipeline_info/pipeline_dag.svg"        }
```

rnavar's own `base.config` defines label-based resources only (`process_medium` = 6 CPU/36 GB/8 h, `process_high` = 12 CPU/72 GB/16 h), so the `withName` overrides above use regex selectors (`'.*:NAME'`) that do not depend on the workflow-name prefix. After the first run, compare the selectors with the process names in `{OUTDIR}/pipeline_info/execution_trace.txt` and adjust if any did not match.

The `resourceLimits` values are literal because the pipeline-level maximum-resource parameters of older nf-core templates are not rnavar 1.3.0 parameters and trigger an invalid-parameter schema warning. `overwrite = true` on the four report scopes is needed because a launch that fails early (for example at parameter validation) has already created the files in `pipeline_info/`, and a rerun would otherwise refuse to overwrite them, leaving an empty trace and no HTML reports.

---

## Step 11 — Generate the params file and the submission script

**Params file.** Set `{PARAMS_YAML}` = the samplesheet file name with `_samplesheet.csv` replaced by `_params.yaml`, in the same directory. Write it in `{CWD}` (first check whether `{PARAMS_YAML}` already exists and, if so, ask (numbered): 1. overwrite · 2. choose another filename — before writing) from this template:

```yaml
# nf-core/rnavar {VERSION} parameters — generated by /nfcore-rnavar-setup
input: "{SAMPLESHEET_CSV}"
outdir: "{OUTDIR}"
multiqc_title: "{MULTIQC_TITLE}"
fasta: "{FASTA_PATH}"
gtf: "{GTF_PATH}"
star_index: "{STAR_INDEX}"
read_length: {READ_LENGTH}
seq_platform: "illumina"
{KNOWN_SITES_PARAMS}
{VARIANT_PARAMS}
{ANNOTATION_PARAMS}
```

Typing rules: `seq_platform: "illumina"` is always written (always written, although it equals the default); paths and strings are double-quoted; numbers and booleans are unquoted. `multiqc_title` is always double-quoted so that a numeric-looking title (for example `260928`) stays a string instead of being parsed as a number. Each `{..._PARAMS}` placeholder is replaced by its complete `key: value` lines; if a placeholder is empty, delete that placeholder line entirely.

**Why a params file:** under Nextflow 26.04, values passed as command-line options (a flag followed by a value) reach the nf-schema validator as strings (observed: `--read_length 151` was rejected as "Value is [string] but should be [number]"). A params file keeps the YAML types, so every pipeline parameter goes into it and the launch line carries only `-params-file`.

**Submission script.** Write `nf-core_rnavar_{VERSION}.sh` in `{CWD}` (first check whether `nf-core_rnavar_{VERSION}.sh` already exists and, if so, ask (numbered): 1. overwrite · 2. choose another filename — before writing):

```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p bcc
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={USER_EMAIL}

module add miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
conda activate {CONDA_ENV}
module add singularity/3.10.4

nextflow run nf-core/rnavar -r {VERSION} -c nextflow.config -profile slurm,singularity -params-file {PARAMS_YAML}
```

Show both files in full and instruct:
```
Params file written: {PARAMS_YAML}
Script written: nf-core_rnavar_{VERSION}.sh
To submit:  sbatch nf-core_rnavar_{VERSION}.sh
```
If any helper script (Step 12) was generated, list the order: helpers first, then the pipeline, so the pipeline never starts before its resources exist. If the STAR helper downloads the FASTA, submit the known-sites helper with `sbatch --dependency=afterok:<star_jobid>` (the known-sites guard reads `{FASTA_PATH}`). Submit the pipeline with `sbatch --dependency=afterok:<star_jobid>:<known_sites_jobid> nf-core_rnavar_{VERSION}.sh`, listing every generated helper.

---

## Step 12 — Helper scripts (only for missing resources)

Generate only what is missing. Each is an `sbatch` script (`#SBATCH -N 1 -p bcc --mail-type=END,FAIL` plus resources scaled as below), run on a compute node — never on the login node. Always use `gunzip -c file.gz > file` (never `gunzip -k`; not available on CentOS 7).

**`build_star_index_rnavar_{REF_TAG}.sh`** — downloads (`wget -c`) and decompresses the FASTA and GTF if absent (never downloaded for a custom reference, whose files already exist; a `.gz` custom FASTA/GTF is decompressed as described in Step 6), then computes the genome length and the STAR suffix-array parameter itself, in shell, at run time (the FASTA may not exist yet when the wizard writes the script, so the wizard never substitutes these two values):

```bash
FASTA="{FASTA_PATH}"
GENOME_LENGTH=$(grep -v '^>' "$FASTA" | tr -d '\n' | wc -c)
SA_INDEX_NBASES=$(awk -v L="$GENOME_LENGTH" 'BEGIN{n=int(log(L)/log(2)/2-1); if(n>14)n=14; if(n<4)n=4; print n}')
```

`SA_INDEX_NBASES` = `min(14, floor(log2(GENOME_LENGTH)/2 - 1))`, clamped to a minimum of 4 (40001 bp gives 6; 3.1e9 bp, as for GRCh38, gives 14). STAR's default `--genomeSAindexNbases 14` is far too large for small genomes and makes indexing fail or blow up memory, so it is always passed explicitly. The thread count is taken from SLURM (`${SLURM_NTASKS:-4}`), so it always matches the `-n` of the chosen tier.

**Resource tier — chosen by the wizard when it writes the script** (these respect the HPC defaults of at most 64 G and 4 h). If the FASTA already exists, do not read it with `grep | tr | wc` on the login node (no heavy work there); approximate `{GENOME_LENGTH}` from the file size (`stat -c %s FASTA`) or, if a `.fai` exists, the sum of its length column, and pick the tier; if it is an Ensembl human or mouse download (FASTA not yet on disk), use the over-1-Gb tier. Genome over 1 Gb: `-n 8 --mem=64G -t 4:00:00`; genome under 100 Mb: `-n 4 --mem=8G -t 00:30:00`; in between: `-n 8 --mem=32G -t 2:00:00`. The other helpers need little: `-n 2 --mem=8G -t 4:00:00` (a full dbSNP download plus `bcftools annotate` plus `tabix` may exceed 2 h).

```bash
module add star/2.7.9a
mkdir -p "{GENOME_DIR}/index/star_rnavar_sjdb{SJDB_OVERHANG}"
STAR \
    --runMode genomeGenerate \
    --genomeDir "{GENOME_DIR}/index/star_rnavar_sjdb{SJDB_OVERHANG}" \
    --genomeFastaFiles "$FASTA" \
    --sjdbGTFfile "{GTF_PATH}" \
    --sjdbOverhang {SJDB_OVERHANG} \
    --genomeSAindexNbases "$SA_INDEX_NBASES" \
    --runThreadN "${SLURM_NTASKS:-4}"
```

**`prepare_known_sites_{REF_TAG}.sh`** — must be safely re-runnable. Download each confirmed Step 7 URL into `{GENOME_DIR}/known_sites/`, but skip the download only when `FILE.vcf.gz` exists AND passes `gzip -t` (this works on bgzip files); an interrupted `wget -c` can otherwise leave a truncated `.vcf.gz` that a re-run would accept. Otherwise download to `FILE.vcf.gz.part`, verify it, and only then move it into place:
Then run `bgzip` on any plain `.vcf` before `tabix -f -p vcf`, and only when `FILE.vcf` exists and `FILE.vcf.gz` does not (`bgzip` deletes the `.vcf`, and refuses to overwrite an existing `.gz`). Always re-index with `tabix -f -p vcf` (`tabix` without `-f` errors when an index exists, and the script has no `set -e`, so a stale index would survive). `bcftools annotate -o` overwrites `NEW`. The schema requires `.vcf.gz`.
```bash
if [ -s "$FILE.vcf" ] || { [ -s "$FILE.vcf.gz" ] && gzip -t "$FILE.vcf.gz"; }; then
  :   # already present (a valid .vcf.gz, or the user's own plain .vcf): no download
else
  wget -c -O "$FILE.vcf.gz.part" "URL" \
    && gzip -t "$FILE.vcf.gz.part" && mv "$FILE.vcf.gz.part" "$FILE.vcf.gz" \
    || { echo "ERROR: download or gzip test failed for $FILE.vcf.gz" >&2; exit 1; }
fi
```
This download branch applies only to files that have a confirmed Step 7 URL (a local file without a URL is never downloaded). If the URL is a plain `.vcf` (not `.gz`), use the same skip test and download to `FILE.vcf.part`, check it is non-empty, and move it to `FILE.vcf`; the bgzip rule below then compresses it:
```bash
wget -c -O "$FILE.vcf.part" "URL" && [ -s "$FILE.vcf.part" ] && mv "$FILE.vcf.part" "$FILE.vcf" \
  || { echo "ERROR: download failed for $FILE.vcf" >&2; exit 1; }
```

**Tools come from a container; never load an htslib module.** This cluster has no htslib or tabix module, and under Lmod a failed `module add` (unknown module) silently makes every later `module add` in the same shell fail (observed: bcftools and singularity were then not found). Take `bcftools`, `tabix` and `bgzip` from one Singularity biocontainer, the same mechanism the pipeline uses. The wizard sets `{BCFTOOLS_SIF}` to `${NXF_SINGULARITY_CACHEDIR:-$HOME/.singularity/cache}/depot.galaxyproject.org-singularity-bcftools-1.20--h8b25389_0.img` if that file exists; otherwise the helper downloads the container to that path before use (the SIF download branch is always emitted (the URL is embedded even when the image already exists), so the URL must be re-verified in the session either way), from `https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0`. The wizard must re-verify the URL with a HEAD request (`curl -sI`) or `WebFetch` in the session before writing it into a script (never embed an unverified URL); if the download fails the helper exits 1 (the pipeline then never starts). `--bind` is needed because `/net/...` paths are not auto-bound. Module order matters and each `module add` in a helper should be checked (`|| exit 1`).
```bash
module add singularity/3.10.4 || { echo "ERROR: cannot load singularity" >&2; exit 1; }
command -v singularity >/dev/null || { echo "ERROR: singularity not on PATH" >&2; exit 1; }
SIF="{BCFTOOLS_SIF}"
if [ ! -s "$SIF" ]; then
  mkdir -p "$(dirname "$SIF")"
  wget -c -O "$SIF.part" "https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0" && mv "$SIF.part" "$SIF" \
    || { echo "ERROR: could not download the bcftools container to $SIF" >&2; exit 1; }
fi
[ -s "$SIF" ] || { echo "ERROR: bcftools container not found: $SIF" >&2; exit 1; }
BIND="{GENOME_DIR}"
[ "$(dirname "{FASTA_PATH}")" = "{GENOME_DIR}" ] || BIND="$BIND,$(dirname "{FASTA_PATH}")"
tabix()    { singularity exec --bind "$BIND" "$SIF" tabix "$@"; }
bgzip()    { singularity exec --bind "$BIND" "$SIF" bgzip "$@"; }
bcftools() { singularity exec --bind "$BIND" "$SIF" bcftools "$@"; }
```
Then run an always-run contig guard on every known-sites VCF (this is the enforcement point; it runs whether or not the wizard could check earlier). For each `FILE`, compare contig names and either continue, rename, or fail fast:
```bash
VCF_CONTIG=$(tabix -l FILE | head -n1)
FASTA_CONTIG=$(grep -m1 '^>' {FASTA_PATH} | cut -d' ' -f1 | sed 's/^>//')
if [ -z "$VCF_CONTIG" ] || [ -z "$FASTA_CONTIG" ]; then
  echo "ERROR: empty contig (VCF: '$VCF_CONTIG', FASTA: '$FASTA_CONTIG') for FILE; the pipeline must not be run" >&2
  exit 1
fi
if [ "$VCF_CONTIG" != "$FASTA_CONTIG" ]; then
  NEW="${FILE%.vcf.gz}.renamed.vcf.gz"
  # MAP.txt (in {GENOME_DIR}/known_sites/) is generated for the needed direction, see below
  bcftools annotate --rename-chrs MAP.txt -O z -o "$NEW" FILE
  tabix -f -p vcf "$NEW"          # then use "$NEW" in place of FILE
  VCF_CONTIG=$(tabix -l "$NEW" | head -n1)
  if [ -z "$VCF_CONTIG" ] || [ "$VCF_CONTIG" != "$FASTA_CONTIG" ]; then
    echo "ERROR: contig mismatch after rename: FILE has '$VCF_CONTIG', FASTA has '$FASTA_CONTIG'; the pipeline must not be run" >&2
    exit 1
  fi
  grep -qF "$NEW" "{CWD}/{PARAMS_YAML}" || { echo "ERROR: renamed to $NEW but {PARAMS_YAML} points elsewhere; update dbsnp/known_indels (+_tbi) and resubmit; the pipeline must not be run" >&2; exit 1; }
fi
```
MAP.txt is a two-column, tab-separated old-name/new-name file written by the helper for the direction that is needed (chr1<->1 ... chrM<->MT for the chromosomes present in the FASTA). If the VCF uses a `chr` prefix and the FASTA does not, it maps `chrN` to `N` for N = 1-22, X, Y and `chrM` to `MT`; if the FASTA uses a `chr` prefix and the VCF does not, it maps the reverse. The order is rename, then bgzip (`-O z`), then `tabix -f -p vcf`. The error message must name the VCF, both contig names and say that the pipeline must not be run; because the pipeline is submitted with `--dependency=afterok`, a failed guard stops it from starting. The guard fails with exit 1 if either contig is empty (an empty VCF, or a FASTA without a header), because an empty string would otherwise compare equal to itself or slip through. Read the VCF contig with `tabix -l FILE | head -n1`, not `zcat | grep | head -1`, because the contig list is tiny and the helper does not set `pipefail` (a broken pipe would otherwise matter). `BIND` lists each directory once, because Singularity prints "destination is already in the mount point list" on every call when the FASTA directory equals `{GENOME_DIR}`.

**Post-guard filenames.** These are the post-guard filenames. `NEW` is always `${FILE%.vcf.gz}.renamed.vcf.gz`. The `dbsnp`/`known_indels` keys (and their `.tbi` companions) in `{PARAMS_YAML}` must point at the files that pass the guard: `NEW` (and `NEW.tbi`) whenever a rename happens, otherwise the original `FILE`. The wizard therefore decides the final filenames before writing the params file (Step 11), writing the `.renamed.vcf.gz` names whenever a rename is expected (for example hg38 resource-bundle `chr`-prefixed VCFs against an Ensembl FASTA). If a rename was not anticipated, the guard's `grep` check exits 1 rather than letting the pipeline run on the original file.

**Final params-file check (end of the helper, after the per-file loop).** Whether or not a rename happened, verify that every file named in the params file exists and is non-empty, so a `.renamed.vcf.gz` name written up front without a rename cannot reach the pipeline:
```bash
for K in dbsnp dbsnp_tbi known_indels known_indels_tbi; do
  P=$(grep "^$K:" "{CWD}/{PARAMS_YAML}" | sed 's/^[^:]*: *"\(.*\)" *$/\1/')
  [ -s "$P" ] || { echo "ERROR: $K in {PARAMS_YAML} points to a missing or empty file: '$P'; the pipeline must not be run" >&2; exit 1; }
done
```

**`prepare_annotation_cache_{TOOL}.sh`** — only if the user has no cache: a script the user runs where internet is available, using the tool's own cache installer (`vep_install` or `snpEff download`) into the directory given as `vep_cache` / `snpeff_cache`, where `{TOOL}` is `snpeff` or `vep`, matching `{ANNOTATION_TOOL}` (one script per chosen tool).

---

## Step 13 — Hand-off note

Print where the results will be, so later analyses can find them:
```
Outputs under {OUTDIR}/ :
  variant_calling/   per sample: {SAMPLE}.haplotypecaller.filtered.vcf.gz (soft-filtered), raw {SAMPLE}.haplotypecaller.vcf.gz (each with .tbi); with `generate_gvcf: true`, gVCFs too
  preprocessing/     per sample: {SAMPLE}.md.bam (duplicate-marked) and {SAMPLE}.recal.bam (with .bai); with `skip_baserecalibration: true` these are the duplicate-marked BAMs and no recal.bam is written
  annotation/        SnpEff / VEP output (only when `tools` is set)
  reports/           FastQC, samtools, Picard and STAR reports; MultiQC at reports/multiqc/{MULTIQC_TITLE}_multiqc_report.html
  pipeline_info/     execution report, timeline, trace, DAG
These VCFs and BAMs are the inputs expected by the ase-pipeline skill (allele-specific expression).
```
Before printing, confirm the actual directory names against the pipeline's `docs/output` page for `{VERSION}` with `WebFetch`, and use the real names.

---

## Notes for the assistant

- **`gh` CLI is not available on this HPC cluster.** Use `WebFetch` for all GitHub API calls.
- **Always present finite-choice questions as numbered lists.** Use open questions only when no reasonable discrete set exists (email, conda env, custom paths).
- Never run heavy computation on the login node; all work through `sbatch`.
- Raw FASTQ/BAM files are read-only; never modify them.
- Omit params-file keys whose value equals the pipeline default.
- Never pass rnavar parameters on the `nextflow run` command line — always the params file.
- Never overwrite an existing `nextflow.config`.
- `read_length` is always written to the params file; a STAR index made for a different read length is never reused.
- Known sites are never assumed — either supply all four files or `skip_baserecalibration: true`.
- Never embed a download URL that was not verified in this session.
- Not rnavar parameters, never emit: `annotation_cache`, `gencode`, `strandedness`.
