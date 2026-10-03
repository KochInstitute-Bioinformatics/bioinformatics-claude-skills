#!/bin/bash
# Static checks for nfcore-rnavar-setup.md. Usage: check_skill.sh <skill.md> [schema.json]
set -u
SKILL=${1:?usage: check_skill.sh <skill.md> [schema.json]}
SCHEMA=${2:-/tmp/rnavar_schema_1.3.0.json}
fail=0

if [ ! -s "$SKILL" ]; then echo "FAIL: skill file missing or empty: $SKILL"; exit 1; fi
if [ ! -s "$SCHEMA" ]; then
  curl -sfL https://raw.githubusercontent.com/nf-core/rnavar/1.3.0/nextflow_schema.json -o "$SCHEMA" \
    || { echo "cannot fetch schema"; exit 2; }
fi

need()   { grep -qF -- "$1" "$SKILL" || { echo "FAIL: missing required text: $1"; fail=1; }; }
forbid() { ! grep -qF -- "$1" "$SKILL" || { echo "FAIL: forbidden text present: $1"; fail=1; }; }

# (a) every --flag in the skill is a schema parameter or an allowlisted non-rnavar flag
schema_names=$(grep -oE '"[a-z_0-9]+": *\{' "$SCHEMA" | sed -E 's/"([a-z_0-9]+)".*/\1/' | sort -u)
allow="bind genomeSAindexNbases rename-chrs dependency runMode genomeDir genomeFastaFiles sjdbGTFfile sjdbOverhang runThreadN outTmpDir mail-type mail-user mem"
for flag in $(grep -oE '(^|[ `(=])--[A-Za-z_][A-Za-z_0-9-]*' "$SKILL" | sed -E 's/^[^-]*--//' | sort -u); do
  if ! echo "$schema_names $allow" | tr ' ' '\n' | grep -qx -- "$flag"; then
    echo "FAIL: flag not in rnavar schema or allowlist: --$flag"; fail=1
  fi
done

# (b)/(c) content checks, appended per task -------------------------------
# --- Task 1
need "# nf-core/rnavar Pipeline Setup Skill"
need "## Step 0"
need "## Step 3"
need "api.github.com/repos/nf-core/rnavar/releases/latest"
need "Do NOT use the \`gh\` CLI"
forbid "--annotation_cache"
forbid "--gencode"
forbid "--strandedness"
# --- end Task 1

# --- Task 2
need "## Step 4"
need "sample,fastq_1,fastq_2"
need "sample,bam,bai"
need "sample,cram,crai"
need "Supplying FASTQ files and a BAM/CRAM file for the same sample"
need "Replace every \`-\` with \`_\`"
need "merged before alignment"
need "## Step 5"
need "most common read length"
need "→ sjdbOverhang {READ_LENGTH − 1}"
forbid "strandedness,"
# --- end Task 2

# --- Task 3
need "## Step 6"
need "star_rnavar_sjdb"
need "GTF source"
need "no flag is emitted"
need "## Step 7"
need "dbsnp:"
need "dbsnp_tbi:"
need "known_indels:"
need "known_indels_tbi:"
need "skip_baserecalibration: true"
need "does not skip base recalibration automatically"
need "resolve the resource URLs at run time"
need 'Store: `{GENOME_DIR}` ='
need "use the highest N unless the user asks otherwise"
need "still apply the STAR index rule below"
need "{STAR_INDEX} is the same path"
# --- end Task 3

# --- Task 4
need "## Step 8"
need "## Step 9"
need "needs internet from compute nodes"
need "not a parameter in the rnavar schema"
# --- end Task 4

# --- Task 5
need "## Step 10"
need "do not overwrite it"
need "resourceLimits = [ cpus: 16"
need "'.*:STAR_ALIGN'"
need "'.*:GATK4_HAPLOTYPECALLER'"
need "'.*:GATK4_BASERECALIBRATOR'"
need "'.*:GATK4_SPLITNCIGARREADS'"
need "## Step 11"
need "nextflow run nf-core/rnavar -r {VERSION}"
need 'seq_platform: "illumina"'
need "read_length: {READ_LENGTH}"
need 'star_index: "{STAR_INDEX}"'
need "## Step 12"
need "build_star_index_rnavar"
need "gunzip -c file.gz > file"
need "## Step 13"
need "## Notes for the assistant"
need "Raw FASTQ/BAM files are read-only"
need 'as `{KNOWN_SITES_PARAMS}`'
need 'as `{VARIANT_PARAMS}`'
need 'as `{ANNOTATION_PARAMS}`'
need 'numbers and booleans are unquoted'
need 'so that a numeric-looking title'
need 'and `{ENS_VERSION}`'
need "where \`{TOOL}\` is"
forbid "{KNOWN_SITES_LINES}"
forbid "{VARIANT_LINES}"
forbid "{ANNOTATION_LINES}"
forbid "{VERSION_ENS}"
# --- end Task 5


# --- Final-review fixes
need "bgzipped \`.vcf.gz\` files with \`.tbi\` indexes"
need "contig-name check"
forbid "zcat FILE.vcf.gz | grep -v '^#' | head -1 | cut -f1"
need "prefer Ensembl-named variation VCFs"
need "bcftools annotate --rename-chrs MAP.txt"
need "chr1<->1 ... chrM<->MT"
need "never proceed with mismatched contigs"
need "run \`bgzip\` on any plain \`.vcf\` before \`tabix -f -p vcf\`"
need "1. these are lanes of the same sample"
need "2. these are different samples"
need "allow the deliberate duplicates"
need "Write \`download_cache: true\` only if the user explicitly chooses it after being warned."
forbid "Use \`snpeff_cache\`, \`vep_cache\` and \`download_cache\` only."
need "--dependency=afterok:<star_jobid>:<known_sites_jobid> nf-core_rnavar_{VERSION}.sh"
need "with \`skip_baserecalibration: true\` these are the duplicate-marked BAMs"
need "check whether \`{SAMPLESHEET_CSV}\` already exists"
need "check whether \`nf-core_rnavar_{VERSION}.sh\` already exists"
need "1. overwrite · 2. choose another filename"
# --- end Final-review fixes

# --- Contig-guard fix
need "always-run contig guard"
need 'VCF_CONTIG=$(tabix -l'
need 'FASTA_CONTIG=$(grep -m1'
need 'exit 1'
need "the pipeline must not be run"
need "chrM\` to \`MT\`"
need "helper's contig guard performs the check after download"
need "Mouse Genomes Project VCFs use Ensembl-style contig names"
forbid "run the check once they are"
forbid "If the Step 7 contig-name check found a mismatch and option (ii) was chosen"
# --- end Contig-guard fix

# --- Params-file rules
# (1) The launch command must carry no rnavar --flag (params file only)
launch=$(awk '/^nextflow run nf-core\/rnavar/{p=1} p{l=l $0; if($0 !~ /\\[ \t]*$/){print l; exit}}' "$SKILL" | sed 's/\\[ \t]*/ /g')
if [ -z "$launch" ]; then echo "FAIL: no 'nextflow run nf-core/rnavar' launch line found"; fail=1
elif echo "$launch" | grep -qE ' --[A-Za-z_]'; then echo "FAIL: launch line carries a --flag (must use -params-file only): $launch"; fail=1
fi
need "-params-file {PARAMS_YAML}"

# (2) Every key in the params-file template (fenced yaml block) must be a schema parameter
yaml_keys=$(awk '/^```yaml/{f=1;next} /^```/{f=0} f' "$SKILL" | grep -oE '^[a-z_0-9]+:' | tr -d ':' | sort -u)
[ -n "$yaml_keys" ] || { echo "FAIL: no fenced yaml params template found"; fail=1; }
for k in $yaml_keys; do
  echo "$schema_names" | grep -qx -- "$k" || { echo "FAIL: params-file key not in rnavar schema: $k"; fail=1; }
done

# --- Task 2 (custom reference + helpers)
need "Custom reference"
need "{REF_TAG}"
need "custom_{WD_NAME}"
need "Other organisms: use option 2"
need "build_star_index_rnavar_{REF_TAG}.sh"
need "prepare_known_sites_{REF_TAG}.sh"
need "prepare_annotation_cache_{TOOL}.sh"
need "{GENOME_LENGTH}"
need "--genomeSAindexNbases 14"
need 'SA_INDEX_NBASES=$(awk'
need '--genomeSAindexNbases "$SA_INDEX_NBASES"'
need '${SLURM_NTASKS:-4}'
forbid "{SA_INDEX_NBASES}"
forbid "{N_THREADS}"
need "resources scaled as below"
need "tabix -l FILE | head -n1"
need "if either contig is empty"
need "**Post-guard filenames.**"
need "the wizard also generates \`prepare_known_sites_{REF_TAG}.sh\`"
forbid "build_star_index_rnavar_{ASSEMBLY}_ens{ENS_VERSION}.sh"
forbid "prepare_known_sites_{ASSEMBLY}"
forbid "-t 8:00:00"
# --- end Task 2 (custom reference + helpers)


# --- Task 3 (config, hand-off)
need "timeline { enabled = true; overwrite = true;"
[ "$(grep -cE '^(timeline|report|trace|dag) +\{ enabled = true; overwrite = true;' "$SKILL")" -eq 4 ] || { echo "FAIL: overwrite = true must appear for timeline, report, trace and dag"; fail=1; }
forbid "params.max_"
forbid "max_cpus"
need "reports/multiqc"
forbid "  multiqc/           MultiQC report"
need "  annotation/        SnpEff"
# --- end Task 3 (config, hand-off)
# --- end Params-file rules


# --- Final fix wave
forbid "module add htslib"
need "singularity exec --bind"
need 'SIF="{BCFTOOLS_SIF}"'
need "depot.galaxyproject.org-singularity-bcftools-1.20--h8b25389_0.img"
need "each \`module add\` in a helper"
need "tabix -f -p vcf"
forbid "tabix -p vcf"
need "safely re-runnable"
need "skip the download only when \`FILE.vcf.gz\` exists AND passes \`gzip -t\`"
need "a custom FASTA may use either style"
need "Treat an empty value as a mismatch."
need "EMPTY contig"
forbid "The FASTA from Step 6 is Ensembl-named (\`1\`, \`2\`, ... \`MT\`), whereas"
need "sbatch --dependency=afterok:<star_jobid>"
need "afterok:<star_jobid>:<known_sites_jobid>"
forbid "afterok:<helper_jobid>"
need ".renamed.vcf.gz"
need "points elsewhere"
need "gunzip -c file.gz > {GENOME_DIR}/<name>"
need "originals are never modified"
need "stat -c %s"
forbid "measure it with the same"
need "-n 2 --mem=8G -t 4:00:00"
forbid "-n 2 --mem=8G -t 2:00:00"
# --- end Final fix wave

# --- Parked-items fixes
# A: bcftools container download branch before the existence check
need 'if [ ! -s "$SIF" ]; then'
need 'wget -c -nv -O "$SIF.part" "https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0" && mv "$SIF.part" "$SIF"'
need 'ERROR: could not download the bcftools container to $SIF'
need "re-verify the URL with a HEAD request (\`curl -sI\`) or \`WebFetch\`"
need "the helper exits 1 (the pipeline then never starts)"
dl_ln=$(grep -nF 'wget -c -nv -O "$SIF.part"' "$SKILL" | head -1 | cut -d: -f1)
ck_ln=$(grep -nF 'ERROR: bcftools container not found' "$SKILL" | head -1 | cut -d: -f1)
{ [ -n "$dl_ln" ] && [ -n "$ck_ln" ] && [ "$dl_ln" -lt "$ck_ln" ]; } || { echo "FAIL: container download branch must precede the SIF existence check"; fail=1; }
# B: FASTA_SOURCE for the login-node contig check
need 'set `{FASTA_SOURCE}`'
need 'FASTA_CONTIG=$(zcat -f {FASTA_SOURCE} | awk'
forbid 'zcat -f {FASTA_PATH}'
need "the original FASTA path the user gave, gzipped or not"
# C: params-file existence check at the end of the helper
need 'for K in dbsnp dbsnp_tbi known_indels known_indels_tbi; do'
need 'points to a missing or empty file'
# D: BIND without a repeated directory
need 'BIND="{GENOME_DIR}"'
need '[ "$(dirname "{FASTA_PATH}")" = "{GENOME_DIR}" ] || BIND="$BIND,$(dirname "{FASTA_PATH}")"'
forbid 'BIND="{GENOME_DIR},$(dirname "{FASTA_PATH}")"'
# E: truncated-download safety
need 'gzip -t "$FILE.vcf.gz"'
need 'wget -c -nv -O "$FILE.vcf.gz.part" "URL"'
need 'gzip -t "$FILE.vcf.gz.part" && mv "$FILE.vcf.gz.part" "$FILE.vcf.gz"'
forbid "Skip \`wget\` when \`FILE.vcf.gz\` already exists"
# F1-F4
need "(always written, although it equals the default)"
need "\`{DBSNP}\` and \`{INDELS}\` are the FINAL post-guard paths"
forbid '2. "What is the base directory where genome files and indexes are stored?"'
need "Ask the base directory only for option 1"
need "because the contig list is tiny and the helper does not set \`pipefail\`"
forbid "which can die of SIGPIPE under"
# F5a: every backticked key: token in Steps 7-9 is a schema parameter
for k in $(awk '/^## Step 7/{f=1} /^## Step 10/{f=0} f' "$SKILL" | grep -oE '`[a-z_0-9]+:' | tr -d '`:' | sort -u); do
  echo "$schema_names" | grep -qx -- "$k" || { echo "FAIL: Step 7-9 backticked key not in rnavar schema: $k"; fail=1; }
done
[ -n "$(awk '/^## Step 7/{f=1} /^## Step 10/{f=0} f' "$SKILL" | grep -oE '`[a-z_0-9]+:')" ] || { echo "FAIL: no backticked keys found in Steps 7-9"; fail=1; }
# F5b: tightened key-form needs (replace the former bare-word needs)
need 'read_length: {READ_LENGTH}`'
need 'remove_duplicates: true'
need 'star_twopass: false'
need 'gatk_hc_call_conf: <int>'
need 'gatk_vf_qd_filter: <number>'
need 'gatk_vf_fs_filter: <number>'
need 'gatk_vf_window_size: <int>'
need 'gatk_vf_cluster_size: <int>'
need 'skip_variantfiltration: true'
need 'generate_gvcf: true'
need 'bam_csi_index: true'
need 'tools: "snpeff"'
need 'tools: "merge"'
need 'snpeff_cache: "{DIR}"'
need 'vep_cache: "{DIR}"'
need 'snpeff_db: "{SNPEFF_DB}"'
need 'vep_genome: "{VEP_GENOME}"'
need 'vep_species: "{VEP_SPECIES}"'
need 'vep_cache_version: "{VEP_CACHE_VERSION}"'
need 'download_cache: true'
need 'The `star_index` key is always written'
# --- end Parked-items fixes

# --- Parked-items fix round 1
need 'if [ -s "$FILE.vcf" ] || { [ -s "$FILE.vcf.gz" ] && gzip -t "$FILE.vcf.gz"; }; then'
need 'wget -c -nv -O "$FILE.vcf.part" "URL" && [ -s "$FILE.vcf.part" ] && mv "$FILE.vcf.part" "$FILE.vcf"'
need "If the URL is a plain \`.vcf\` (not \`.gz\`)"
need "applies only to files that have a confirmed Step 7 URL"
need "Step 6 option 1 binds \`{FASTA_SOURCE}\` = \`{FASTA_PATH}\`"
need "skip the wizard-side check and rely on the helper's contig guard"
need "always emitted (the URL is embedded even when the image already exists)"
# --- end Parked-items fix round 1

# --- Real-data fix round (GM12878 / HG001 test, 2026-10-02)
README="$(dirname "$SKILL")/README.md"
FIXDIR="$(cd "$(dirname "$0")" && pwd)/fixtures"
needr()   { grep -qF -- "$1" "$README" 2>/dev/null || { echo "FAIL: README missing required text: $1"; fail=1; }; }
forbidr() { ! grep -qF -- "$1" "$README" 2>/dev/null || { echo "FAIL: README has forbidden text: $1"; fail=1; }; }
# D1: a working resolution path for the GATK hg38 bundle
need "https://gatk.broadinstitute.org/hc/en-us/articles/360035890811"
need "**Human: when the GATK page cannot be read.**"
need 'https://storage.googleapis.com/storage/v1/b/gcp-public-data--broad-references/o?prefix=hg38/v0/&fields=items(name,size),nextPageToken'
need '&pageToken=<nextPageToken>'
need "The download URL of an object is https://storage.googleapis.com/gcp-public-data--broad-references/ followed by its listed name"
need "show the exact URLs and sizes to the user and download only after they confirm"
# D2: which bundle files
need '`Homo_sapiens_assembly38.dbsnp138.vcf.gz` (dbSNP, 1,560,889,937 bytes)'
need '`Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` (known indels, 20,685,880 bytes)'
need 'Do NOT use the plain `Homo_sapiens_assembly38.dbsnp138.vcf`'
# D8: head-job resources
need '#SBATCH -N 1 -n 2 --mem=8G -t 2-00:00:00 -p bcc'
forbid '#SBATCH -n 32'
need 'a deliberate exception to the usual 4 h default'
sub_block=$(awk '/^```bash/{f=1;b="";next} /^```/{if(f&&b~/nextflow run nf-core\/rnavar/)print b; f=0;next} f{b=b $0 "\n"}' "$SKILL")
{ echo "$sub_block" | grep -qE '^#SBATCH .*-t [0-9]' && echo "$sub_block" | grep -qE '^#SBATCH .*--mem=' \
  && echo "$sub_block" | grep -qE '^#SBATCH .*-n 2( |$)' && ! echo "$sub_block" | grep -qE '^#SBATCH .*-n ([013-9]|[0-9]{2,})( |$)'; } \
  || { echo "FAIL: the pipeline submission script must request a time (-t), memory (--mem=) and exactly -n 2"; fail=1; }
# D3, D4
need '`_1.` (SRA/ENA names: `SRR5665260_1.fastq.gz` gives `SRR5665260`)'
need 'When `{SEQ_DATE}` equals `{TODAY_YYMMDD}` (no date prefix was found, so it fell back to today), options 1 and 2 are the same: offer only 1. `{TODAY_YYMMDD}_{WD_NAME}_samplesheet.csv` · 2. Custom.'
need 'For both questions, when `{SEQ_DATE}` equals `{TODAY_YYMMDD}`, options 1 and 2 are the same: offer only 1. `{TODAY_YYMMDD}_{WD_NAME}` · 2. Custom.'
# D5: every module line is checked
need 'module add star/2.7.9a || { echo "ERROR: cannot load star" >&2; exit 1; }'
need 'module add miniconda3/v4 || { echo "ERROR: cannot load miniconda3" >&2; exit 1; }'
bad_mod=$(grep -nE '^[[:space:]]*module add ' "$SKILL" | grep -vE '\|\| (exit 1|\{ echo "ERROR: [^"]*" >&2; exit 1; \})$')
[ -z "$bad_mod" ] || { echo "FAIL: unchecked module add line(s) (only '|| exit 1' or '|| { echo \"ERROR: ...\" >&2; exit 1; }' passes): $bad_mod"; fail=1; }
# D6: quiet downloads
need 'downloads (`wget -c -nv`'
bad_wget=$(grep -nE 'wget -[A-Za-z]+ ' "$SKILL" | grep -vF -- '-nv')
[ -z "$bad_wget" ] || { echo "FAIL: wget line(s) without -nv: $bad_wget"; fail=1; }
# D9: STAR temp dir next to the index, removed by literal path only
need '--outTmpDir "{GENOME_DIR}/index/_STARtmp_rnavar_sjdb{SJDB_OVERHANG}"'
[ "$(grep -cxF 'rm -rf -- "{GENOME_DIR}/index/_STARtmp_rnavar_sjdb{SJDB_OVERHANG}"' "$SKILL")" -eq 2 ] \
  || { echo "FAIL: the STAR helper must remove its temp dir by literal path before and after STAR (2 lines)"; fail=1; }
need 'STAR_RC=$?'
need 'ERROR: STAR genomeGenerate failed (exit $STAR_RC)'
bad_rm=$(grep -nE '(^|[^A-Za-z])rm (.* )?(-[A-Za-z]*[rR][A-Za-z]*|--recursive)( .*)?\$' "$SKILL")
[ -z "$bad_rm" ] || { echo "FAIL: recursive rm on a shell variable: $bad_rm"; fail=1; }
# D7: generic defaults kept, described as a fallback
need "// Fallback for unlabelled processes only"
awk '/Fallback for unlabelled processes only/{f=1} f&&/withName/{exit} f&&/^ +cpus = 2$/{c=1} f&&/^ +memory = .8 GB.$/{m=1} f&&/^ +time = .4h.$/{t=1} END{exit !(c&&m&&t)}' "$SKILL" \
  || { echo "FAIL: the fallback cpus = 2 / memory = '8 GB' / time = '4h' must follow the fallback comment"; fail=1; }
# Tiers: exact values per selector (measured on one dataset)
sel_vals() { awk -v s="withName: '$1' {" 'index($0,s){f=1;next} f&&/^ *}$/{exit} f{gsub(/^ +/,""); printf "%s;", $0}' "$SKILL"; }
chk_sel() {
  [ "$(grep -cF "withName: '$1' {" "$SKILL")" -eq 1 ] || { echo "FAIL: selector $1 must appear exactly once (Nextflow applies the last block)"; fail=1; }
  [ "$(sel_vals "$1")" = "$2" ] || { echo "FAIL: selector $1 values are '$(sel_vals "$1")', expected '$2'"; fail=1; }; }
# values for attempt 1; memory and time scale with task.attempt (rnavar's retry), cpus fixed
chk_sel '.*:STAR_ALIGN'             "cpus = 8;memory = { 64.GB * task.attempt };time = { 8.h * task.attempt };"
chk_sel '.*:PICARD_MARKDUPLICATES'  "cpus = 2;memory = { 48.GB * task.attempt };time = { 8.h * task.attempt };"
chk_sel '.*:GATK4_SPLITNCIGARREADS' "cpus = 4;memory = { 24.GB * task.attempt };time = { 4.h * task.attempt };"
chk_sel '.*:GATK4_BASERECALIBRATOR' "cpus = 2;memory = { 8.GB * task.attempt };time = { 4.h * task.attempt };"
chk_sel '.*:GATK4_HAPLOTYPECALLER'  "cpus = 2;memory = { 8.GB * task.attempt };time = { 4.h * task.attempt };"
# every memory/time line inside any withName block is an attempt-scaled closure (also for selectors added later)
bad_sel=$(awk '/withName: /{f=1;next} f&&/^ *}$/{f=0} f&&/(memory|time) =/&&!/^ +(memory = \{ [0-9]+\.GB|time = \{ [0-9]+\.h) \* task\.attempt \}$/{print NR": "$0}' "$SKILL")
[ -z "$bad_sel" ] || { echo "FAIL: withName memory/time must be '{ N.GB * task.attempt }' / '{ N.h * task.attempt }': $bad_sel"; fail=1; }
need "resourceLimits = [ cpus: 16, memory: '64 GB', time: '24h' ]"
need "// Values are for attempt 1; a retried task (time or memory kill) gets them x2, capped by resourceLimits."
need "the values above are for attempt 1 and are doubled on the retry; \`cpus\` stays fixed."
need "At attempt 2: STAR_ALIGN asks for 128 GB / 16 h and gets 64 GB / 16 h (the same memory, so a memory kill of STAR_ALIGN is not helped by the retry)"
need "If a task still fails after the retry with exit 140/143 (time) or 137 (memory), raise that selector"
need "Tiers measured on one human sample (39.6 M read pairs)"
need "MarkDuplicates, which had no selector, peaked at 27.8 GB of its label's 36 GB (77%); the 48 GB is a headroom judgment for deeper libraries, not a measured need"
forbid "hence 48 GB of headroom for deeper libraries"
need 'Genome over 1 Gb: `-n 8 --mem=64G -t 4:00:00` (measured for GRCh38 in the real-data test: 42 min 27 s and at least 43.1 GB'
# every withName selector matches a process of the real run (fixture from its execution_trace.txt)
PROCS="$FIXDIR/realdata_gm12878_processes.txt"
if [ ! -s "$PROCS" ]; then echo "FAIL: fixture missing: $PROCS"; fail=1
else
  sels=$(grep -oE "withName: '[^']+'" "$SKILL" | sed -E "s/withName: '([^']+)'/\1/")
  [ -n "$sels" ] || { echo "FAIL: no withName selectors found"; fail=1; }
  set -f
  for s in $sels; do
    grep -v '^#' "$PROCS" | grep -qxE -- "$s" || { echo "FAIL: selector '$s' matches no process of the real-data run"; fail=1; }
  done
  set +f
fi
# Hand-off note: RNA editing and evaluation limits (measured on one dataset)
need "measured on ONE human dataset: GM12878/HG001"
need "78% of the false-positive SNVs were A>G/T>C (2,300 of 2,932), and 88.5% of those"
need "(2,036) sit on known REDIportal editing sites, against 1.37% of the true-positive SNVs."
need "Mask known editing sites (for example REDIportal) before treating A>G/T>C calls as genomic"
need "would rise from 0.885 to about 0.96, a counterfactual that was not re-run, at a cost of"
need "about 1.3 points of SNV recall (308 true SNVs also sit on REDIportal positions)."
# bucket listing: status check, no silent fall-through
need "OUT=\$(curl -s -w '\\nHTTP_STATUS=%{http_code}' \"https://storage.googleapis.com/storage/v1/b/gcp-public-data--broad-references/o?prefix=hg38/v0/&fields=items(name,size),nextPageToken\")"
need 'if [ "$CODE" != 200 ] || [ -z "${BODY//[[:space:]]/}" ]; then'
need 'echo "ERROR: bucket listing failed (HTTP $CODE, ${#BODY} bytes): stop and tell the user; never guess the URLs"'
need "On the ERROR line, stop: do not fall back to paging or to URLs from memory"
forbid 'curl -s "https://storage.googleapis.com/storage/v1/b/gcp-public-data--broad-references/o?prefix=hg38/v0/&fields=items(name,size),nextPageToken" | grep'
need "76.5% of the missed truth variants were heterozygous"
need "SnpCluster (91% of the filtered ones)"
need "32.9 Mb, about 1.3% of the GIAB confident genome. No accuracy claim is made outside it;"
# README: real-data test record and honesty statements
needr "## Real-data test (one human dataset)"
needr "It is ONE dataset, one sample and one run; the numbers below hold for that dataset only."
needr "SRR5665260 (BioProject PRJNA389940, study accession from the ENA filereport API for SRR5665260, fetched 2026-10-02), GM12878 = GIAB HG001, a whole human RNA-seq sample"
forbidr "whole-transcriptome"
while IFS= read -r row; do needr "$row"; done <<'ROWS'
| exons_all_dp10 | allele | 0.891 | 0.923 | 0.907 | 0.733 | 0.886 | 0.802 |
| exons_all_dp20 | genotype | 0.924 | 0.922 | 0.923 | 0.743 | 0.856 | 0.795 |
| exons_all_dp20 | allele | 0.927 | 0.925 | 0.926 | 0.800 | 0.922 | 0.857 |
| exons_pc_dp10 | genotype | 0.905 | 0.935 | 0.919 | 0.675 | 0.821 | 0.741 |
| exons_pc_dp10 | allele | 0.909 | 0.940 | 0.924 | 0.738 | 0.897 | 0.810 |
| exons_pc_dp20 | genotype | 0.937 | 0.940 | 0.939 | 0.738 | 0.866 | 0.797 |
| exons_pc_dp20 | allele | 0.940 | 0.943 | 0.942 | 0.796 | 0.935 | 0.860 |
| PASS | 0.885 | 0.917 | 0.901 | 0.674 | 0.815 | 0.738 |
ROWS
needr "consistent with RNA editing (88.5% of those sit on known editing sites)"
forbidr "they are RNA editing"
needr "at a cost of about 1.3 points of SNV recall (308 true SNVs also sit on REDIportal positions); both are estimates"
needr "for attempt 1: memory and time double on the automatic retry"
needr "Memory and time double on the one automatic retry, up to 64 GB and 24 h"
needr "The 48 GB for MarkDuplicates is a headroom judgment, not a measured need"
needr "### Changes made after the real-data test"
forbidr "fix round"
forbidr "Defects found and fixed in this round"
forbidr "this round"
d_ids=$(grep -nE '(^|[^A-Za-z0-9])D[1-9]([^0-9]|$)' "$README" 2>/dev/null)
[ -z "$d_ids" ] || { echo "FAIL: README uses internal defect IDs (D1-D9): $d_ids"; fail=1; }
needr "Pipeline wall time **3 h 10 min 42 s**; **22.9 CPU-hours** used by all tasks"
needr "| exons_all_dp10 | genotype | 0.885 | 0.917 | 0.901 | 0.674 | 0.815 | 0.738 |"
needr "| all records | 0.766 | 0.957 | 0.851 | 0.534 | 0.823 | 0.648 |"
needr "85,799 calls, 84.6% PASS; ts/tv 2.92 (PASS SNVs, autosomes); het/hom-alt 1.15; 72.9% of PASS calls in dbSNP138"
needr "No accuracy claim is made outside GIAB confident regions ∩ exons ∩ depth ≥ 10 (or ≥ 20)"
needr "about 67 GB under the test folder"
needr "STAR index 30 GB and \`known_sites/\` 2.9 GB"
needr "### Resources: requested vs observed (that run)"
needr "| Head job (Nextflow) | 32 cores, no memory or time | about 1 core |"
needr "### Not verified"
needr "A cluster that enforces memory limits."
needr "A second run for reproducibility (one run only)."
needr "have not been run on real data; they were checked statically and with stubs."
needr "The resource tiers are measured on one human dataset"
needr "deliberate exception to the usual 4 h default"
forbidr "Process-name selectors were verified against real tasks only for the four listed"
forbidr "validated on real data"
forbidr "validated on real genomes"
forbidr "fully validated"
forbidr "production-ready"
forbidr "validated genome-wide"
# --- end Real-data fix round

[ $fail -eq 0 ] && echo "PASS" || exit 1
