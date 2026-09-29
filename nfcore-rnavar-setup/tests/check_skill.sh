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
allow="bind genomeSAindexNbases rename-chrs dependency runMode genomeDir genomeFastaFiles sjdbGTFfile sjdbOverhang runThreadN mail-type mail-user mem"
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
need 'wget -c -O "$SIF.part" "https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0" && mv "$SIF.part" "$SIF"'
need 'ERROR: could not download the bcftools container to $SIF'
need "re-verify the URL with a HEAD request (\`curl -sI\`) or \`WebFetch\`"
need "the helper exits 1 (the pipeline then never starts)"
dl_ln=$(grep -nF 'wget -c -O "$SIF.part"' "$SKILL" | head -1 | cut -d: -f1)
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
need 'wget -c -O "$FILE.vcf.gz.part" "URL"'
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
need 'wget -c -O "$FILE.vcf.part" "URL" && [ -s "$FILE.vcf.part" ] && mv "$FILE.vcf.part" "$FILE.vcf"'
need "If the URL is a plain \`.vcf\` (not \`.gz\`)"
need "applies only to files that have a confirmed Step 7 URL"
need "Step 6 option 1 binds \`{FASTA_SOURCE}\` = \`{FASTA_PATH}\`"
need "skip the wizard-side check and rely on the helper's contig guard"
need "always emitted (the URL is embedded even when the image already exists)"
# --- end Parked-items fix round 1

[ $fail -eq 0 ] && echo "PASS" || exit 1
