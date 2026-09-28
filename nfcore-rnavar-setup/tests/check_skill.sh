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
allow="rename-chrs dependency runMode genomeDir genomeFastaFiles sjdbGTFfile sjdbOverhang runThreadN mail-type mail-user mem"
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
need "--read_length"
need "most common read length"
need "sjdbOverhang"
forbid "strandedness,"
# --- end Task 2

# --- Task 3
need "## Step 6"
need "star_rnavar_sjdb"
need "GTF source"
need "no flag is emitted"
need "## Step 7"
need "--dbsnp"
need "--dbsnp_tbi"
need "--known_indels"
need "--known_indels_tbi"
need "--skip_baserecalibration"
need "does not skip base recalibration automatically"
need "resolve the resource URLs at run time"
need "--star_index"
need "Store:"
need "use the highest N unless the user asks otherwise"
need "still apply the STAR index rule below"
need "{STAR_INDEX} is the same path"
# --- end Task 3

# --- Task 4
need "## Step 8"
need "--remove_duplicates"
need "--star_twopass"
need "--gatk_hc_call_conf"
need "--gatk_vf_qd_filter"
need "--gatk_vf_fs_filter"
need "--gatk_vf_window_size"
need "--gatk_vf_cluster_size"
need "--skip_variantfiltration"
need "--generate_gvcf"
need "--bam_csi_index"
need "## Step 9"
need "--tools"
need "--snpeff_cache"
need "--vep_cache"
need "--snpeff_db"
need "--vep_genome"
need "--vep_species"
need "--vep_cache_version"
need "--download_cache"
need "needs internet from compute nodes"
need "not a parameter in the rnavar schema"
# --- end Task 4

# --- Task 5
need "## Step 10"
need "do not overwrite it"
need "resourceLimits"
need "'.*:STAR_ALIGN'"
need "'.*:GATK4_HAPLOTYPECALLER'"
need "'.*:GATK4_BASERECALIBRATOR'"
need "'.*:GATK4_SPLITNCIGARREADS'"
need "## Step 11"
need "nextflow run nf-core/rnavar -r {VERSION}"
need "--seq_platform illumina"
need "--read_length {READ_LENGTH}"
need "--star_index '{STAR_INDEX}'"
need "## Step 12"
need "build_star_index_rnavar"
need "gunzip -c"
need "## Step 13"
need "## Notes for the assistant"
need "read-only"
need "never leave a bare"
need "already ends in"
need "{ENS_VERSION}"
need "where \`{TOOL}\` is"
forbid '{KNOWN_SITES_LINES}\'
forbid '{VARIANT_LINES}\'
forbid '{ANNOTATION_LINES}\'
forbid "{VERSION_ENS}"
# --- end Task 5


# --- Final-review fixes
need "bgzipped \`.vcf.gz\` files with \`.tbi\` indexes"
need "contig-name check"
need "zcat FILE.vcf.gz | grep -v '^#' | head -1 | cut -f1"
need "grep -m1 '^>' {FASTA_PATH} | cut -d' ' -f1 | sed 's/^>//'"
need "prefer Ensembl-named variation VCFs"
need "bcftools annotate --rename-chrs MAP.txt"
need "chr1<->1 ... chrM<->MT"
need "never proceed with mismatched contigs"
need "run \`bgzip\` on any plain \`.vcf\` before \`tabix -p vcf\`"
need "1. these are lanes of the same sample"
need "2. these are different samples"
need "allow the deliberate duplicates"
need "Emit \`--download_cache\` only if the user explicitly chooses it after being warned."
forbid "Use \`--snpeff_cache\`, \`--vep_cache\` and \`--download_cache\` only."
need "--dependency=afterok:<helper_jobid> nf-core_rnavar_{VERSION}.sh"
need "with \`--skip_baserecalibration\` these are the duplicate-marked BAMs"
need "check whether \`{SAMPLESHEET_CSV}\` already exists"
need "check whether \`nf-core_rnavar_{VERSION}.sh\` already exists"
need "1. overwrite · 2. choose another filename"
# --- end Final-review fixes

# --- Contig-guard fix
need "always-run contig guard"
need 'VCF_CONTIG=$(zcat'
need 'FASTA_CONTIG=$(grep -m1'
need 'exit 1'
need "the pipeline must not be run"
need "chrM\` to \`MT\`"
need "helper's contig guard performs the check after download"
need "Mouse Genomes Project VCFs use Ensembl-style contig names"
forbid "run the check once they are"
forbid "If the Step 7 contig-name check found a mismatch and option (ii) was chosen"
# --- end Contig-guard fix

[ $fail -eq 0 ] && echo "PASS" || exit 1
