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
allow="runMode genomeDir genomeFastaFiles sjdbGTFfile sjdbOverhang runThreadN mail-type mail-user mem"
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
# --- end Task 3

[ $fail -eq 0 ] && echo "PASS" || exit 1
