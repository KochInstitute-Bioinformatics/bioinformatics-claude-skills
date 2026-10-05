#!/bin/bash
# Test of the module-detection block of the skill (Step 3, anchor "**Module detection.**"): the block is cut out of the skill
# byte for byte, its placeholders {OUTDIR_RNASPLICE} and {CONTRASTS_CSV} are replaced by paths (as the wizard writes them out
# literally), and it runs under env -i on fixture trees built in a temp dir from fixtures/real_testdata. Only light shell runs.
# The FASTQ paths of the fixture sample sheet are never opened. Prints "DETECT PASS" or exits 1.
set -u
here=$(cd "$(dirname "$0")" && pwd)
skill=${1:-$here/../nfcore-rnasplice-report.md}
fx=$here/fixtures/real_testdata
tmp=$(mktemp -d "${TMPDIR:-/tmp}/rsr_detect.XXXXXX") || exit 1
trap 'case "$tmp" in */rsr_detect.??????) rm -rf -- "$tmp" ;; esac' EXIT
fail=0
check() { # name expected actual
  if [ "$2" != "$3" ]; then echo "FAIL $1: expected '$2' got '$3'"; fail=1; else echo "ok   $1"; fi
}
blk=$(bash "$here/cut_block.sh" "$skill" '**Module detection.**') || { echo "FAIL: the detection block is missing from the skill"; exit 1; }
[ -n "$blk" ] || { echo "FAIL: the detection block is empty"; exit 1; }
# run <tree> <contrasts csv>: the block with its placeholders filled in, under env -i; the output goes to $tmp/last
run() {
  local code=${blk//\{OUTDIR_RNASPLICE\}/$1}; code=${code//\{CONTRASTS_CSV\}/$2}
  env -i PATH=/usr/bin:/bin /bin/bash --noprofile --norc -c "$code" > "$tmp/last" 2>&1
  grep -c '^MODULES=' "$tmp/last" > "$tmp/nmod"
}
modules() { grep '^MODULES=' "$tmp/last"; }
# tree <dir> <flags>: copy real fixtures into the layout of a run (r rmats, s suppa, d dexseq_exon, e edger_exon, t dtu)
tree() {
  local d=$1 f=$2; mkdir -p "$d"
  case $f in *r*) mkdir -p "$d/star/rmats/GBR_vs_YRI"; cp "$fx"/star/rmats/GBR_vs_YRI/SE.MATS.JC.txt "$d/star/rmats/GBR_vs_YRI/" ;; esac
  case $f in *d*) mkdir -p "$d/star/dexseq_exon/results"; cp "$fx"/star/dexseq_exon/results/DEXSeqResults.GBR_vs_YRI.csv "$d/star/dexseq_exon/results/" ;; esac
  case $f in *e*) mkdir -p "$d/star/edger"; cp "$fx"/star/edger/contrast_GBR_vs_YRI.usage.exon.csv "$d/star/edger/" ;; esac
  case $f in *s*) mkdir -p "$d/salmon/suppa/diffsplice/per_local_event"; cp "$fx"/salmon/suppa/diffsplice/per_local_event/local_GBR-YRI.dpsi "$d/salmon/suppa/diffsplice/per_local_event/" ;; esac
  case $f in *t*) mkdir -p "$d/salmon/dexseq_dtu/results/dexseq"; cp "$fx"/salmon/dexseq_dtu/results/dexseq/DEXSeqResults.GBR-YRI.tsv "$d/salmon/dexseq_dtu/results/dexseq/" ;; esac
}
csv=$fx/261003_a_contrasts.csv

# the fixture facts the block must cope with: rMATS tables have the column ID twice
ids=$(head -n1 "$fx/star/rmats/GBR_vs_YRI/SE.MATS.JC.txt" | tr '\t' '\n' | grep -c '^ID$')
check fixture-rmats-ID-twice 2 "$ids"
check skill-says-ID-twice 1 "$(grep -c 'column `ID` twice' "$skill" | awk '{print ($1 > 0)}')"

# 1. all modules
tree "$tmp/all" rsdet; run "$tmp/all" "$csv"
check all "MODULES=rmats suppa dexseq_exon edger_exon dtu" "$(modules)"
check all-one-line 1 "$(cat "$tmp/nmod")"
# 2. rMATS only (the duplicate ID header must not trip the detection)
tree "$tmp/rm" r; run "$tmp/rm" "$csv"
check rmats-only "MODULES=rmats" "$(modules)"
# 3. only the JCEC file of rMATS present
mkdir -p "$tmp/jcec/star/rmats/GBR_vs_YRI"; cp "$fx"/star/rmats/GBR_vs_YRI/SE.MATS.JCEC.txt "$tmp/jcec/star/rmats/GBR_vs_YRI/"; run "$tmp/jcec" "$csv"
check rmats-jcec-only "MODULES=rmats" "$(modules)"
# 3b. paired run: rMATS tables only in star/rmats/{CONTRAST}_paired/ (rmats_paired_stats: true) count as found; header-only there does not
mkdir -p "$tmp/pr/star/rmats/GBR_vs_YRI_paired"; cp "$fx"/star/rmats/GBR_vs_YRI/SE.MATS.JC.txt "$tmp/pr/star/rmats/GBR_vs_YRI_paired/"; run "$tmp/pr" "$csv"
check rmats-paired-dir "MODULES=rmats" "$(modules)"
mkdir -p "$tmp/prh/star/rmats/GBR_vs_YRI_paired"; head -n1 "$fx"/star/rmats/GBR_vs_YRI/SE.MATS.JC.txt > "$tmp/prh/star/rmats/GBR_vs_YRI_paired/SE.MATS.JC.txt"; run "$tmp/prh" "$csv"
check rmats-paired-header-only "MODULES=" "$(modules)"
# 4. header-only files must NOT count (every module)
tree "$tmp/hdr" rsdet
for f in star/rmats/GBR_vs_YRI/SE.MATS.JC.txt star/dexseq_exon/results/DEXSeqResults.GBR_vs_YRI.csv star/edger/contrast_GBR_vs_YRI.usage.exon.csv \
  salmon/suppa/diffsplice/per_local_event/local_GBR-YRI.dpsi salmon/dexseq_dtu/results/dexseq/DEXSeqResults.GBR-YRI.tsv; do
  head -n1 "$fx/${f}" > "$tmp/hdr/$f"
done
run "$tmp/hdr" "$csv"
check header-only "MODULES=" "$(modules)"
# 4b. one header-only module among real ones: only that one drops out
tree "$tmp/hdr1" rsdet; head -n1 "$fx/star/edger/contrast_GBR_vs_YRI.usage.exon.csv" > "$tmp/hdr1/star/edger/contrast_GBR_vs_YRI.usage.exon.csv"
run "$tmp/hdr1" "$csv"
check header-only-one "MODULES=rmats suppa dexseq_exon dtu" "$(modules)"
# 4c. an empty file does not count
tree "$tmp/empty" r; : > "$tmp/empty/star/rmats/GBR_vs_YRI/SE.MATS.JC.txt"; run "$tmp/empty" "$csv"
check empty-file "MODULES=" "$(modules)"
# 5. BAM-input run: no salmon directory (the sample sheet schema is sample,condition,genome_bam,strandedness,single_end; detection never reads the sheet)
tree "$tmp/bam" rde
run "$tmp/bam" "$csv"
check bam-run "MODULES=rmats dexseq_exon edger_exon" "$(modules)"
check bam-run-says-no-salmon 1 "$(grep -c 'no salmon directory' "$tmp/last" | awk '{print ($1 > 0)}')"
# 6. wrong contrast name: nothing found, and the message names the missing path
printf 'contrast,treatment,control\nGBR_vs_CEU,GBR,CEU\n' > "$tmp/wrong.csv"
run "$tmp/all" "$tmp/wrong.csv"
check wrong-contrast "MODULES=" "$(modules)"
check wrong-contrast-names-path 1 "$(grep -cF "$tmp/all/star/rmats/GBR_vs_CEU" "$tmp/last" | awk '{print ($1 > 0)}')"
check wrong-contrast-names-suppa-path 1 "$(grep -cF "$tmp/all/salmon/suppa/diffsplice/per_local_event/local_GBR-CEU.dpsi" "$tmp/last" | awk '{print ($1 > 0)}')"
# 7. a contrast with - and _ in its group names: names come from the contrasts file, file names are never split
printf 'contrast,treatment,control\nGBR-1_vs_YRI_2,GBR-1,YRI_2\n' > "$tmp/odd.csv"
mkdir -p "$tmp/odd/star/rmats/GBR-1_vs_YRI_2" "$tmp/odd/salmon/suppa/diffsplice/per_local_event"
cp "$fx"/star/rmats/GBR_vs_YRI/SE.MATS.JC.txt "$tmp/odd/star/rmats/GBR-1_vs_YRI_2/"
cp "$fx"/salmon/suppa/diffsplice/per_local_event/local_GBR-YRI.dpsi "$tmp/odd/salmon/suppa/diffsplice/per_local_event/local_GBR-1-YRI_2.dpsi"
run "$tmp/odd" "$tmp/odd.csv"
check odd-names "MODULES=rmats suppa" "$(modules)"
# 8. no results at all: the block stops with a message and no MODULES line is a found module
mkdir -p "$tmp/none"; run "$tmp/none" "$csv"
check none "MODULES=" "$(modules)"
# 9. a missing contrasts file stops the block with a message
run "$tmp/all" "$tmp/does_not_exist.csv"
check missing-csv-no-modules-line 0 "$(cat "$tmp/nmod")"
check missing-csv-message 1 "$(grep -c 'ERROR' "$tmp/last" | awk '{print ($1 > 0)}')"
[ $fail -eq 0 ] && echo "DETECT PASS" || exit 1
