#!/bin/bash
# Test of the block that builds the expected file and column names from the contrasts file (Step 4, anchor
# "**Names built from the contrasts file.**"). A contrast whose treatment and control contain - and _ (GBR-1 and YRI_2) must give
# the exact names below; the names are never obtained by splitting file names. It also tests the sample sheet block of Step 4. Light shell only;
# prints "NAMES PASS" or exits 1.
set -u
here=$(cd "$(dirname "$0")" && pwd)
skill=${1:-$here/../nfcore-rnasplice-report.md}
tmp=$(mktemp -d "${TMPDIR:-/tmp}/rsr_names.XXXXXX") || exit 1
trap 'case "$tmp" in */rsr_names.??????) rm -rf -- "$tmp" ;; esac' EXIT
fail=0
blk=$(bash "$here/cut_block.sh" "$skill" '**Names built from the contrasts file.**') || { echo "FAIL: the names block is missing from the skill"; exit 1; }
printf 'contrast,treatment,control\nGBR-1_vs_YRI_2,GBR-1,YRI_2\n' > "$tmp/c.csv"
code=${blk//\{CONTRASTS_CSV\}/$tmp/c.csv}
got=$(env -i PATH=/usr/bin:/bin /bin/bash --noprofile --norc -c "$code" 2>&1)
want='CONTRAST=GBR-1_vs_YRI_2
TREATMENT=GBR-1
CONTROL=YRI_2
RMATS_DIR=star/rmats/GBR-1_vs_YRI_2
DEXSEQ_EXON_FILE=star/dexseq_exon/results/DEXSeqResults.GBR-1_vs_YRI_2.csv
EDGER_EXON_FILE=star/edger/contrast_GBR-1_vs_YRI_2.usage.exon.csv
SUPPA_FILE=salmon/suppa/diffsplice/per_local_event/local_GBR-1-YRI_2.dpsi
DTU_FILE=salmon/dexseq_dtu/results/dexseq/DEXSeqResults.GBR-1-YRI_2.tsv
DEXSEQ_COLUMN=log2fold_YRI_2_GBR-1
SUPPA_COLUMN=local_GBR-1-local_YRI_2_dPSI'
if [ "$got" = "$want" ]; then echo "ok   names"; else echo "FAIL names: expected"; echo "$want"; echo "got"; echo "$got"; fail=1; fi
# two contrasts: one block per contrast, in file order, no mixing
printf 'contrast,treatment,control\nA_vs_B,A,B\nC-d_vs_E,C-d,E\n' > "$tmp/c2.csv"
code=${blk//\{CONTRASTS_CSV\}/$tmp/c2.csv}
got=$(env -i PATH=/usr/bin:/bin /bin/bash --noprofile --norc -c "$code" 2>&1)
n=$(printf '%s\n' "$got" | grep -c '^CONTRAST=')
[ "$n" = 2 ] && echo "ok   two-contrasts" || { echo "FAIL two-contrasts: $n CONTRAST lines"; fail=1; }
printf '%s\n' "$got" | grep -qx 'DEXSEQ_COLUMN=log2fold_E_C-d' && echo "ok   second-column" || { echo "FAIL second-column"; fail=1; }
# a header-only or missing contrasts file stops with a message
code=${blk//\{CONTRASTS_CSV\}/$tmp/nope.csv}
got=$(env -i PATH=/usr/bin:/bin /bin/bash --noprofile --norc -c "$code" 2>&1)
printf '%s\n' "$got" | grep -q 'ERROR' && echo "ok   missing-file" || { echo "FAIL missing-file"; fail=1; }
# The sample sheet block (Step 4, anchor "Read the pipeline's own sample sheet"): it reads the sheets only and finds the condition
# column by name, for the FASTQ schema (the real fixture) and for the genome BAM schema; the files named in the sheets never exist here.
sblk=$(bash "$here/cut_block.sh" "$skill" "Read the pipeline's own sample sheet") || { echo "FAIL: the sample sheet block is missing from the skill"; exit 1; }
runsheet() { local code=${sblk//\{SAMPLESHEET\}/$1}; code=${code//\{CONTRASTS_CSV\}/$2}; env -i PATH=/usr/bin:/bin /bin/bash --noprofile --norc -c "$code" 2>&1; }
fx=$here/fixtures/real_testdata
got=$(runsheet "$fx/261003_a_samplesheet.csv" "$fx/261003_a_contrasts.csv")
for want in 'condition GBR 2' 'condition YRI 2' 'ok GBR' 'ok YRI'; do
  printf '%s\n' "$got" | grep -qx "$want" && echo "ok   sheet-fastq: $want" || { echo "FAIL sheet-fastq: missing line '$want'"; fail=1; }
done
printf 'sample,condition,genome_bam,strandedness,single_end\nS1,A,/nonexistent/a.bam,unstranded,0\nS2,B,/nonexistent/b.bam,unstranded,0\n' > "$tmp/bam.csv"
printf 'contrast,treatment,control\nA_vs_C,A,C\n' > "$tmp/c3.csv"
got=$(runsheet "$tmp/bam.csv" "$tmp/c3.csv")
printf '%s\n' "$got" | grep -qx 'condition A 1' && echo "ok   sheet-bam: condition A 1" || { echo "FAIL sheet-bam: condition A 1"; fail=1; }
printf '%s\n' "$got" | grep -qx 'MISSING in the sample sheet: C' && echo "ok   sheet-bam: group missing" || { echo "FAIL sheet-bam: group C not reported"; fail=1; }
[ $fail -eq 0 ] && echo "NAMES PASS" || exit 1
