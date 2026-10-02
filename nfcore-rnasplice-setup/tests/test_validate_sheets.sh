#!/bin/bash
# Usage: test_validate_sheets.sh <skill.md>
# Runs validate_rnasplice_sheets from the skill (block after "**Sheet validation.**") on good and bad sheets.
# Prints "VALIDATE PASS" or exits 1.
set -u
source "$(cd "$(dirname "$0")" && pwd)/test_env.sh" "$@"
SKILL=${1:?usage: test_validate_sheets.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd); FIX="$HERE/fixtures"
gv() { awk -F'\t' -v k="$1" '$1 == k {print $2; exit}' "$FIX/gate_values.tsv"; }
CODE=$(bash "$HERE/cut_block.sh" "$SKILL" "**Sheet validation.**") || { echo "FAIL: sheet validation block not found"; exit 1; }
eval "$CODE"
declare -F validate_rnasplice_sheets >/dev/null || { echo "FAIL: validate_rnasplice_sheets is not defined by the block"; exit 1; }
S="$TEST_TMP/s.csv"; C="$TEST_TMP/c.csv"; P="$TEST_TMP/p.csv"; fail=0
# ok <case> <source> <paired 0|1> [pairs file]; bad <case> <source> <paired 0|1> <expected error text> [pairs file]
ok()  { local out; out=$(validate_rnasplice_sheets "$S" "$C" "$2" "$3" ${4:+"$4"}); [ $? -eq 0 ] && [ "$out" = "SHEETS OK" ] || { echo "FAIL: case $1: expected SHEETS OK, got: $out"; fail=1; }; }
bad() { local out; out=$(validate_rnasplice_sheets "$S" "$C" "$2" "$3" ${5:+"$5"}); [ $? -ne 0 ] && grep -qF -- "$4" <<< "$out" || { echo "FAIL: case $1: expected an error containing '$4', got: $out"; fail=1; }; }
good_fastq() { printf '%s\n' "sample,fastq_1,fastq_2,strandedness,condition" "WT_1,a_R1.fastq.gz,a_R2.fastq.gz,reverse,WT" "WT_2,b_R1.fastq.gz,b_R2.fastq.gz,reverse,WT" "KO_1,c_R1.fastq.gz,c_R2.fastq.gz,reverse,KO" "KO_2,d_R1.fastq.gz,d_R2.fastq.gz,reverse,KO" > "$S"; }
good_con() { printf '%s\n' "contrast,treatment,control" "KO_vs_WT,KO,WT" > "$C"; }
pairs() { printf '%s\n' "sample,subject" "$@" > "$P"; }

cp "$FIX/test_samplesheet.csv" "$S"; cp "$FIX/test_contrastsheet.csv" "$C"; ok upstream_test_sheets fastq 0
good_fastq; good_con; ok good fastq 0
good_fastq; printf '%s\n' "contrast,treatment,control" "KO_vs_WT,K0,WT" > "$C"; bad label_typo_treatment fastq 0 'treatment "K0" of contrast KO_vs_WT is not a value of the condition column'
good_fastq; printf '%s\n' "contrast,treatment,control" "KO_vs_WT,KO,wt" > "$C"; bad label_typo_control fastq 0 'control "wt" of contrast KO_vs_WT is not a value of the condition column'
good_fastq; printf '%s\n' "KO_vs_WT,KO,WT" > "$C"; bad contrasts_no_header fastq 0 'contrasts header is'
good_fastq; printf '%s\n' "contrast,treatment,control" "KO_vs_KO,KO,KO" > "$C"; bad self_contrast fastq 0 'compares KO with itself'
good_fastq; printf '%s\n' "contrast,treatment,control" "KO_vs_WT,KO,WT" "KO_vs_WT,WT,KO" > "$C"; bad duplicate_contrast fastq 0 'contrast name KO_vs_WT is used twice'
good_fastq; printf '%s\n' "contrast,treatment,control" "KO_vs_WT,KO,WT" "KO_WT,KO,WT" > "$C"; bad duplicate_pair fastq 0 'contrasts KO_vs_WT and KO_WT both compare KO with WT'
good_fastq; printf '%s\n' "contrast,treatment,control" > "$C"; bad no_contrast_rows fastq 0 'has no contrast rows'
good_fastq; sed -i '2s/,reverse,/,forward,/' "$S"; good_con; bad mixed_strandedness fastq 0 'mixed strandedness'
good_fastq; sed -i 's/,reverse,/,auto,/' "$S"; good_con; bad strandedness_auto fastq 0 'strandedness "auto"'
good_fastq; sed -i '3s/,b_R2.fastq.gz,/,,/' "$S"; good_con; bad mixed_layout fastq 0 'mixed single-end and paired-end'
good_fastq; sed -i '5d' "$S"; good_con; bad one_sample_condition fastq 0 'condition KO has 1 sample'
good_fastq; sed -i '2s/^WT_1/WT-1/' "$S"; good_con; bad dash_in_name fastq 0 'sample name "WT-1"'
good_fastq; sed -i '2s/^WT_1/1WT/' "$S"; good_con; bad digit_first fastq 0 'sample name "1WT"'
good_fastq; sed -i '2s/^WT_1/NA/' "$S"; good_con; bad reserved_word fastq 0 'sample name "NA" is an R reserved word'
good_fastq; sed -i '2s/,WT$/,wild-type/' "$S"; good_con; bad bad_condition fastq 0 'condition "wild-type"'
good_fastq; printf '%s\n' "WT_1,a2_R1.fastq.gz,a2_R2.fastq.gz,reverse,WT" >> "$S"; good_con; ok tech_replicates fastq 0
good_fastq; printf '%s\n' "WT_1,a2_R1.fastq.gz,a2_R2.fastq.gz,reverse,KO" >> "$S"; good_con; bad tech_replicate_conflict fastq 0 'rows of sample WT_1 have different conditions'
good_fastq; sed -i '3s/^WT_2/WT_1/' "$S"; good_con; bad replicates_hide_size fastq 0 'condition WT has 1 sample'
good_fastq; sed -i '2s/a_R1.fastq.gz/a_R1.fq/' "$S"; good_con; bad not_gz fastq 0 'does not end in .fastq.gz or .fq.gz'
good_fastq; sed -i '2s/,a_R2.fastq.gz,/,"a_R2.fastq.gz",/' "$S"; good_con; bad quoted fastq 0 'contains a double quote'
good_fastq; good_con; sed -i 's/$/\r/' "$S" "$C"; ok crlf fastq 0
good_fastq; sed -i '1s/^sample,/Sample,/' "$S"; good_con; bad samplesheet_header fastq 0 'samplesheet header is'
# Paired design: a subject per sample; within each condition the i-th sample has the same subject.
good_fastq; good_con; pairs WT_1,m1 WT_2,m2 KO_1,m1 KO_2,m2; ok paired_two_by_two fastq 1 "$P"
good_fastq; good_con; bad paired_no_pairs_file fastq 1 'a paired design needs the pairs file'
good_fastq; good_con; pairs WT_1,m1 WT_2,m2 KO_1,m2 KO_2,m1; bad paired_order_mismatch fastq 1 'sort the rows of each condition by subject' "$P"
good_fastq; sed -i '4{h;d};5G' "$S"; good_con; pairs WT_1,m1 WT_2,m2 KO_1,m2 KO_2,m1; ok paired_order_sorted fastq 1 "$P"
good_fastq; good_con; pairs WT_1,m1 WT_2,m2 KO_1,m1; bad paired_missing_subject fastq 1 'sample KO_2 has no subject in the pairs file' "$P"
good_fastq; good_con; pairs WT_1,m1 WT_2,m2 KO_1,m1 KO_2,m2 KO_9,m3; bad paired_unknown_sample fastq 1 'pairs file names sample KO_9, which is not in the samplesheet' "$P"
good_fastq; good_con; pairs WT_1,m1 WT_2,m2 KO_1,m1 KO_2,m2 KO_2,m3; bad paired_sample_twice fastq 1 'sample KO_2 appears twice in the pairs file' "$P"
good_fastq; good_con; pairs WT_1,m1 WT_2,m1 KO_1,m1 KO_2,m2; bad paired_subject_repeated fastq 1 'subject m1 has 2 samples in condition WT' "$P"
good_fastq; good_con; pairs WT_1,m-1 WT_2,m2 KO_1,m-1 KO_2,m2; bad paired_bad_subject fastq 1 'subject "m-1"' "$P"
good_fastq; good_con; printf '%s\n' "sample,pair" WT_1,m1 WT_2,m2 KO_1,m1 KO_2,m2 > "$P"; bad paired_pairs_header fastq 1 'pairs header is' "$P"
good_fastq; printf '%s\n' "WT_1,a2_R1.fastq.gz,a2_R2.fastq.gz,reverse,WT" >> "$S"; good_con; pairs WT_1,m1 WT_2,m2 KO_1,m1 KO_2,m2; bad paired_split_replicates fastq 1 'rows of sample WT_1 are not next to each other' "$P"
good_fastq; sed -i '2a WT_1,a2_R1.fastq.gz,a2_R2.fastq.gz,reverse,WT' "$S"; good_con; pairs WT_1,m1 WT_2,m2 KO_1,m1 KO_2,m2; ok paired_adjacent_replicates fastq 1 "$P"
good_fastq; printf '%s\n' "HET_1,e_R1.fastq.gz,e_R2.fastq.gz,reverse,HET" "HET_2,f_R1.fastq.gz,f_R2.fastq.gz,reverse,HET" >> "$S"; good_con; bad paired_three_conditions fastq 1 'exactly two conditions'
good_fastq; printf '%s\n' "KO_3,e_R1.fastq.gz,e_R2.fastq.gz,reverse,KO" >> "$S"; good_con; bad paired_unequal fastq 1 'same number of samples in both conditions'
: > "$S"; good_con; bad empty_samplesheet fastq 0 'missing or empty'
# Review fixes: reserved condition labels, argument checks, clearer messages
good_fastq; sed -i 's/,KO$/,NA/' "$S"; printf '%s\n' "contrast,treatment,control" "NA_vs_WT,NA,WT" > "$C"; bad reserved_condition fastq 0 'condition "NA" of KO_1 is an R reserved word'
good_fastq; good_con; bad bad_source FASTQ 0 'source must be fastq or genome_bam'
good_fastq; good_con; bad bad_paired_flag fastq true 'paired design must be 0 or 1'
good_fastq; sed -i '1d' "$S"; good_con; out=$(validate_rnasplice_sheets "$S" "$C" fastq 0)
grep -qF 'samplesheet header is' <<< "$out" && ! grep -qF 'has 1 sample' <<< "$out" || { echo "FAIL: case header_missing_no_size_noise: got: $out"; fail=1; }
good_fastq; good_con; pairs WT_1,m1 WT_2,m2 KO_1,m1 KO_2,m3; bad subject_unpartnered fastq 1 'subject m3 has no sample in condition WT' "$P"
# genome-BAM sheets with the header recorded by the verification run
BH=$(gv BAM_SHEET_HEADER)
bam_rows() { awk -F',' -v h="$BH" 'BEGIN { n = split(h, c, ","); print h }
  { for (i = 1; i <= n; i++) { v = (c[i] == "sample") ? $1 : (c[i] == "condition") ? $2 : (c[i] == "genome_bam") ? $3 : (c[i] == "strandedness") ? "reverse" : (c[i] == "single_end") ? "false" : "x"; printf "%s%s", v, (i < n ? "," : "\n") } }'; }
printf '%s\n' "WT_1,WT,a.bam" "WT_2,WT,b.bam" "KO_1,KO,c.bam" "KO_2,KO,d.bam" | bam_rows > "$S"; good_con; ok bam_good genome_bam 0
printf '%s\n' "WT_1,WT,a.bam" "WT_1,WT,b.bam" "KO_1,KO,c.bam" "KO_2,KO,d.bam" | bam_rows > "$S"; good_con; bad bam_duplicate genome_bam 0 'appears twice'
printf '%s\n' "WT_1,WT,a.cram" "WT_2,WT,b.bam" "KO_1,KO,c.bam" "KO_2,KO,d.bam" | bam_rows > "$S"; good_con; bad bam_not_bam genome_bam 0 'does not end in .bam'
printf '%s\n' "WT_1,WT,a.bam" "WT_2,WT,b.bam" "KO_2,KO,d.bam" "KO_1,KO,c.bam" | bam_rows > "$S"; good_con; pairs WT_1,m1 WT_2,m2 KO_1,m1 KO_2,m2; bad bam_paired_order genome_bam 1 'sort the rows of each condition by subject' "$P"
printf '%s\n' "WT_1,WT,a.bam" "WT_2,WT,b.bam" "KO_1,KO,c.bam" "KO_2,KO,d.bam" | bam_rows > "$S"; good_con; ok bam_paired genome_bam 1 "$P"
# install_rnasplice_sheets (same block): validate, check the input files, move into the current directory, always remove the scratch dir
W="$TEST_TMP/work"; mkdir -p "$W"; for f in a b c d; do echo x > "$W/${f}_R1.fastq.gz"; echo x > "$W/${f}_R2.fastq.gz"; done
declare -F install_rnasplice_sheets >/dev/null || { echo "FAIL: install_rnasplice_sheets is not defined by the block"; fail=1; }
# inst <case> <expected rc> <expected text> <arguments after the scratch dir>; installs copies of $S, $C (and $P when non-empty)
inst() { local name=$1 want=$2 text=$3 d out rc; shift 3
  d=$(mktemp -d "$TEST_TMP/scr.XXXXXX"); cp "$S" "$d/samplesheet.csv"; cp "$C" "$d/contrasts.csv"; [ -s "$P" ] && cp "$P" "$d/pairs.csv"
  rm -f "$W/out_s.csv" "$W/out_c.csv" "$W/out_p.csv"
  out=$(cd "$W" && install_rnasplice_sheets "$d" "$@"); rc=$?
  [ ! -e "$d" ] || { echo "FAIL: case $name: scratch directory left behind"; fail=1; }
  [ $rc -eq "$want" ] && grep -qF -- "$text" <<< "$out" || { echo "FAIL: case $name: expected rc $want and '$text', got rc $rc: $out"; fail=1; }; }
good_fastq; good_con; : > "$P"; inst install_ok 0 "WROTE out_s.csv out_c.csv" fastq 0 out_s.csv out_c.csv
cmp -s "$S" "$W/out_s.csv" && cmp -s "$C" "$W/out_c.csv" || { echo "FAIL: case install_ok: sheets not moved unchanged"; fail=1; }
good_fastq; printf '%s\n' "contrast,treatment,control" "KO_vs_WT,K0,WT" > "$C"; inst install_invalid 1 'is not a value of the condition column' fastq 0 out_s.csv out_c.csv
[ ! -e "$W/out_s.csv" ] && [ ! -e "$W/out_c.csv" ] || { echo "FAIL: case install_invalid: a sheet was moved although validation failed"; fail=1; }
good_fastq; good_con; mv "$W/d_R2.fastq.gz" "$W/d_R2.keep"; inst install_missing_input 1 'input file d_R2.fastq.gz is missing or empty' fastq 0 out_s.csv out_c.csv; mv "$W/d_R2.keep" "$W/d_R2.fastq.gz"
[ ! -e "$W/out_s.csv" ] || { echo "FAIL: case install_missing_input: a sheet was moved although an input file is missing"; fail=1; }
good_fastq; good_con; pairs WT_1,m1 WT_2,m2 KO_1,m1 KO_2,m2; inst install_paired_ok 0 "WROTE out_s.csv out_c.csv out_p.csv" fastq 1 out_s.csv out_c.csv out_p.csv
cmp -s "$P" "$W/out_p.csv" || { echo "FAIL: case install_paired_ok: pairs file not moved"; fail=1; }
good_fastq; good_con; pairs WT_1,m1 WT_2,m2 KO_1,m2 KO_2,m1; inst install_paired_order 1 'sort the rows of each condition by subject' fastq 1 out_s.csv out_c.csv out_p.csv
good_fastq; good_con; inst install_paired_no_name 1 'a paired design needs the pairs file name' fastq 1 out_s.csv out_c.csv
good_fastq; good_con; : > "$P"; inst install_bad_target 1 'must be a plain file name' fastq 0 sub/out_s.csv out_c.csv
# A custom prefix with a space ("KO vs WT"): quoted, the name check refuses it; unquoted, the split names exceed the argument
# count. Either way nothing is moved and the scratch directory is removed.
good_fastq; good_con; : > "$P"; inst install_space_quoted 1 'must be a plain file name in the current directory (letters, digits, ., _ and - only)' fastq 0 "KO vs WT_samplesheet.csv" "KO vs WT_contrasts.csv"
good_fastq; good_con; : > "$P"; inst install_space_split 1 'arguments, at most 6: a file name was split at a space' fastq 0 KO vs WT_samplesheet.csv KO vs WT_contrasts.csv
[ -z "$(ls -A "$W" | grep -v -e '_R[12]\.fastq\.gz$' -e '^scr\.')" ] || { echo "FAIL: case install_space: a file was moved: $(ls -A "$W")"; fail=1; }
good_fastq; good_con; : > "$P"; inst install_quote 1 'letters, digits, ., _ and - only' fastq 0 "KO\"x_samplesheet.csv" out_c.csv
forbidden_ran && fail=1
[ $fail -eq 0 ] && echo "VALIDATE PASS" || exit 1
