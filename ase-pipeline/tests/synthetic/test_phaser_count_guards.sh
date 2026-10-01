#!/bin/bash
# Usage: bash test_phaser_count_guards.sh <skill.md> <finished project dir> <work dir (must not exist)>
# Negative and safety tests of phaser_count.sh (Step 18) on real data. <finished project dir> is the CWD of a finished UNPHASED
# outbred smoke run (PHASED_GT 0) that holds values.tsv and genotype_map.tsv (as written for assemble_outbred_from_skill.sh) and
# whose RESULTS_DIR holds bam/, genotypes/, reference/genes_span.bed and phaser/. The project is copied to <work dir>/proj (the
# source is never modified), its scripts are regenerated from the skill with the assembler, and every case runs array task 1
# through `sbatch -p bcc --wait` (phASER, samtools and bcftools only run in jobs; this script only edits text and calls sbatch).
# Cases (task 1 = first data row of the sheet, sample S1):
#   1. s1.redo outputs survive: finished outputs of a sample "S1.redo" survive a successful re-run of S1 (no glob removal)
#   2. contig mismatch with stale outputs: gene spans renamed -> "lies on a contig of the BAM", no S1 output left
#   3. false phased claim: PHASED_GT=1 with unphased VCFs -> "PHASED_GT=1 but only 0 of"
#   4. contig name with '_': S1's het VCF on contig "<contig>_x" -> "contigs whose names contain '_'"
#   5. install marker with only line 1 -> "no complete phASER installation"
#   6. sample name "x/../../bam" -> "must not contain '/'", and the BAM directory is untouched
#   7. PHASED_GT=yes -> "PHASED_GT must be 0 or 1"
#   8. single-end BAM -> "holds no paired reads"
#   9. partial phasing under PHASED_GT=0: S1's het VCF with every other genotype written 0|1 and a PS FORMAT tag (as GATK
#      HaplotypeCaller / rnavar write physical phasing); phASER must get a copy without phase, so the run exits 0, its gene
#      tables identical to the unphased run (case 1) on the same BAM, the genotype VCF unchanged and no tmp/phaser_S1 left
#  10. PS tag under PHASED_GT=1: every genotype 0|1 (passes the 90 % rule) with a PS FORMAT tag -> "carries a PS FORMAT tag"
# Case 1 also checks that the gene table (every column but the BAM path) is equal to the source run's table of S1 (the
# unphased copy that PHASED_GT=0 now gives phASER changes nothing for a VCF without '|').
# Every failing case must exit 1 with its message and leave none of S1's phASER outputs. Ends with "ALL GUARD TESTS OK" (exit 0)
# or "SOME GUARD TESTS BAD" (exit 1).
set -uo pipefail
SKILL=${1:?usage}; SRC=${2:?usage}; WORK=${3:?usage}; HERE=$(cd "$(dirname "$0")" && pwd)
die() { echo "ERROR: $*" >&2; exit 1; }
case "$WORK" in /*) ;; *) die "work dir must be an absolute path" ;; esac
[ ! -e "$WORK" ] || die "$WORK exists; give a new work dir"
[ -s "$SRC/values.tsv" ] && [ -s "$SRC/genotype_map.tsv" ] || die "$SRC has no values.tsv or genotype_map.tsv"
mkdir -p "$WORK" && cp -r "$SRC" "$WORK/proj" || die "cannot copy $SRC"
P=$WORK/proj
awk -F'\t' -v a="$SRC" -v b="$P" 'BEGIN { OFS = "\t" } { if (index($2, a) == 1) $2 = b substr($2, length(a) + 1); print }' "$SRC/values.tsv" > "$P/values.tsv" || die "values.tsv"
val() { awk -F'\t' -v k="$1" '$1 == k {print $2; f = 1} END {exit !f}' "$P/values.tsv" || die "values.tsv has no $1"; }
RD=$(val RESULTS_DIR) || exit 1; CSV=$(val SAMPLES_CSV) || exit 1; PH=$(val PHASER_HOME) || exit 1
[ "$(val PHASED_GT)" = 0 ] || die "the source project must be the unphased run (PHASED_GT 0)"
case "$RD" in "$P"/*) ;; *) die "RESULTS_DIR $RD is not inside the copied project" ;; esac
bash "$HERE/assemble_outbred_from_skill.sh" "$SKILL" "$P/values.tsv" "$P/genotype_map.tsv" > /dev/null || die "assembler failed"
SC=$RD/scripts/phaser_count.sh; V=$WORK/variants; mkdir -p "$V" || exit 1
S1=$(awk -F, 'NR == 2 {print $1}' "$CSV"); I1=$(awk -F, 'NR == 2 {print $6}' "$CSV" | tr -d '\r')
[ -n "$S1" ] && [ -n "$I1" ] || die "no first data row in $CSV"
SFX="haplotypic_counts.txt haplotypes.txt allelic_counts.txt allele_config.txt variant_connections.txt vcf.gz vcf.gz.tbi gene_ae.txt gene_ae.txt.part"
bad=0
# leftover count, taken on a compute node: the login node's NFS attribute cache can still show files that a job has just removed
# (seen in the first run of this test: case 2 reported 8 outputs left that were already gone)
NL=0
left() { NL=$((NL + 1)); helper "left_$NL" "n=0; for x in $SFX; do [ -e '$RD/phaser/$S1.'\$x ] && n=\$((n + 1)); done; echo \$n > '$WORK/left_$NL.txt'" > /dev/null 2>&1; LEFT=$(cat "$WORK/left_$NL.txt" 2>/dev/null); LEFT=${LEFT:-unknown}; }
job() {   # $1 script: runs array task 1 and waits; sets J (job id) and ST (exit code); prints the task log path in LOGF
  J=$(sbatch -p bcc --parsable --wait --array=1 "$1"); ST=$?; J=${J%%;*}; LOGF=$RD/logs/phaser_count_${J}_1.out; }
expect_fail() {   # $1 case name, $2 script, $3 message
  job "$2"; left; local n=$LEFT
  if [ "$ST" -ne 0 ] && grep -qF -- "$3" "$LOGF" 2>/dev/null && [ "$n" -eq 0 ]; then echo "OK   $1: job $J exit $ST, '$3', no $S1 output left"
  else echo "BAD  $1: job $J exit $ST, $n $S1 outputs left; log: $(tail -2 "$LOGF" 2>/dev/null | paste -sd' ')"; bad=1; fi; }
helper() {   # $1 name, $2 commands (run in a job with the phASER environment's tools on PATH)
  printf '#!/bin/bash\n#SBATCH -N 1 -n 1 --mem=2G -t 0:15:00\n#SBATCH -o %s/%s_%%j.out\nset -uo pipefail\nexport PATH="%s/env/bin:$PATH"\n%s\n' \
    "$WORK" "$1" "$PH" "$2" > "$V/$1.sh"
  sbatch -p bcc --parsable --wait "$V/$1.sh" > /dev/null || { echo "BAD  helper $1 failed (see $WORK/$1_*.out)"; bad=1; return 1; }; }

# 1. s1.redo outputs survive a successful re-run of S1
for x in gene_ae.txt haplotypic_counts.txt; do echo "finished $S1.redo" > "$RD/phaser/$S1.redo.$x"; done
job "$SC"
if [ "$ST" -eq 0 ] && [ -s "$RD/phaser/$S1.gene_ae.txt" ] && [ ! -e "$RD/phaser/$S1.gene_ae.txt.part" ] \
   && [ "$(cat "$RD/phaser/$S1.redo.gene_ae.txt" "$RD/phaser/$S1.redo.haplotypic_counts.txt")" = "finished $S1.redo
finished $S1.redo" ]; then echo "OK   1 s1.redo outputs survive a re-run of s1: job $J exit 0, $S1.gene_ae.txt written, $S1.redo.* unchanged"
else echo "BAD  1 s1.redo outputs survive: job $J exit $ST; redo files: $(ls "$RD/phaser/$S1.redo".* 2>/dev/null | wc -l)"; bad=1; fi
# 1b. the gene table of case 1 (every column but the BAM path, column 12) is equal to the source run's; it is kept as the
#     unphased reference of case 9. Compared in a job (login-node NFS cache, above).
SRC_RD="$SRC${RD#"$P"}"
helper cmp_src "[ -s '$SRC_RD/phaser/$S1.gene_ae.txt' ] && cp '$RD/phaser/$S1.gene_ae.txt' '$V/gene_ae.case1' && if [ \"\$(cut -f1-11 '$SRC_RD/phaser/$S1.gene_ae.txt' | md5sum)\" = \"\$(cut -f1-11 '$V/gene_ae.case1' | md5sum)\" ]; then echo same; else echo differ; fi > '$WORK/cmp_src.txt'" \
  && [ "$(cat "$WORK/cmp_src.txt" 2>/dev/null)" = same ] && echo "OK   1b gene table of $S1 equal to the source run's (columns 1-11; the source run used the VCF itself)" \
  || { echo "BAD  1b gene table of $S1 differs from the source run's: $(cat "$WORK/cmp_src.txt" 2>/dev/null)"; bad=1; }
# 2. contig mismatch, with the outputs of case 1 as stale outputs
BED=$RD/reference/genes_span.bed; cp "$BED" "$V/genes_span.bed.orig" || exit 1
awk 'BEGIN { OFS = "\t" } { $1 = "X" $1; print }' "$V/genes_span.bed.orig" > "$BED"
left; [ "$LEFT" -gt 0 ] || { echo "BAD  2: no stale output from case 1"; bad=1; }
expect_fail "2 contig mismatch with stale outputs" "$SC" "lies on a contig of the BAM"
cp "$V/genes_span.bed.orig" "$BED" || exit 1
# 3. false phased claim
sed 's/^\(PHASER_HOME=.*; PHASED_GT=\)0$/\11/' "$SC" > "$V/phased_claim.sh"; grep -q 'PHASED_GT=1$' "$V/phased_claim.sh" || { echo "BAD  3: no PHASED_GT line"; bad=1; }
expect_fail "3 false phased claim" "$V/phased_claim.sh" "PHASED_GT=1 but only 0 of"
# 4. contig name with '_' in S1's het VCF (copy kept and restored)
G=$RD/genotypes/$I1.het.vcf.gz; cp "$G" "$V/het.vcf.gz.orig" && cp "$G.tbi" "$V/het.vcf.gz.tbi.orig" || exit 1
helper rename_contig "bcftools view '$V/het.vcf.gz.orig' | awk 'BEGIN { OFS = \"\\t\" } /^##contig=<ID=/ { sub(/<ID=[^,>]+/, \"&_x\"); print; next } /^#/ { print; next } { \$1 = \$1 \"_x\"; print }' | bgzip -c > '$G' && tabix -f -p vcf '$G'" \
  && expect_fail "4 contig name with '_'" "$SC" "contigs whose names contain '_'"
cp "$V/het.vcf.gz.orig" "$G" && cp "$V/het.vcf.gz.tbi.orig" "$G.tbi" || exit 1
# 5. install marker with only line 1
F=$WORK/fake_phaser_home; mkdir -p "$F/env/bin" && ln -s "$PH/env/bin/python" "$F/env/bin/python" && head -1 "$PH/install_ok.txt" > "$F/install_ok.txt" || exit 1
sed "s#^PHASER_HOME=\"[^\"]*\"#PHASER_HOME=\"$F\"#" "$SC" > "$V/marker_line1.sh"
expect_fail "5 install marker with only line 1" "$V/marker_line1.sh" "no complete phASER installation"
# 6. sample name with path characters: the BAM directory must stay untouched
NB=$(ls "$RD/bam" | wc -l)
awk -F, 'BEGIN { OFS = "," } NR == 2 { $1 = "x/../../bam" } { print }' "$CSV" > "$V/bad_sample.csv"
awk -v a="$CSV" -v b="$V/bad_sample.csv" '{ while ((i = index($0, a)) > 0) $0 = substr($0, 1, i - 1) b substr($0, i + length(a)); print }' "$SC" > "$V/bad_sample.sh"
expect_fail "6 sample name x/../../bam (must not contain '/')" "$V/bad_sample.sh" "must not contain '/'"
[ "$(ls "$RD/bam" | wc -l)" -eq "$NB" ] && [ -s "$RD/bam/$S1.bam" ] && echo "OK   6b BAM directory untouched ($NB files)" || { echo "BAD  6b BAM directory changed"; bad=1; }
# 7. PHASED_GT neither 0 nor 1
sed 's/^\(PHASER_HOME=.*; PHASED_GT=\)0$/\1yes/' "$SC" > "$V/phased_yes.sh"
expect_fail "7 PHASED_GT=yes" "$V/phased_yes.sh" "PHASED_GT must be 0 or 1"
# 8. single-end BAM (pairing flags 1, 2, 8, 32, 64, 128 cleared; copy kept and restored)
B=$RD/bam/$S1.bam; cp "$B" "$V/bam.orig" && cp "$B.bai" "$V/bam.bai.orig" || exit 1
helper single_end "samtools view -h '$V/bam.orig' | awk 'BEGIN { OFS = \"\\t\" } /^@/ { print; next } { f = \$2; n = split(\"1 2 8 32 64 128\", b, \" \"); for (i = 1; i <= n; i++) if (int(f / b[i]) % 2 == 1) f -= b[i]; \$2 = f; print }' | samtools view -b -o '$B' - && samtools index '$B'" \
  && expect_fail "8 single-end BAM" "$SC" "holds no paired reads"
cp "$V/bam.orig" "$B" && cp "$V/bam.bai.orig" "$B.bai" || exit 1
# 9. partial phasing (every other genotype 0|1, with a PS FORMAT tag) under PHASED_GT=0: phASER must get a copy without phase,
#    so the gene table equals the unphased run of case 1 on the same BAM, and the genotype VCF (STAR's file) is not changed
# awk programs for the two VCFs: every other genotype (case 9) or every genotype (case 10) written 0|1, with a PS FORMAT tag
AWK_HDR='BEGIN { OFS = "\t" } /^##/ { print; next } /^#CHROM/ { print "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">"; print; next }'
printf '%s\n' "$AWK_HDR" '{ n++; if (n % 2 == 1) { $9 = $9 ":PS"; sub(/^0\/1/, "0|1", $10); $10 = $10 ":" $2 } print }' > "$V/partial.awk" || exit 1
printf '%s\n' "$AWK_HDR" '{ $9 = $9 ":PS"; sub(/^0\/1/, "0|1", $10); $10 = $10 ":1"; print }' > "$V/ps_all.awk" || exit 1
helper partial_vcf "bcftools view '$V/het.vcf.gz.orig' | awk -F'\\t' -f '$V/partial.awk' | bgzip -c > '$G' && tabix -f -p vcf '$G' && md5sum < '$G' > '$WORK/partial_md5_before.txt' && bcftools query -f '[%GT]\\n' '$G' | grep -c '|' > '$WORK/partial_nph.txt'" || true
NPH=$(cat "$WORK/partial_nph.txt" 2>/dev/null); NPH=${NPH:-0}
[ "$NPH" -gt 0 ] && echo "OK   9a partial VCF built: $NPH phased genotypes with PS" || { echo "BAD  9a the partial VCF has no phased genotype"; bad=1; }
job "$SC"
helper check9 "{ if [ \"\$(cut -f1-11 '$RD/phaser/$S1.gene_ae.txt' | md5sum)\" = \"\$(cut -f1-11 '$V/gene_ae.case1' | md5sum)\" ]; then echo ga_same; else echo ga_differ; fi
  if [ \"\$(md5sum < '$G')\" = \"\$(cat '$WORK/partial_md5_before.txt')\" ]; then echo vcf_same; else echo vcf_changed; fi
  if [ -e '$RD/tmp/phaser_$S1' ]; then echo tmp_left; else echo tmp_gone; fi
  grep -c 'phASER gets a copy without phase' '$LOGF'; } > '$WORK/check9.txt'" || true
C9=$(paste -sd' ' "$WORK/check9.txt" 2>/dev/null)
if [ "$ST" -eq 0 ] && [ "$C9" = "ga_same vcf_same tmp_gone 1" ]; then
  echo "OK   9 partial phasing under PHASED_GT=0: job $J exit 0, gene tables identical to the unphased run (case 1), genotype VCF unchanged, temporary copy removed"
else echo "BAD  9 partial phasing under PHASED_GT=0: job $J exit $ST, checks: ${C9:-none} (expected: ga_same vcf_same tmp_gone 1); log: $(tail -2 "$LOGF" 2>/dev/null | paste -sd' ')"; bad=1; fi
# 10. a PS FORMAT tag under PHASED_GT=1: every genotype phased (the 90 % rule passes), so only the PS stop can fire
helper ps_vcf "bcftools view '$V/het.vcf.gz.orig' | awk -F'\\t' -f '$V/ps_all.awk' | bgzip -c > '$G' && tabix -f -p vcf '$G'" \
  && expect_fail "10 PS tag under PHASED_GT=1" "$V/phased_claim.sh" "carries a PS FORMAT tag"
cp "$V/het.vcf.gz.orig" "$G" && cp "$V/het.vcf.gz.tbi.orig" "$G.tbi" || exit 1
if [ "$bad" -eq 0 ]; then echo "ALL GUARD TESTS OK"; else echo "SOME GUARD TESTS BAD"; fi
exit $bad
