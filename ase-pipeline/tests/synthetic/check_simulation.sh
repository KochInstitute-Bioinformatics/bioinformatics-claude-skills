#!/bin/bash
# Usage: check_simulation.sh <outdir>   (run after simulate_ase_data.R)
set -u; D=${1:?usage: check_simulation.sh <outdir>}; fail=0
chk() { [ -s "$1" ] || { echo "FAIL: missing/empty $1"; fail=1; }; }
for m in f1 outbred; do
  chk $D/$m/samples.csv; chk $D/$m/truth_genes.tsv; chk $D/$m/truth_snps.tsv
  for s in $(tail -n +2 $D/$m/samples.csv 2>/dev/null | cut -d, -f1); do chk $D/$m/${s}_1.fastq.gz; chk $D/$m/${s}_2.fastq.gz; done
done
chk $D/f1/parental_snps.vcf
[ -e $D/genome/genome.fa ] && { chk $D/genome/genome.fa; chk $D/genome/genome.gtf; }
for i in $(tail -n +2 $D/outbred/samples.csv 2>/dev/null | cut -d, -f6 | sort -u); do chk $D/outbred/$i.het.vcf; chk $D/outbred/$i.all.vcf; done
# truth: at least 2 planted imbalanced genes and 4 null genes; alt fraction values in [0,1]
awk -v m=${MIN_NULL:-4} -F'\t' 'NR>1{c[$4]++; if($3<0||$3>1) bad=1} END{ if(c["imbalanced"]<2||c["null"]<m||bad) exit 1 }' $D/f1/truth_genes.tsv \
  || { echo "FAIL: f1 truth_genes.tsv needs >=2 imbalanced and >=${MIN_NULL:-4} null genes with p_alt in [0,1]"; fail=1; }
# reads: 4 lines per record, equal read counts in mates, read length 100
for f in $D/f1/*_1.fastq.gz; do
  n1=$(zcat $f | awk 'END{print NR/4}'); n2=$(zcat ${f%_1.fastq.gz}_2.fastq.gz | awk 'END{print NR/4}')
  [ "$n1" = "$n2" ] && [ "$n1" -ge 500 ] || { echo "FAIL: $f read counts ($n1 vs $n2)"; fail=1; }
  zcat $f | awk 'NR%4==2{print length($0)}' | sort -u | grep -qx 100 || { echo "FAIL: $f read length not 100"; fail=1; }
done
# determinism marker
[ -s $D/seed.txt ] || { echo "FAIL: missing seed.txt"; fail=1; }
[ $fail -eq 0 ] && echo PASS || exit 1
