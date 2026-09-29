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
# reads: equal read counts in mates (>=500), every read of every mate (f1 and outbred) has length exactly 100
for m in f1 outbred; do
  for f in $D/$m/*_1.fastq.gz; do
    [ -e "$f" ] || continue
    n1=$(zcat $f | awk 'END{print NR/4}'); n2=$(zcat ${f%_1.fastq.gz}_2.fastq.gz | awk 'END{print NR/4}')
    [ "$n1" = "$n2" ] && [ "$n1" -ge 500 ] || { echo "FAIL: $f read counts ($n1 vs $n2)"; fail=1; }
    for g in $f ${f%_1.fastq.gz}_2.fastq.gz; do
      [ "$(zcat $g | awk 'NR%4==2{print length($0)}' | sort -u)" = 100 ] || { echo "FAIL: $g read length not exactly 100"; fail=1; }
    done
  done
done
# full-mode quotas (synthetic genome present and MIN_NULL unset): 4 imbalanced {0.70,0.70,0.85,0.30}, 1 bias, 7 null, bias/null p_alt 0.5
if [ -d $D/genome ] && [ -z "${MIN_NULL:-}" ]; then
  for m in f1 outbred; do
    t=$D/$m/truth_genes.tsv
    imb=$(awk -F'\t' 'NR>1&&$4=="imbalanced"{printf "%.2f\n",$3}' $t | sort | paste -sd' ')
    nb=$(awk -F'\t' 'NR>1&&$4=="bias"' $t | wc -l); nn=$(awk -F'\t' 'NR>1&&$4=="null"' $t | wc -l)
    [ "$imb" = "0.30 0.70 0.70 0.85" ] && [ "$nb" = 1 ] && [ "$nn" = 7 ] \
      || { echo "FAIL: $t quotas (imbalanced p_alt '$imb', bias $nb, null $nn; want 0.30 0.70 0.70 0.85 / 1 / 7)"; fail=1; }
    awk -F'\t' 'NR>1&&($4=="bias"||$4=="null")&&$3!=0.5{b=1} END{exit b}' $t || { echo "FAIL: $t bias/null genes must have p_alt 0.5"; fail=1; }
  done
fi
# VCF content
if [ -s $D/genome/genome.fa ]; then
  awk 'NR==FNR{ if(/^>/) next; g=g $0; next } !/^#/{ if(substr(g,$2,1)!=$4){ print "  POS " $2 " REF " $4 " genome " substr(g,$2,1); bad=1 } } END{exit bad}' \
    $D/genome/genome.fa $D/f1/parental_snps.vcf || { echo "FAIL: parental_snps.vcf REF differs from genome base"; fail=1; }
else echo "SKIP: no $D/genome/genome.fa, not checking parental VCF REF against the genome"; fi
vcf_rows() { grep -v '^#' "$1" | awk -F'\t' '{print $1"\t"$2"\t"$4"\t"$5}' | sort; }
truth_rows() { tail -n +2 "$1" | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4}' | sort; }
vcfs="$D/f1/parental_snps.vcf"
for i in $(tail -n +2 $D/outbred/samples.csv 2>/dev/null | cut -d, -f6 | sort -u); do vcfs="$vcfs $D/outbred/$i.all.vcf $D/outbred/$i.het.vcf"; done
for v in $vcfs; do [ -s $v ] && { grep -q '^#CHROM' $v || { echo "FAIL: $v has no #CHROM header line"; fail=1; }; }; done
[ -s $D/f1/parental_snps.vcf ] && { diff <(vcf_rows $D/f1/parental_snps.vcf) <(truth_rows $D/f1/truth_snps.tsv) >/dev/null \
  || { echo "FAIL: parental_snps.vcf POS/REF/ALT differ from f1/truth_snps.tsv"; fail=1; }; }
for i in $(tail -n +2 $D/outbred/samples.csv 2>/dev/null | cut -d, -f6 | sort -u); do
  a=$D/outbred/$i.all.vcf; h=$D/outbred/$i.het.vcf; [ -s $a ] && [ -s $h ] || continue
  diff <(vcf_rows $a) <(truth_rows $D/outbred/truth_snps.tsv) >/dev/null || { echo "FAIL: $a differs from outbred/truth_snps.tsv"; fail=1; }
  diff <(grep -v '^#' $a | awk -F'\t' '$10=="0/1"') <(grep -v '^#' $h) >/dev/null || { echo "FAIL: $h is not exactly the 0/1 rows of $a"; fail=1; }
done
# determinism marker
[ -s $D/seed.txt ] || { echo "FAIL: missing seed.txt"; fail=1; }
[ $fail -eq 0 ] && echo PASS || exit 1
