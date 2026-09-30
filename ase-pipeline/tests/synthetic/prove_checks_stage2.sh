#!/bin/bash
# Usage: prove_checks_stage2.sh <good outdir> <scratch dir>: the check must PASS on the good data and FAIL on each injected defect
set -u; G=${1:?}; W=${2:?}; HERE=$(dirname "$0"); fail=0
bash $HERE/check_simulation_stage2.sh $G >/dev/null || { echo "FAIL: check does not pass on the good data"; exit 1; }
defect() {   # $1 = name, $2 = shell command run inside the copy ($C)
  rm -rf $W/copy; cp -r $G $W/copy; C=$W/copy; eval "$2"
  if out=$(bash $HERE/check_simulation_stage2.sh $C 2>&1); then echo "FAIL: defect '$1' not detected"; fail=1; else echo "ok   defect '$1' detected: $(echo "$out" | grep "^FAIL" | cut -c1-140 | paste -sd"|")"; fi
}
defect "AxB sample relabelled BxA"     'sed -i "0,/,AxB,/s//,BxA,/" $C/f1_recip/samples.csv'
defect "maternal b1 set to 0"          'awk -F"\t" "BEGIN{OFS=FS} NR>1 && \$3==\"maternal\" && !x {\$5=0; x=1} 1" $C/f1_recip/truth_genes.tsv > $C/t && mv $C/t $C/f1_recip/truth_genes.tsv'
defect "REF/ALT swapped in parental VCF" 'awk -F"\t" "BEGIN{OFS=FS} !/^#/ && !x {t=\$4; \$4=\$5; \$5=t; x=1} 1" $C/f1_recip/parental_snps.vcf > $C/t && mv $C/t $C/f1_recip/parental_snps.vcf'
defect "one read pair dropped from mate 1" 's=$(tail -n +2 $C/f1_recip/samples.csv | head -1 | cut -d, -f1); zcat $C/f1_recip/${s}_1.fastq.gz | tail -n +5 | gzip > $C/t && mv $C/t $C/f1_recip/${s}_1.fastq.gz'
defect "direct count ref+alt != total"  's=$(tail -n +2 $C/f1_recip/samples.csv | head -1 | cut -d, -f1); awk -F"\t" "BEGIN{OFS=FS} NR==2{\$6=\$6+5} 1" $C/f1_recip/direct_counts/$s.table > $C/t && mv $C/t $C/f1_recip/direct_counts/$s.table'
defect "outbred individual unpaired"   'sed -i "0,/,ind1$/s//,ind9/" $C/outbred_diff/samples.csv'
defect "cell fraction inconsistent"     'awk -F"\t" "BEGIN{OFS=FS} NR>1 && \$3==\"strain_A_high\" && !x {\$7=0.5; x=1} 1" $C/f1_recip/truth_genes.tsv > $C/t && mv $C/t $C/f1_recip/truth_genes.tsv'
defect "diff_phase phases not mixed"   'awk -F"\t" "BEGIN{OFS=FS} NR==FNR{ if (\$3==\"diff_phase\") g[\$1]=1; next } FNR>1 && (\$2 in g) {\$3=1} 1" $C/outbred_diff/truth_genes.tsv $C/outbred_diff/truth_individual_genes.tsv > $C/t && mv $C/t $C/outbred_diff/truth_individual_genes.tsv'
defect "direct F1 fractions altered"    'for f in $C/f1_recip/direct_counts/*.table; do awk -F"\t" "BEGIN{OFS=FS} NR>1{\$6=int(\$8*0.9); \$7=\$8-\$6} 1" $f > $C/t && mv $C/t $f; done'
defect "seed.txt removed"              'rm -f $C/seed.txt'
# Stage 2 fix round: data-vs-truth rules (sample labels, outbred sign)
defect "F1 condition labels swapped within a direction (s01 ctrl <-> s04 treat)" 'awk -F, "BEGIN{OFS=FS} NR>1 && (\$1==\"f1r_s01\" || \$1==\"f1r_s04\"){\$4=(\$4==\"ctrl\")?\"treat\":\"ctrl\"} 1" $C/f1_recip/samples.csv > $C/t && mv $C/t $C/f1_recip/samples.csv'
defect "F1 cross_direction swapped between two samples (s01 AxB <-> s07 BxA)" 'awk -F, "BEGIN{OFS=FS} NR>1 && (\$1==\"f1r_s01\" || \$1==\"f1r_s07\"){\$5=(\$5==\"AxB\")?\"BxA\":\"AxB\"} 1" $C/f1_recip/samples.csv > $C/t && mv $C/t $C/f1_recip/samples.csv'
defect "outbred REF/ALT counts swapped in one table" 'awk -F"\t" "BEGIN{OFS=FS} NR>1{t=\$6; \$6=\$7; \$7=t} 1" $C/outbred_diff/direct_counts/ob_s1.table > $C/t && mv $C/t $C/outbred_diff/direct_counts/ob_s1.table; cp $C/outbred_diff/direct_counts/ob_s1.table $C/outbred_diff/direct_counts/ob_s1.unfiltered.table'
defect "outbred p_alt read as REF fraction (truth flipped p -> 1-p consistently)" 'awk -F"\t" "BEGIN{OFS=FS} NR>1{\$4=1-\$4; \$5=1-\$5} 1" $C/outbred_diff/truth_genes.tsv > $C/t && mv $C/t $C/outbred_diff/truth_genes.tsv; awk -F"\t" "BEGIN{OFS=FS} NR>1 && \$4>0{\$5=1-\$5; \$6=1-\$6} 1" $C/outbred_diff/truth_individual_genes.tsv > $C/t && mv $C/t $C/outbred_diff/truth_individual_genes.tsv'
defect "outbred diff_consistent expect_diff up (plan convention: down)" 'sed -i "s/\tdown\$/\tup/; " $C/outbred_diff/truth_genes.tsv'
rm -rf $W/copy
[ $fail -eq 0 ] && echo "PROOFS PASS" || exit 1
