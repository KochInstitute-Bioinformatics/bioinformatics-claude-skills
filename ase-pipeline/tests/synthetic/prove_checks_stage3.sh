#!/bin/bash
# Usage: prove_checks_stage3.sh <good outdir> <scratch dir>: the check must PASS on the good data and FAIL on each injected defect
set -u; G=${1:?}; W=${2:?}; HERE=$(dirname "$0"); fail=0
bash $HERE/check_simulation_stage3.sh $G >/dev/null || { echo "FAIL: check does not pass on the good data"; exit 1; }
defect() {   # $1 = name, $2 = shell command run inside the copy ($C)
  rm -rf $W/copy; cp -r $G $W/copy; C=$W/copy/outbred_phase; eval "$2"
  if out=$(bash $HERE/check_simulation_stage3.sh $W/copy 2>&1); then echo "FAIL: defect '$1' not detected"; fail=1
  else echo "ok   defect '$1' detected: $(echo "$out" | grep "^FAIL" | cut -c1-140 | paste -sd"|")"; fi
}
defect "phased GT of one het site swapped"   'sed -i "0,/\t0|1\$/s//\t1|0/" $C/ind1.all.phased.vcf'
defect "class quota changed"                 'awk -F"\t" "BEGIN{OFS=FS} NR>1 && \$2==\"null\" && !x {\$2=\"null_linked\"; x=1} 1" $C/truth_genes.tsv > $C/t && mv $C/t $C/truth_genes.tsv'
defect "hap_lowdepth SNP depth 12"           'awk -F"\t" "BEGIN{OFS=FS} NR==FNR{ if (\$2==\"hap_lowdepth\") g[\$1]=1; next } FNR>1 && (\$2 in g) && !x {\$7=12; x=1} 1" $C/truth_genes.tsv $C/truth_fragments.tsv > $C/t && mv $C/t $C/truth_fragments.tsv'
defect "phased direct table a/b swapped"     'awk -F"\t" "BEGIN{OFS=FS} NR>1 && \$5!=\$6 && !x {t=\$5; \$5=\$6; \$6=t; x=1} 1" $C/direct_gene_ae_phased/obp_s1.gene_ae.txt > $C/t && mv $C/t $C/direct_gene_ae_phased/obp_s1.gene_ae.txt'
defect "two clusters moved together"         'awk -F"\t" "BEGIN{OFS=FS} NR==FNR{ if (\$2==\"two_block\" && !g) g=\$1; next } FNR>1 && \$5==g && \$7==2 && !x {\$6=300; x=1} 1" $C/truth_genes.tsv $C/truth_snps.tsv > $C/t && mv $C/t $C/truth_snps.tsv'
defect "one read pair dropped from mate 1"   'zcat $C/obp_s1_1.fastq.gz | tail -n +5 | gzip > $C/t && mv $C/t $C/obp_s1_1.fastq.gz'
defect "het.vcf misses a het row"            'awk "/^#/ || !x {if (!/^#/) {x=1; next}} 1" $C/ind1.het.vcf > $C/t && mv $C/t $C/ind1.het.vcf'
defect "seed.txt removed"                    'rm -f $W/copy/seed.txt'
defect "planted gene homozygous at one SNP"  'awk -F"\t" "BEGIN{OFS=FS} NR==FNR{ if (\$2==\"hap_strong\") g[\$1]=1; next } FNR>1 && (\$4 in g) && !x {\$5=\"0/0\"; \$6=\"NA\"; x=1} 1" $C/truth_genes.tsv $C/truth_phase.tsv > $C/t && mv $C/t $C/truth_phase.tsv'
defect "direct count ref+alt != total"       'awk -F"\t" "BEGIN{OFS=FS} NR==2{\$6=\$6+5} 1" $C/direct_counts/obp_s1.table > $C/t && mv $C/t $C/direct_counts/obp_s1.table'
defect "h1 flipped for one planted cell"     'awk -F"\t" "BEGIN{OFS=FS} NR==FNR{ if (\$2==\"hap_strong\") g[\$1]=1; next } FNR>1 && (\$2 in g) && !x {\$3=1-\$3; x=1} 1" $C/truth_genes.tsv $C/truth_individual_genes.tsv > $C/t && mv $C/t $C/truth_individual_genes.tsv'
defect "unphased one-block counts altered"   'awk -F"\t" "BEGIN{OFS=FS} NR==FNR{ if (\$1==\"obp_s1\" && \$3==1) g[\$2]=1; next } FNR>1 && (\$4 in g) && !x {\$5=\$5+3; \$7=\$7+3; x=1} 1" $C/direct_flips.tsv $C/direct_gene_ae_unphased/obp_s1.gene_ae.txt > $C/t && mv $C/t $C/direct_gene_ae_unphased/obp_s1.gene_ae.txt'
defect "direct counts ref/alt labels swapped (whole sample)" 'for f in $C/direct_counts/obp_s2.table $C/direct_counts/obp_s2.unfiltered.table; do awk -F"\t" "BEGIN{OFS=FS} NR>1{t=\$6; \$6=\$7; \$7=t} 1" $f > $C/t && mv $C/t $f; done'
defect "direct counts ref/alt swapped in one hap_moderate gene" 'for f in $C/direct_counts/obp_s2.table $C/direct_counts/obp_s2.unfiltered.table; do awk -F"\t" "BEGIN{OFS=FS} NR==FNR{ if (\$2==\"hap_moderate\" && !g) g=\$1; next } FILENAME ~ /truth_snps/ { if (\$5==g) gp[\$2]=1; next } FNR>1 && (\$2 in gp){t=\$6; \$6=\$7; \$7=t} 1" $C/truth_genes.tsv $C/truth_snps.tsv $f > $C/t && mv $C/t $f; done'
rm -rf $W/copy
[ $fail -eq 0 ] && echo "PROOFS PASS" || exit 1
