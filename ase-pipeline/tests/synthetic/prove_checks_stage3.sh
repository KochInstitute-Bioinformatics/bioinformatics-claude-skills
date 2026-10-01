#!/bin/bash
# Usage: prove_checks_stage3.sh <good outdir> <scratch dir>   (optional env PHASE_RUN="singularity exec ... Rscript" also proves the read-level check)
# The check must PASS on the good data; every injected defect must make the SPECIFIC check named by its tag FAIL (the tag is a fixed
# substring of the check's FAIL line). A defect fails the proof when its tag is absent, i.e. when no check or only another check fires.
# Other checks may fire in addition (a changed truth file also breaks the VCF comparison, for example); only the tag is asserted.
set -u; G=${1:?}; W=${2:?}; HERE=$(dirname "$0"); fail=0
bash $HERE/check_simulation_stage3.sh $G >/dev/null || { echo "FAIL: check does not pass on the good data"; exit 1; }
defect() {   # $1 = name, $2 = shell command run inside the copy ($C), $3 = tag that the check must print in a FAIL line
  rm -rf $W/copy; cp -r $G $W/copy; C=$W/copy/outbred_phase; eval "$2"
  out=$(bash $HERE/check_simulation_stage3.sh $W/copy 2>&1)
  if echo "$out" | grep "^FAIL" | grep -qF -- "$3"; then echo "ok   defect '$1' fired '$3': $(echo "$out" | grep "^FAIL" | cut -c1-150 | paste -sd"|")"
  else echo "FAIL: defect '$1': expected check '$3' did not fire; got: $(echo "$out" | grep "^FAIL" | cut -c1-150 | paste -sd"|")"; fail=1; fi
}
# single-cell ref/alt swap in the direct counts of the planted cell of class $1 with the SMALLEST swapped |z| (worst case)
swapcell() {
  read ind gene zz < <(STAT=1 bash $HERE/check_simulation_stage3.sh $G | awk -v c=$1 '$1=="STAT" && $4==c { a=($6<0)?-$6:$6; if (m=="" || a<m) { m=a; i=$2; g=$3 } } END{print i, g, m}')
  s=obp_s${ind#ind}; echo "     (class $1: weakest cell $ind $gene, |z| after swap $zz)"
  for f in $C/direct_counts/$s.table $C/direct_counts/$s.unfiltered.table; do
    awk -F"\t" -v g=$gene 'BEGIN{OFS=FS} NR==FNR{ if ($5==g) gp[$2]=1; next } FNR>1 && ($2 in gp){t=$6; $6=$7; $7=t} 1' $C/truth_snps.tsv $f > $C/t && mv $C/t $f
  done
}
defect "phased GT of one het site swapped"   'sed -i "0,/\t0|1\$/s//\t1|0/" $C/ind1.all.phased.vcf' "ind1.all.phased.vcf GT differs from truth_phase"
defect "class quota changed"                 'awk -F"\t" "BEGIN{OFS=FS} NR>1 && \$2==\"null\" && !x {\$2=\"null_linked\"; x=1} 1" $C/truth_genes.tsv > $C/t && mv $C/t $C/truth_genes.tsv' "truth_genes classes"
defect "hap_lowdepth SNP depth 12"           'awk -F"\t" "BEGIN{OFS=FS} NR==FNR{ if (\$2==\"hap_lowdepth\") g[\$1]=1; next } FNR>1 && (\$2 in g) && !x {\$7=12; x=1} 1" $C/truth_genes.tsv $C/truth_fragments.tsv > $C/t && mv $C/t $C/truth_fragments.tsv' "hap_lowdepth cells"
defect "phased direct table a/b swapped"     'awk -F"\t" "BEGIN{OFS=FS} NR>1 && \$5!=\$6 && !x {t=\$5; \$5=\$6; \$6=t; x=1} 1" $C/direct_gene_ae_phased/obp_s1.gene_ae.txt > $C/t && mv $C/t $C/direct_gene_ae_phased/obp_s1.gene_ae.txt' "obp_s1 phased direct gene table differs from the truth cover counts"
defect "two clusters moved together"         'awk -F"\t" "BEGIN{OFS=FS} NR==FNR{ if (\$2==\"two_block\" && !g) g=\$1; next } FNR>1 && \$5==g && \$7==2 && !x {\$6=300; x=1} 1" $C/truth_genes.tsv $C/truth_snps.tsv > $C/t && mv $C/t $C/truth_snps.tsv' "SNP layouts"
defect "one read pair dropped from mate 1"   'zcat $C/obp_s1_1.fastq.gz | tail -n +5 | gzip > $C/t && mv $C/t $C/obp_s1_1.fastq.gz' "obp_s1 read pairs"
defect "het.vcf misses a het row"            'awk "/^#/ || !x {if (!/^#/) {x=1; next}} 1" $C/ind1.het.vcf > $C/t && mv $C/t $C/ind1.het.vcf' "ind1.het.vcf is not the 0/1 rows of ind1.all.vcf"
defect "seed.txt removed"                    'rm -f $W/copy/seed.txt' "seed.txt"
defect "planted gene homozygous at one SNP"  'awk -F"\t" "BEGIN{OFS=FS} NR==FNR{ if (\$2==\"hap_strong\") g[\$1]=1; next } FNR>1 && (\$4 in g) && !x {\$5=\"0/0\"; \$6=\"NA\"; x=1} 1" $C/truth_genes.tsv $C/truth_phase.tsv > $C/t && mv $C/t $C/truth_phase.tsv' "truth_phase (planted genes all het"
defect "direct count ref+alt != total"       'awk -F"\t" "BEGIN{OFS=FS} NR==2{\$6=\$6+5} 1" $C/direct_counts/obp_s1.table > $C/t && mv $C/t $C/direct_counts/obp_s1.table' "obp_s1 direct counts header or ref+alt/total/rawDepth"
defect "h1 flipped for one planted cell"     'awk -F"\t" "BEGIN{OFS=FS} NR==FNR{ if (\$2==\"hap_strong\") g[\$1]=1; next } FNR>1 && (\$2 in g) && !x {\$3=1-\$3; x=1} 1" $C/truth_genes.tsv $C/truth_individual_genes.tsv > $C/t && mv $C/t $C/truth_individual_genes.tsv' "planted cells: realised haplotype-1 fractions do not follow h1"
defect "unphased one-block counts altered"   'awk -F"\t" "BEGIN{OFS=FS} NR==FNR{ if (\$1==\"obp_s1\" && \$3==1) g[\$2]=1; next } FNR>1 && (\$4 in g) && !x {\$5=\$5+3; \$7=\$7+3; x=1} 1" $C/direct_flips.tsv $C/direct_gene_ae_unphased/obp_s1.gene_ae.txt > $C/t && mv $C/t $C/direct_gene_ae_unphased/obp_s1.gene_ae.txt' "obp_s1 unphased direct gene table"
defect "direct counts ref/alt labels swapped (whole sample)" 'for f in $C/direct_counts/obp_s2.table $C/direct_counts/obp_s2.unfiltered.table; do awk -F"\t" "BEGIN{OFS=FS} NR>1{t=\$6; \$6=\$7; \$7=t} 1" $f > $C/t && mv $C/t $f; done' "direct counts label rule"
defect "ref/alt swapped in ONE hap_lowdepth cell"  'swapcell hap_lowdepth' "direct counts label rule"
defect "ref/alt swapped in ONE hap_moderate cell"  'swapcell hap_moderate' "direct counts label rule"
defect "ref/alt swapped in ONE two_block cell"     'swapcell two_block' "direct counts label rule"
defect "gene row duplicated in phased gene_ae (replaces another gene)" 'awk -F"\t" "BEGIN{OFS=FS} NR==3{\$4=prev} {if (NR==2) prev=\$4} 1" $C/direct_gene_ae_phased/obp_s1.gene_ae.txt > $C/t && mv $C/t $C/direct_gene_ae_phased/obp_s1.gene_ae.txt' "obp_s1 gene_ae_phased duplicate, unknown or missing gene rows"
defect "direct-count rows duplicated"        'for f in $C/direct_counts/obp_s1.table $C/direct_counts/obp_s1.unfiltered.table; do l=$(sed -n 5p $f); echo "$l" >> $f; done' "obp_s1 direct counts duplicate rows"
defect "chr2 contig in a direct-count row"   'for f in $C/direct_counts/obp_s1.table $C/direct_counts/obp_s1.unfiltered.table; do awk -F"\t" "BEGIN{OFS=FS} NR==3{\$1=\"chr2\"} 1" $f > $C/t && mv $C/t $f; done' "obp_s1 direct counts contig is not"
defect "chr2 contig in a phased gene_ae row" 'awk -F"\t" "BEGIN{OFS=FS} NR==3{\$1=\"chr2\"} 1" $C/direct_gene_ae_phased/obp_s1.gene_ae.txt > $C/t && mv $C/t $C/direct_gene_ae_phased/obp_s1.gene_ae.txt' "obp_s1 gene_ae_phased contig is not"
defect "VCF sample column renamed"           'sed -i "s/^\(#CHROM.*\)ind1\$/\1ind2/" $C/ind1.all.phased.vcf' "ind1.all.phased.vcf sample column name"
defect "samples.csv points obp_s1 mate 1 at another sample's file" 'sed -i "s#obp_s1_1.fastq.gz#obp_s2_1.fastq.gz#" $C/samples.csv' "obp_s1 samples.csv fastq_1 path"
rm -rf $W/copy

# read-level check (needs R: PHASE_RUN="singularity exec --bind ... sif Rscript")
pdefect() {   # $1 = name, $2 = truth_phase awk program body edit, $3 = tag required in the output
  rm -rf $W/pcopy; mkdir -p $W/pcopy/outbred_phase; cp -r $G/genome $W/pcopy/
  for f in $G/outbred_phase/*; do ln -s $f $W/pcopy/outbred_phase/; done; rm $W/pcopy/outbred_phase/truth_phase.tsv
  awk -F"\t" "BEGIN{OFS=FS} $2 1" $G/outbred_phase/truth_phase.tsv > $W/pcopy/outbred_phase/truth_phase.tsv
  out=$($PHASE_RUN $HERE/check_phase_reads_stage3.R $W/pcopy 2>&1)
  if echo "$out" | grep -qF -- "$3" && echo "$out" | grep -q "PHASE READS FAIL"; then echo "ok   read defect '$1' fired '$3': $(echo "$out" | grep -F FAIL | cut -c1-170 | paste -sd'|')"
  else echo "FAIL: read defect '$1': expected '$3' did not fire; got: $(echo "$out" | tail -8 | cut -c1-170 | paste -sd'|')"; fail=1; fi
}
if [ -n "${PHASE_RUN:-}" ]; then
  g=$(awk -F'\t' '$2=="hap_strong"{print $1; exit}' $G/outbred_phase/truth_genes.tsv)
  pdefect "alt_on swapped at one SNP of $g in ind1 (partition broken)" '$1=="ind1" && $4=="'$g'" && !x {$6=($6=="1")?"2":"1"; x=1}' "partition FAIL"
  pdefect "haplotype 1/2 swapped everywhere in ind1 (partition intact, orientation reversed)" '$1=="ind1" && $6!="NA" {$6=($6=="1")?"2":"1"}' "orientation FAIL"
  pdefect "haplotype 1/2 swapped everywhere in ind3" '$1=="ind3" && $6!="NA" {$6=($6=="1")?"2":"1"}' "orientation FAIL"
  rm -rf $W/pcopy
fi
[ $fail -eq 0 ] && echo "PROOFS PASS" || exit 1
