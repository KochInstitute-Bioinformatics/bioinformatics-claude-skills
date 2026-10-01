#!/bin/bash
# Usage: check_simulation_stage3.sh <outdir>   (after simulate_ase_stage3.R): exact quotas, layouts, phase truth, VCFs, direct tables
set -u; D=${1:?usage: check_simulation_stage3.sh <outdir>}; fail=0
bad() { echo "FAIL: $*"; fail=1; }
chk() { [ -s "$1" ] || bad "missing/empty $1"; }
export LC_ALL=C
O=$D/outbred_phase
CH=$(grep "^>" $D/genome/genome.fa 2>/dev/null | head -1 | sed "s/^>//; s/ .*//")   # the only contig; table contigs must equal it exactly (no chr-prefix drift)
for f in $D/seed.txt $D/genome/genome.fa $D/genome/genome.gtf $O/samples.csv $O/truth_genes.tsv $O/truth_snps.tsv \
         $O/truth_individual_genes.tsv $O/truth_phase.tsv $O/truth_fragments.tsv $O/direct_flips.tsv; do chk $f; done
[ "$(head -1 $O/samples.csv 2>/dev/null)" = "sample,fastq_1,fastq_2,condition,cross_direction,individual" ] || bad "samples.csv header"
[ "$(tail -n +2 $O/samples.csv 2>/dev/null | awk -F, '{print $1":"$4":"$5":"$6}' | paste -sd' ')" = \
  "obp_s1:ctrl:NA:ind1 obp_s2:ctrl:NA:ind2 obp_s3:ctrl:NA:ind3 obp_s4:ctrl:NA:ind4 obp_s5:ctrl:NA:ind5 obp_s6:ctrl:NA:ind6" ] || bad "samples.csv design"
# classes: exact quota with h_abs, n_snps, lambda and layout
keys=$(awk -F'\t' 'NR>1{print $2":"$3":"$4":"$5":"$6}' $O/truth_genes.tsv | sort | uniq -c | awk '{print $2"="$1}' | paste -sd' ')
want="hap_lowdepth:0.4:8:30:spread=3 hap_moderate:0.17:4:450:cluster=4 hap_strong:0.3:4:450:cluster=4 null:0:4:450:random=61 null_linked:0:4:450:cluster=6 two_block:0.3:6:450:two_cluster=2"
[ "$keys" = "$want" ] || bad "truth_genes classes: got '$keys'"
# SNP layouts in transcript offsets: cluster within 300 bases; spread every 120; two clusters >= 800 apart, cluster 1 within the
# first 150 bases; random >= 40 apart
awk -F'\t' 'NR==FNR{ if (FNR>1) { lay[$1]=$6; ns[$1]=$4 } next }
  FNR>1{ g=$5; k=++n[g]; o[g,k]=$6; c[g,k]=$7 }
  END{ for (g in lay) {
         if (n[g]!=ns[g]) { print "  " g ": " n[g] " SNPs, truth says " ns[g]; e=1; continue }
         span=o[g,n[g]]-o[g,1]
         if (lay[g]=="cluster" && span>=300) { print "  " g ": cluster span " span; e=1 }
         if (lay[g]=="spread") for (i=2;i<=n[g];i++) if (o[g,i]-o[g,i-1]!=120) { print "  " g ": spread spacing"; e=1 }
         if (lay[g]=="two_cluster" && !(c[g,1]==1 && c[g,3]==1 && c[g,4]==2 && c[g,6]==2 && o[g,3]<=150 && o[g,4]-o[g,3]>=800)) { print "  " g ": two clusters"; e=1 }
         if (lay[g]!="spread") for (i=2;i<=n[g];i++) if (o[g,i]-o[g,i-1]<40) { print "  " g ": SNPs closer than 40"; e=1 }
       } exit e }' $O/truth_genes.tsv $O/truth_snps.tsv || bad "SNP layouts"
# REF = genome base
awk 'NR==FNR{ if(/^>/) next; g=g $0; next } FNR>1{ if (substr(g,$2,1)!=$3) e=1 } END{exit e}' $D/genome/genome.fa $O/truth_snps.tsv || bad "truth_snps REF differs from the genome"
# phase truth: planted and null_linked genes heterozygous at every SNP in every individual; alt_on exactly for 0/1;
# per individual at least 2 cis and 2 trans genes among the cluster genes (hap_strong, hap_moderate, null_linked)
awk -F'\t' 'NR==FNR{ if (FNR>1) cl[$1]=$2; next }
  FNR>1{ if (($5=="0/1") != ($6=="1" || $6=="2")) e=1
         if (cl[$4]!="null" && $5!="0/1") e=1
         if (cl[$4]=="hap_strong" || cl[$4]=="hap_moderate" || cl[$4]=="null_linked") { k=$1 SUBSEP $4; a[k]=a[k] $6 } }
  END{ for (k in a) { split(k, p, SUBSEP); if (a[k] ~ /1/ && a[k] ~ /2/) tr[p[1]]++; else ci[p[1]]++ }
       for (i=1;i<=6;i++) { id="ind" i; if (tr[id]<2 || ci[id]<2) e=1 } exit e }' $O/truth_genes.tsv $O/truth_phase.tsv \
  || bad "truth_phase (planted genes all het, alt_on only for 0/1, >= 2 cis and >= 2 trans cluster genes per individual)"
# VCFs against the truth
vcf_rows() { grep -v '^#' "$1" | awk -F'\t' '{print $1"\t"$2"\t"$4"\t"$5}' | sort; }
truth_rows() { tail -n +2 "$1" | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4}' | sort; }
for i in 1 2 3 4 5 6; do
  I=ind$i
  for f in all het all.phased het.phased; do chk $O/$I.$f.vcf; [ "$(grep "^#CHROM" $O/$I.$f.vcf 2>/dev/null | cut -f10)" = "$I" ] || bad "$I.$f.vcf sample column name is not $I"; done
  diff <(vcf_rows $O/$I.all.vcf) <(truth_rows $O/truth_snps.tsv) >/dev/null || bad "$I.all.vcf sites differ from truth_snps.tsv"
  diff <(vcf_rows $O/$I.all.phased.vcf) <(truth_rows $O/truth_snps.tsv) >/dev/null || bad "$I.all.phased.vcf sites differ from truth_snps.tsv"
  awk -F'\t' -v ind=$I 'NR==FNR{ if ($1==ind) g[$3]=$5; next } !/^#/{ if (g[$2]!=$10) e=1 } END{exit e}' $O/truth_phase.tsv $O/$I.all.vcf \
    || bad "$I.all.vcf GT differs from truth_phase"
  awk -F'\t' -v ind=$I 'NR==FNR{ if ($1==ind) { g[$3]=$5; a[$3]=$6 } next }
    !/^#/{ w = (g[$2]=="0/0") ? "0|0" : (g[$2]=="1/1") ? "1|1" : (a[$2]=="2") ? "0|1" : "1|0"; if ($10!=w) e=1 } END{exit e}' \
    $O/truth_phase.tsv $O/$I.all.phased.vcf || bad "$I.all.phased.vcf GT differs from truth_phase (left allele = haplotype 1)"
  diff <(grep -v '^#' $O/$I.all.vcf | awk -F'\t' '$10=="0/1"') <(grep -v '^#' $O/$I.het.vcf) >/dev/null || bad "$I.het.vcf is not the 0/1 rows of $I.all.vcf"
  diff <(grep -v '^#' $O/$I.all.phased.vcf | awk -F'\t' '$10=="0|1" || $10=="1|0"') <(grep -v '^#' $O/$I.het.phased.vcf) >/dev/null \
    || bad "$I.het.phased.vcf is not the het rows of $I.all.phased.vcf"
done
# reads, direct counts and direct gene tables per sample
for s in $(tail -n +2 $O/samples.csv | cut -d, -f1); do
  ind=$(awk -F, -v s=$s '$1==s{print $6}' $O/samples.csv)
  for f in ${s}_1.fastq.gz ${s}_2.fastq.gz direct_counts/$s.table direct_counts/$s.unfiltered.table direct_counts/$s.wasp_stats.tsv \
           direct_gene_ae_phased/$s.gene_ae.txt direct_gene_ae_unphased/$s.gene_ae.txt; do chk $O/$f; done
  n1=$(zcat $O/${s}_1.fastq.gz 2>/dev/null | awk 'END{print NR/4}'); n2=$(zcat $O/${s}_2.fastq.gz 2>/dev/null | awk 'END{print NR/4}')
  nf=$(awk -F'\t' -v s=$s '$1==s{t+=$3+$4} END{print t+0}' $O/truth_fragments.tsv)
  [ "$n1" = "$n2" ] && [ "$n1" = "$nf" ] || bad "$s read pairs ($n1 / $n2) differ from truth_fragments ($nf)"
  awk -F'\t' 'NR==1{ if ($1!="contig"||$6!="refCount"||$7!="altCount"||$8!="totalCount") exit 1; next } { if ($6+$7!=$8 || $11!=$8 || $8<1) exit 1 }' \
    $O/direct_counts/$s.table || bad "$s direct counts header or ref+alt/total/rawDepth"
  awk -F'\t' 'NR==FNR{ if (!/^#/ && $10=="0/1") h[$2]=1; next } FNR>1 && !($2 in h){e=1} END{exit e}' $O/$ind.het.vcf $O/direct_counts/$s.table \
    || bad "$s direct counts at positions that are not het in $ind"
  for k in 1 2; do
    p=$(awk -F, -v s=$s -v c=$((k+1)) '$1==s{print $c}' $O/samples.csv)
    { [ "$(basename "$p")" = "${s}_$k.fastq.gz" ] && [ -s "$O/$(basename "$p")" ]; } || bad "$s samples.csv fastq_$k path ($p) is not this sample's file"
  done
  awk -F'\t' 'FNR>1{ if (seen[$1 SUBSEP $2]++) e=1 } END{exit e}' $O/direct_counts/$s.table || bad "$s direct counts duplicate rows"
  awk -F'\t' -v c=$CH 'FNR>1 && $1!=c{e=1} END{exit e}' $O/direct_counts/$s.table || bad "$s direct counts contig is not $CH"
  for m in phased unphased; do
    t=$O/direct_gene_ae_$m/$s.gene_ae.txt
    awk -F'\t' 'NR==FNR{ if (FNR>1) g[$1]=1; next } FNR>1{ n++; if (seen[$4]++ || !($4 in g)) e=1 } END{ if (n!=80) e=1; exit e }' $O/truth_genes.tsv $t \
      || bad "$s gene_ae_$m duplicate, unknown or missing gene rows (80 distinct genes expected)"
    awk -F'\t' -v c=$CH 'FNR>1 && $1!=c{e=1} END{exit e}' $t || bad "$s gene_ae_$m contig is not $CH"
  done
  cmp -s $O/direct_counts/$s.table $O/direct_counts/$s.unfiltered.table || bad "$s unfiltered.table must equal the table"
  awk -F'\t' -v s=$s 'NR==FNR{ if ($1==s) { c1[$2]=$5; c2[$2]=$6 } next }
    FNR==1{ if ($4!="name"||$5!="aCount"||$6!="bCount"||$7!="totalCount"||$9!="n_variants"||$11!="gw_phased"||$12!="bam") e=1; next }
    { n++; if ($5!=c1[$4] || $6!=c2[$4] || $7!=$5+$6 || $11!=1 || $12!=s) e=1; if ($7==0 && $8!="inf") e=1 } END{ if (n!=80) e=1; exit e }' \
    $O/truth_fragments.tsv $O/direct_gene_ae_phased/$s.gene_ae.txt || bad "$s phased direct gene table differs from the truth cover counts"
  awk -F'\t' -v s=$s 'FILENAME ~ /truth_genes/ { if (FNR>1) cl[$1]=$2; next }
    FILENAME ~ /truth_fragments/ { if ($1==s) { c1[$2]=$5; c2[$2]=$6 } next }
    FILENAME ~ /direct_flips/ { if ($1==s) nb[$2]=$3; next }
    FNR>1 { n++; t=c1[$4]+c2[$4]; if ($7>t || $5+$6!=$7) e=1
            if (nb[$4]==1 && !(($5==c1[$4] && $6==c2[$4]) || ($5==c2[$4] && $6==c1[$4]))) e=1
            if (cl[$4]=="two_block" && $9>3) e=1
            if ($7>0 && $11!=0) e=1 } END{ if (n!=80) e=1; exit e }' \
    $O/truth_genes.tsv $O/truth_fragments.tsv $O/direct_flips.tsv $O/direct_gene_ae_unphased/$s.gene_ae.txt || bad "$s unphased direct gene table"
done
# hap_lowdepth cells: every SNP below depth 10 and at least 18 covering fragments, in every sample (6 x 3 = 18 cells);
# and the direct counts at those SNPs agree (total < 10)
awk -F'\t' 'NR==FNR{ if (FNR>1) cl[$1]=$2; next } FNR>1 && cl[$2]=="hap_lowdepth" { k++; if ($7>=10 || $5+$6<18) e=1 } END{ if (k!=18) e=1; exit e }' \
  $O/truth_genes.tsv $O/truth_fragments.tsv || bad "hap_lowdepth cells (max SNP depth < 10, >= 18 covering fragments, 18 cells)"
awk -F'\t' 'FILENAME ~ /truth_genes/ { if (FNR>1 && $2=="hap_lowdepth") lg[$1]=1; next } FILENAME ~ /truth_snps/ { if (FNR>1 && ($5 in lg)) lp[$2]=1; next }
  FNR>1 && ($2 in lp) && $8>=10 { e=1 } END{exit e}' $O/truth_genes.tsv $O/truth_snps.tsv $O/direct_counts/*.table || bad "direct counts: a hap_lowdepth SNP reaches depth 10"
# planted cells: realised haplotype-1 fractions follow h1 (z within +-5; sum of z^2 over the 78 planted cells <= 130)
awk -F'\t' 'FILENAME ~ /truth_genes/ { if (FNR>1) cl[$1]=$2; next } FILENAME ~ /truth_individual_genes/ { if (FNR>1) h[$1 SUBSEP $2]=$3; next }
  FNR>1 && cl[$2]!="null" && cl[$2]!="null_linked" { ind="ind" substr($1, 6); n=$5+$6; p=h[ind SUBSEP $2]
     if (n<1) { e=1; next } z=($5-n*p)/sqrt(n*p*(1-p)*(1+(n-1)*0.005)); k++; s2+=z*z; if (z>5 || z<-5) { print "  " $1, $2, z; e=1 } }
  END{ if (k!=78 || s2>130) { print "  cells " k ", sum z^2 " s2; e=1 } exit e }' \
  $O/truth_genes.tsv $O/truth_individual_genes.tsv $O/truth_fragments.tsv || bad "planted cells: realised haplotype-1 fractions do not follow h1"
# planted h1: |h1 - 0.5| = h_abs, both signs present per class; nulls exactly 0.5
awk -F'\t' 'NR==FNR{ if (FNR>1) { ha[$1]=$3; cl[$1]=$2 } next } FNR>1{ d=$3-0.5; a=(d<0)?-d:d; if ((a-ha[$2])^2>1e-10) e=1; if (d>0) up[cl[$2]]++; if (d<0) dn[cl[$2]]++ }
  END{ split("hap_strong hap_moderate hap_lowdepth two_block", c, " "); for (i in c) if (!up[c[i]] || !dn[c[i]]) e=1; exit e }' \
  $O/truth_genes.tsv $O/truth_individual_genes.tsv || bad "truth_individual_genes h1 (|h1 - 0.5| = h_abs, both signs in every planted class)"
# direct counts against the truth (labels): the haplotype-1 side of the counts (altCount where alt_on = 1, refCount where alt_on = 2)
# must follow h1 in the planted cells (hap_strong, hap_moderate, hap_lowdepth, two_block: 78 cells).
# y = sign(h1 - 0.5) * (H - T h1) / sqrt(T h1 (1-h1) D); D = the largest number of SNPs of the gene inside any 350-bp transcript window
# (a fragment of <= 350 bp is counted at most D times). On correct data y ~ N(0,1); a ref/alt swap turns y into a large NEGATIVE value.
# Bounds (model): a cell fails when y < -4.4 (P = 5.4e-6 per cell, 4.2e-4 per 78 cells); the dataset fails when sum y^2 > 126 = qchisq(0.9995, 78)
# (5e-4); false-fail risk per dataset < 1e-3. Observed on correct data, 5 seeds (20261001,03,05,06,07): max |y| 2.6-3.7, sum y^2 58-104
# (mean 78, as for chi-square(78)). Smallest |y| after a single-cell swap over the same seeds: hap_lowdepth 6.4, two_block 7.9,
# hap_moderate 4.5 (4.54, 4.58, 5.04, 5.32, 5.48), hap_strong 11.7: all above 4.4. STAT=1 prints the cell statistics.
ZB=4.4; ZS=126
cells=""
for s in $(tail -n +2 $O/samples.csv | cut -d, -f1); do
  ind=$(awk -F, -v s=$s '$1==s{print $6}' $O/samples.csv)
  cells="$cells
$(awk -F'\t' -v ind=$ind 'FILENAME ~ /truth_genes/ { if (FNR>1) cl[$1]=$2; next } FILENAME ~ /truth_individual_genes/ { if ($1==ind) h[$2]=$3; next }
    FILENAME ~ /truth_snps/ { if (FNR>1) { gg[$2]=$5; n[$5]++; off[$5,n[$5]]=$6 } next } FILENAME ~ /truth_phase/ { if ($1==ind) ao[$3]=$6; next }
    FNR>1 { g=gg[$2]; if (!(g in h) || cl[g]=="null" || cl[g]=="null_linked") next; T[g]+=$8; H[g]+=(ao[$2]=="1") ? $7 : $6 }
    END{ for (g in T) { D=0; for (i=1;i<=n[g];i++) { c=0; for (j=i;j<=n[g];j++) if (off[g,j]-off[g,i]<350) c++; if (c>D) D=c }
           p=h[g]; sd=sqrt(T[g]*p*(1-p)*D); sg=(p>0.5)?1:-1; print ind, g, cl[g], sg*(H[g]-T[g]*p)/sd, sg*(T[g]-H[g]-T[g]*p)/sd } }' \
    $O/truth_genes.tsv $O/truth_individual_genes.tsv $O/truth_snps.tsv $O/truth_phase.tsv $O/direct_counts/$s.table)"
done
if [ "${STAT:-0}" = 1 ]; then echo "$cells" | awk 'NF==5{print "STAT", $0}'; fi
echo "$cells" | awk -v zb=$ZB -v zs=$ZS 'NF==5{ k++; s2+=$4*$4; if ($4<-zb) { print "  " $1, $2, $4; e=1 } }
  END{ if (k!=78) { print "  cells " k; e=1 } if (s2>zs) { print "  sum z^2 " s2; e=1 } exit e }' \
  || bad "direct counts label rule: haplotype-1 side does not follow h1 (cell y < -4.4 or sum y^2 > 126; ref/alt labels vs alt_on)"
# the unphased emulation flips labels both ways (>= 10 each over all covered cells)
awk -F'\t' 'NR>1 && $5=="TRUE"{t++} NR>1 && $5=="FALSE"{f++} END{exit !(t>=10 && f>=10)}' $O/direct_flips.tsv || bad "direct_flips: fewer than 10 flipped or unflipped cells"
[ $fail -eq 0 ] && echo PASS || exit 1
