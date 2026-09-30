#!/bin/bash
# Usage: check_simulation_stage2.sh <outdir>   (after simulate_ase_stage2.R); exact quotas, designs and consistency
# Truth conventions: F1 p_A_* = strain-A (REF) fractions; outbred f_*/p_alt_* = ALT-side fractions (REF fraction = 1 - p_alt); see simulate_ase_stage2.R header.
set -u; D=${1:?usage: check_simulation_stage2.sh <outdir>}; fail=0
bad() { echo "FAIL: $*"; fail=1; }
chk() { [ -s "$1" ] || bad "missing/empty $1"; }
export LC_ALL=C
F=$D/f1_recip; O=$D/outbred_diff
chk $D/seed.txt; chk $D/genome/genome.fa; chk $D/genome/genome.gtf; chk $F/parental_snps.vcf; chk $F/truth_fragments.tsv; chk $O/truth_fragments.tsv
for m in $F $O; do
  chk $m/samples.csv; chk $m/truth_genes.tsv; chk $m/truth_snps.tsv
  [ "$(head -1 $m/samples.csv 2>/dev/null)" = "sample,fastq_1,fastq_2,condition,cross_direction,individual" ] || bad "$m/samples.csv header"
  for s in $(tail -n +2 $m/samples.csv 2>/dev/null | cut -d, -f1); do
    chk $m/${s}_1.fastq.gz; chk $m/${s}_2.fastq.gz; chk $m/direct_counts/$s.table
    n1=$(zcat $m/${s}_1.fastq.gz 2>/dev/null | awk 'END{print NR/4}'); n2=$(zcat $m/${s}_2.fastq.gz 2>/dev/null | awk 'END{print NR/4}')
    nf=$(awk -F'\t' -v s=$s '$1==s{t+=$3+$4} END{print t+0}' $m/truth_fragments.tsv)
    [ "$n1" = "$n2" ] && [ "$n1" = "$nf" ] || bad "$s read pairs ($n1 / $n2) differ from truth_fragments ($nf)"
    for g in $m/${s}_1.fastq.gz $m/${s}_2.fastq.gz; do
      [ "$(zcat $g 2>/dev/null | awk 'NR%4==2{print length($0)}' | sort -u)" = 100 ] || bad "$g read length not exactly 100"
    done
    awk -F'\t' 'NR==1{ if ($1!="contig"||$6!="refCount"||$7!="altCount"||$8!="totalCount") exit 1; next }
                { if ($6+$7!=$8 || $11!=$8 || $8<1) exit 1 }' $m/direct_counts/$s.table || bad "$m/direct_counts/$s.table header or ref+alt/total/rawDepth"
  done
done
# F1 design: exactly 3 samples per direction x condition; individual = sample
got=$(tail -n +2 $F/samples.csv | awk -F, '{print $5"_"$4}' | sort | uniq -c | awk '{print $2":"$1}' | paste -sd' ')
[ "$got" = "AxB_ctrl:3 AxB_treat:3 BxA_ctrl:3 BxA_treat:3" ] || bad "f1_recip design '$got'"
tail -n +2 $F/samples.csv | awk -F, '$6!=$1{e=1} END{exit e}' || bad "f1_recip: individual must equal sample"
# outbred design: 4 individuals, each exactly once in ctrl and once in treat, cross_direction NA
tail -n +2 $O/samples.csv | awk -F, '{c[$6","$4]++; i[$6]++; if ($5!="NA") e=1}
  END{ n=0; for (k in i) { n++; if (c[k",ctrl"]!=1 || c[k",treat"]!=1) e=1 } if (n!=4 || NR!=8) e=1; exit e }' || bad "outbred_diff design (4 paired individuals, 8 rows, cross_direction NA)"
# F1 truth: exact planted classes and values; cell fractions = plogis(b0 + b1 d + bc t)
keys=$(awk -F'\t' 'function pl(x){return 1/(1+exp(-x))}
  NR>1{ print $3":"sprintf("%.2f",pl($4))":"sprintf("%.2f",pl($5))":"sprintf("%.2f",pl($6))":"$11":"$12":"$13
        for (j=0;j<4;j++){ d=(j<2)?1:-1; t=j%2; p=pl($4+$5*d+$6*t); if ((p-$(7+j))^2>1e-8) print "BADCELL:"$1 } }' $F/truth_genes.tsv | sort | uniq -c | awk '{print $2"="$1}' | paste -sd' ')
want="bias:0.50:0.50:0.50:none:none:none=1 dense_null:0.50:0.50:0.50:none:none:none=2 diff_down:0.50:0.50:0.25:any:none:down=1 diff_on_strain:0.70:0.50:0.30:any:none:down=1 diff_up:0.50:0.50:0.75:any:none:up=2 maternal:0.50:0.85:0.50:none:maternal:none=1 maternal:0.50:0.95:0.50:none:maternal:none=1 null:0.50:0.50:0.50:none:none:none=44 paternal:0.50:0.10:0.50:none:paternal:none=1 paternal:0.50:0.15:0.50:none:paternal:none=1 strain_A_high:0.70:0.50:0.50:A_higher:none:none=1 strain_A_high:0.75:0.50:0.50:A_higher:none:none=1 strain_B_high:0.25:0.50:0.50:B_higher:none:none=1 strain_B_high:0.30:0.50:0.50:B_higher:none:none=1 strain_and_maternal:0.70:0.80:0.50:A_higher:maternal:none=1"
[ "$keys" = "$want" ] || bad "f1_recip truth classes/values: got '$keys'"
awk -F'\t' 'NR>1 && (($3=="dense_null" && $2!=8) || ($3=="bias" && $2!=5) || ($3!="dense_null" && $3!="bias" && $2!=4)){e=1} END{exit e}' $F/truth_genes.tsv || bad "f1_recip n_snps (8 dense, 5 bias, 4 otherwise)"
# dense_null genes: 8 SNPs within 300 transcript bases; bias gene: 5 SNPs within 60
awk -F'\t' 'NR==FNR{ if (FNR>1) c[$1]=$3; next } FNR>1{ k=$5; if (!(k in lo) || $6<lo[k]) lo[k]=$6; if (!(k in hi) || $6>hi[k]) hi[k]=$6 }
  END{ for (k in c) { if (c[k]=="dense_null" && hi[k]-lo[k]>=300) e=1; if (c[k]=="bias" && hi[k]-lo[k]>=60) e=1 } exit e }' $F/truth_genes.tsv $F/truth_snps.tsv || bad "dense/bias SNP spans"
# outbred truth: classes and values
keys=$(awk -F'\t' 'NR>1{print $3":"sprintf("%.2f",$4)":"sprintf("%.2f",$5)":"$6":"$7}' $O/truth_genes.tsv | sort | uniq -c | awk '{print $2"="$1}' | paste -sd' ')
want="base_imbalanced:0.75:0.75:FALSE:none=3 bias:0.50:0.50:FALSE:none=1 diff_consistent:0.50:0.75:TRUE:down=2 diff_phase:0.50:0.80:FALSE:mixed=3 null:0.50:0.50:FALSE:none=51"
[ "$keys" = "$want" ] || bad "outbred_diff truth classes/values: got '$keys'"
# outbred per-individual truth: p_alt follows phase; every planted gene has >= 3 individuals with a het SNP; diff_phase phases mixed
awk -F'\t' 'NR==FNR{ if (FNR>1){ c[$1]=$3; fc[$1]=$4; ft[$1]=$5; cons[$1]=$6 } next }
  FNR>1{ g=$2; if (cons[g]=="TRUE") { ec=fc[g]; et=ft[g] } else if ($3==1) { ec=fc[g]; et=ft[g] } else { ec=1-fc[g]; et=1-ft[g] }
         if ($4>0 && ((($5-ec)^2>1e-8) || (($6-et)^2>1e-8))) e=1
         if ($4>0) { nh[g]++; if ($3==1) up[g]=1; else dn[g]=1 } }
  END{ for (g in c) { if (c[g]!="null" && c[g]!="bias" && nh[g]<3) e=1; if (c[g]=="diff_phase" && !(up[g] && dn[g])) e=1 } exit e }' \
  $O/truth_genes.tsv $O/truth_individual_genes.tsv || bad "outbred_diff truth_individual_genes (phase, p_alt, >= 3 het individuals, mixed phases)"
# VCFs: REF = genome base; parental VCF = truth SNPs; het VCF = exactly the 0/1 rows of the all VCF
awk 'NR==FNR{ if(/^>/) next; g=g $0; next } !/^#/{ if(substr(g,$2,1)!=$4) e=1 } END{exit e}' $D/genome/genome.fa $F/parental_snps.vcf || bad "parental_snps.vcf REF differs from genome"
vcf_rows() { grep -v '^#' "$1" | awk -F'\t' '{print $1"\t"$2"\t"$4"\t"$5}' | sort; }
truth_rows() { tail -n +2 "$1" | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4}' | sort; }
diff <(vcf_rows $F/parental_snps.vcf) <(truth_rows $F/truth_snps.tsv) >/dev/null || bad "parental_snps.vcf differs from truth_snps.tsv"
for i in ind1 ind2 ind3 ind4; do
  chk $O/$i.all.vcf; chk $O/$i.het.vcf
  diff <(vcf_rows $O/$i.all.vcf) <(truth_rows $O/truth_snps.tsv) >/dev/null || bad "$i.all.vcf differs from truth_snps.tsv"
  diff <(grep -v '^#' $O/$i.all.vcf | awk -F'\t' '$10=="0/1"') <(grep -v '^#' $O/$i.het.vcf) >/dev/null || bad "$i.het.vcf is not the 0/1 rows of $i.all.vcf"
done
# direct counts: positions are truth SNPs (outbred: that individual's het SNPs); F1 cell fractions near the planted values
for s in $(tail -n +2 $O/samples.csv | cut -d, -f1); do
  ind=$(awk -F, -v s=$s '$1==s{print $6}' $O/samples.csv)
  awk -F'\t' 'NR==FNR{ if (!/^#/ && $10=="0/1") h[$2]=1; next } FNR>1 && !($2 in h){e=1} END{exit e}' $O/$ind.het.vcf $O/direct_counts/$s.table || bad "$s direct counts at non-het positions"
  cmp -s $O/direct_counts/$s.table $O/direct_counts/$s.unfiltered.table || bad "$s unfiltered.table must equal the table"
  [ "$(head -1 $O/direct_counts/$s.wasp_stats.tsv 2>/dev/null)" = "$(printf 'vW\tvA\tn')" ] || bad "$s wasp_stats.tsv header"
done
awk -F'\t' 'FNR==1{ f++; next } f==1{ gene[$2]=$5; next }                       # truth_snps: pos -> gene
  f==2{ split($0,a,","); if (a[1]!="sample") cell[a[1]]=a[5]"_"a[4]; next }     # samples.csv read as tab file: whole line in $0
  f==3{ p[$1"|AxB_ctrl"]=$7; p[$1"|AxB_treat"]=$8; p[$1"|BxA_ctrl"]=$9; p[$1"|BxA_treat"]=$10; next }
  { s=FILENAME; sub(/.*\//,"",s); sub(/\.table$/,"",s); k=gene[$2]"|"cell[s]; A[k]+=$6; N[k]+=$8 }
  END{ for (k in N) if (N[k]>=200 && (A[k]/N[k]-p[k])^2>0.10^2) { print "  " k, A[k]/N[k], p[k]; e=1 } exit e }' \
  $F/truth_snps.tsv $F/samples.csv $F/truth_genes.tsv $F/direct_counts/*.table || bad "F1 direct-count strain-A fractions differ from truth by more than 0.10"
# Per-sample z rules (data against truth; catches label swaps that keep the design quotas, and sign/REF-ALT errors).
# Variance of a pooled REF (or strain-A) fraction A/N at a gene in one sample, p = planted fraction:
#   p(1-p) * (PHI_BIO + 1.5/N). PHI_BIO = 0.005 is the beta-binomial (biological) dispersion; the binomial term is 1.5/N
#   because N counts a fragment once per covered SNP (mean weight ~1.3, so sum(w^2)/N^2 ~ 1.5/N for fragments covering 1-2 SNPs).
# FAIL on any |z| > 5; F1: also on sum z^2 over the 13 planted genes of a sample > ZSUM (calibrated below, see task-2-report.md).
ZSUM=${ZSUM:-45}; ZREPORT=${ZREPORT:-0}
awk -F'\t' -v zsum=$ZSUM -v rep=$ZREPORT 'FNR==1{ f++; next } f==1{ gene[$2]=$5; next }
  f==2{ split($0,a,","); cell[a[1]]=a[5]"_"a[4]; next }
  f==3{ if ($3!="null" && $3!="dense_null" && $3!="bias") pl[$1]=1; P[$1"|AxB_ctrl"]=$7; P[$1"|AxB_treat"]=$8; P[$1"|BxA_ctrl"]=$9; P[$1"|BxA_treat"]=$10; next }
  { s=FILENAME; sub(/.*\//,"",s); sub(/\.table$/,"",s); k=s"|"gene[$2]; A[k]+=$6; N[k]+=$8 }
  END{ for (k in N) { split(k,b,"|"); g=b[2]; if (!(g in pl) || N[k]<50) continue
         p=P[g"|"cell[b[1]]]; z=(A[k]/N[k]-p)/sqrt(p*(1-p)*(0.005+1.5/N[k])); if (z<0) z=-z
         ss[b[1]]+=z*z; ng[b[1]]++; if (z>mz[b[1]]) mz[b[1]]=z }
       n=0; for (s in ss) { n++; if (rep) print "F1z", s, ng[s], ss[s], mz[s]; if (mz[s]>5 || ss[s]>zsum || ng[s]<10) e=1 }
       if (n!=12) e=1; exit e }' \
  $F/truth_snps.tsv $F/samples.csv $F/truth_genes.tsv $F/direct_counts/*.table || bad "F1 per-sample z rule (|z| > 5 or sum z^2 > $ZSUM over planted genes; sample labels vs data)"
outs=$(for s in $(tail -n +2 $O/samples.csv | cut -d, -f1); do echo $O/direct_counts/$s.table; done)
awk -F'\t' -v rep=$ZREPORT 'FNR==1{ f++; next } f==1{ gene[$2]=$5; next }
  f==2{ split($0,a,","); ind[a[1]]=a[6]; cnd[a[1]]=a[4]; next }
  f==3{ pc[$1"|"$2]=$5; pt[$1"|"$2]=$6; next }
  { s=FILENAME; sub(/.*\//,"",s); sub(/\.table$/,"",s); k=s"|"gene[$2]; R[k]+=$6; N[k]+=$8 }
  END{ for (k in N) { if (N[k]<100) continue; split(k,b,"|"); q=b[1]; g=b[2]
         pa=(cnd[q]=="ctrl") ? pc[ind[q]"|"g] : pt[ind[q]"|"g]
         z=(R[k]/N[k]-(1-pa))/sqrt(pa*(1-pa)*(0.005+1.5/N[k])); if (z<0) z=-z
         nt++; if (pa!=0.5) np++; if (z>mz) mz=z; if (z>5) e=1 }
       if (rep) print "OBz tested", nt, "nonnull", np, "max|z|", mz
       if (nt<100 || np<30) e=1; exit e }' \
  $O/truth_snps.tsv $O/samples.csv $O/truth_individual_genes.tsv $outs || bad "outbred per-individual REF fraction vs 1 - p_alt (|z| > 5, or fewer than 100 tested / 30 non-null individual-gene-condition units)"
[ $fail -eq 0 ] && echo PASS || exit 1
