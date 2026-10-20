# Test fixtures

The files in `real_testdata/` are copies of the output of one nf-core/rnasplice run on the nf-core test data (4 paired-end human
chrX samples; contrast `GBR_vs_YRI`, treatment GBR, control YRI), taken on 2026-10-03. All are unmodified except one:
`salmon/dexseq_dtu/results/dexseq/DEXSeqDataSet.GBR-YRI.rds` was reduced to its first two rows (see the table below).

Run directory: `/path/to/work/rnasplice_accept4/a/`. The results directory of the run is `results/2026-10-03_a/`;
every path below is relative to it unless it says "run directory". `real_testdata/` mirrors the same relative layout.

Size rule: no file over 2 MB is stored; a larger file would be cut to its header, its first 300 rows and every row with
FDR/padj < 0.5 (SUPPA2: first 300 rows plus rows with p-val < 0.01). All files below are smaller than 2 MB, so none was trimmed by this rule.

| Fixture (under `real_testdata/`) | Source path | Date | Trimming | Bytes |
|---|---|---|---|---|
| `star/rmats/GBR_vs_YRI/{SE,A5SS,A3SS,MXE,RI}.MATS.JC.txt` | same path | 2026-10-03 | none | 3441 to 32060 |
| `star/rmats/bamlist/{GBR,YRI}_bamlist.txt` | same path | 2026-10-03 | none | 52, 52 |
| `star/rmats/GBR_vs_YRI/SE.MATS.JCEC.txt` | same path | 2026-10-03 | none | 33665 |
| `salmon/suppa/diffsplice/per_local_event/local_GBR-YRI.dpsi` | same path | 2026-10-03 | none | 386098 |
| `salmon/suppa/diffsplice/per_isoform/transcript_GBR-YRI.dpsi` | same path | 2026-10-03 | none | 225465 |
| `star/dexseq_exon/results/DEXSeqResults.GBR_vs_YRI.csv` | same path | 2026-10-03 | none | 1736559 |
| `star/dexseq_exon/results/perGeneQValue.GBR_vs_YRI.csv` | same path | 2026-10-03 | none | 9507 |
| `star/edger/contrast_GBR_vs_YRI.usage.{exon,gene,simes}.csv` | same path | 2026-10-03 | none | 217147, 17232, 13199 |
| `salmon/dexseq_dtu/results/dexseq/DEXSeqResults.GBR-YRI.tsv` | same path | 2026-10-03 | none | 40247 |
| `salmon/dexseq_dtu/results/dexseq/perGeneQValue.GBR-YRI.tsv` | same path | 2026-10-03 | none | 1999 |
| `salmon/dexseq_dtu/results/stager/getAdjustedPValues.GBR-YRI.tsv` | same path | 2026-10-03 | none | 9114 |
| `salmon/dexseq_dtu/results/stager/stageRTx.GBR-YRI.rds` | same path | 2026-10-03 | none | 6437 |
| `salmon/dexseq_dtu/results/dexseq/DEXSeqDataSet.GBR-YRI.rds` | same path | 2026-10-03 | reduced to its first two rows (`dxd[1:2, ]`): the sample annotation (sample, condition) is complete, the counts are not | 61701 |
| `261003_a_samplesheet.csv` | `261003_a_samplesheet.csv` in the run directory | 2026-10-03 | none | 806 |
| `261003_a_contrasts.csv` | `261003_a_contrasts.csv` in the run directory | 2026-10-03 | none | 46 |

Format notes taken from these files:
- rMATS tables have the column `ID` twice (first column, and again after the coordinate columns, before `IJC_SAMPLE_1`).
- SUPPA2 files have two columns (`..._dPSI`, `..._p-val`) with the row name `GENE;EVENT_ID` as row label; untested rows are `nan`.
- The samplesheet of this run is the FASTQ schema (`sample,fastq_1,fastq_2,strandedness,condition`).

## real_hipscm/: trimmed tables of a real run

Trimmed copies of the output of one nf-core/rnasplice run on real data: human hiPSC-derived cardiomyocytes, contrast `KO_vs_WT`
(treatment KO, control WT), 3 versus 3 samples, reverse-stranded, STAR and Salmon, all five modules; finished on 2026-10-06.
Run directory: `/path/to/work/realdata_test/rnasplice_hipscm_rerun/project/`; results directory
`results/2026-10-06_project/`; every path below is relative to the results directory unless it says "run directory", and
`real_hipscm/` mirrors that layout. The source files were only read. Used by `tests/r/test_real_hipscm.R`.

Trimming. Every header line is kept byte for byte and the rows keep their order in the source file. A data row is kept when it is:
one of the first 200 data rows; a row of one of the named genes below (all rows of those genes, so that their direction can be
recomputed as on the whole run); a row whose FDR, padj or (SUPPA2) p-val is exactly 0, the first 200 such rows in file order; and,
for rMATS, one of 100 rows with low read support (the smaller group mean of IJC + SJC below 10) drawn among the rows not kept
otherwise with reservoir sampling in awk after `srand(42)`. The gene-level tables keep the rows of the genes (groups) of the kept
rows. For DTU whole genes are kept: the genes of the first 200 rows of the DEXSeq table, the named genes and the genes of the first
200 stageR rows with a transcript value of 0, with every row of those genes in the DEXSeq table, the per-gene table and the stageR table.

Named genes: RBPMS2 ENSG00000166831, FGFR2 ENSG00000066468, TFPI ENSG00000003436, PI4KB ENSG00000143393, TEAD3 ENSG00000007866
(rMATS SE 199981 and MXE 42852, both FDR 0), CLSTN1 ENSG00000171603, MTMR3 ENSG00000100330, NQO2 ENSG00000124588, ENSG00000103187,
ENSG00000186480, ENSG00000129422, ENSG00000163931, ENSG00000204681, ENSG00000006831, ENSG00000164754, ENSG00000187244,
ENSG00000197565, ENSG00000006704, ENSG00000055813 (ENST00000894230), ENSG00000091513, ENSG00000163362, ENSG00000180694,
ENSG00000183087, and the three genes of the DEXSeq group `ENSG00000282278+ENSG00000134853+ENSG00000145216` (bins with no fold
change and padj down to 1e-99). A DEXSeq group or a SUPPA2 locus `A_and_B` is kept when one of its genes is named; the SUPPA2 rows
with p-val 0 include many `_and_` loci, and the first 200 SUPPA2 rows include `nan` rows. DEXSeq exon usage and edgeR have no
row with padj or FDR exactly 0 in this run.

| Fixture (under `real_hipscm/`) | Source path | Date | Trimming | Data rows | Bytes |
|---|---|---|---|---|---|
| `star/rmats/KO_vs_WT/SE.MATS.JC.txt` | same path | 2026-10-06 | rules above, of 112465 rows | 675 | 127593 |
| `star/rmats/KO_vs_WT/A5SS.MATS.JC.txt` | same path | 2026-10-06 | rules above, of 24010 rows | 356 | 65194 |
| `star/rmats/KO_vs_WT/A3SS.MATS.JC.txt` | same path | 2026-10-06 | rules above, of 49348 rows | 465 | 86428 |
| `star/rmats/KO_vs_WT/MXE.MATS.JC.txt` | same path | 2026-10-06 | rules above, of 31068 rows | 456 | 95409 |
| `star/rmats/KO_vs_WT/RI.MATS.JC.txt` | same path | 2026-10-06 | rules above, of 15397 rows | 334 | 62927 |
| `star/rmats/bamlist/{KO,WT}_bamlist.txt` | same path | 2026-10-06 | none | 1, 1 | 48, 48 |
| `salmon/suppa/diffsplice/per_local_event/local_KO-WT.dpsi` | same path | 2026-10-06 | rules above (p-val 0), of 1271489 rows | 2374 | 291358 |
| `salmon/suppa/diffsplice/per_isoform/transcript_KO-WT.dpsi` | same path | 2026-10-06 | rules above (p-val 0), of 414380 rows | 913 | 47901 |
| `salmon/suppa/psi_per_local_event/suppa_tpm.psi` | same path | 2026-10-06 | the rows of the kept local events | 2374 | 364316 |
| `star/dexseq_exon/results/DEXSeqResults.KO_vs_WT.csv` | same path | 2026-10-06 | rules above, of 1176187 rows | 1850 | 326764 |
| `star/dexseq_exon/results/perGeneQValue.KO_vs_WT.csv` | same path | 2026-10-06 | the groups of the kept bins | 29 | 1041 |
| `star/edger/contrast_KO_vs_WT.usage.exon.csv` | same path | 2026-10-06 | rules above, of 380597 rows | 1095 | 126577 |
| `star/edger/contrast_KO_vs_WT.usage.gene.csv` | same path | 2026-10-06 | the genes of the kept exons | 111 | 9127 |
| `star/edger/contrast_KO_vs_WT.usage.simes.csv` | same path | 2026-10-06 | the genes of the kept exons | 111 | 7229 |
| `salmon/dexseq_dtu/results/dexseq/DEXSeqResults.KO-WT.tsv` | same path | 2026-10-06 | whole genes (above), of 34488 rows | 465 | 92270 |
| `salmon/dexseq_dtu/results/dexseq/perGeneQValue.KO-WT.tsv` | same path | 2026-10-06 | the same genes | 184 | 5976 |
| `salmon/dexseq_dtu/results/stager/getAdjustedPValues.KO-WT.tsv` | same path | 2026-10-06 | the same genes | 465 | 25126 |
| `salmon/dexseq_dtu/results/stager/stageRTx.KO-WT.rds` | same path | 2026-10-06 | none (the object of the whole run) | - | 1070000 |
| `261006_project_{samplesheet,contrasts}.csv`, `261006_project_params.yaml` | same names in the run directory | 2026-10-06 | none | 6, 1, - | 1294, 42, 2595 |

Total 2.8 MB. Not copied: the DTU `DEXSeqDataSet.KO-WT.rds` (35 MB), so the usage report draws the DEXSeq coefficient plot for
the top DTU genes, and the JCEC rMATS tables. The `stageRTx` object is that of the whole run: its alpha (0.05) equals the default
FDR, so at FDR 0.05 the trimmed stageR table is used as it is; at another FDR the report re-runs stageR from the object over the
whole run, which no longer matches the trimmed tables (the test uses FDR 0.05 only). SUPPA2 adjusted p-values are recomputed by the
report (Benjamini-Hochberg over the tested rows of each file), so on the trimmed files they differ from those of the whole run;
the test recomputes them from the trimmed files.

The commands that built `real_hipscm/` (bash and GNU Awk 4.0.2; the low-support rMATS rows depend on the `rand()` sequence of
that awk after `srand(42)`). Save the two blocks, one after the other, as one file, for example `build_real_hipscm.sh` in a scratch directory, and run
`bash build_real_hipscm.sh <new directory>`: it only reads the run directory and writes into `<new directory>` (and the list of DTU
genes, `dtu_genes`, next to the script). Compared with `cmp` against this directory, every file was identical when these commands were
run again on 2026-10-06.

```bash
#!/bin/bash
# Builds tests/fixtures/real_hipscm/ from the real nf-core/rnasplice output (read only). Trimming rules:
# header kept byte for byte; first 200 data rows; every row of the named genes; rows with FDR/padj (SUPPA2: p-val) exactly 0,
# the first 200 in file order; rMATS: plus 100 low-support rows (min group mean of IJC+SJC < 10) drawn with awk srand(42)
# from the rows not kept otherwise. Rows keep their file order.
set -u
RUN=/path/to/work/realdata_test/rnasplice_hipscm_rerun/project
SRC=$RUN/results/2026-10-06_project
OUT=${1:?out dir}; SCR=$(cd "$(dirname "$0")" && pwd)
GENES="ENSG00000166831 ENSG00000066468 ENSG00000003436 ENSG00000143393 ENSG00000007866 ENSG00000171603 ENSG00000100330 ENSG00000124588 ENSG00000103187 ENSG00000186480 ENSG00000129422 ENSG00000163931 ENSG00000204681 ENSG00000006831 ENSG00000164754 ENSG00000187244 ENSG00000197565 ENSG00000006704 ENSG00000055813 ENSG00000091513 ENSG00000163362 ENSG00000180694 ENSG00000183087 ENSG00000282278 ENSG00000134853 ENSG00000145216 ENSG00000166831"
mkdir -p "$OUT/star/rmats/KO_vs_WT" "$OUT/star/rmats/bamlist" "$OUT/salmon/suppa/diffsplice/per_local_event" "$OUT/salmon/suppa/diffsplice/per_isoform" \
  "$OUT/salmon/suppa/psi_per_local_event" "$OUT/star/dexseq_exon/results" "$OUT/star/edger" "$OUT/salmon/dexseq_dtu/results/dexseq" "$OUT/salmon/dexseq_dtu/results/stager"
# common awk: g[] = named genes; ingene(id, sep) true when a part of id (split on sep) is a named gene
AWKLIB='BEGIN { n = split(genes, a, " "); for (i = 1; i <= n; i++) g[a[i]] = 1 }
function ingene(id, sep,   p, k, i) { gsub(/"/, "", id); k = split(id, p, sep); for (i = 1; i <= k; i++) if (p[i] in g) return 1; return 0 }'
# rMATS: columns by header name (the ID column appears twice; GeneID, FDR and the counts once)
for t in SE A5SS A3SS MXE RI; do
  f=$SRC/star/rmats/KO_vs_WT/$t.MATS.JC.txt
  awk -F'\t' -v genes="$GENES" "$AWKLIB"'
    function gm(s,   v, k, i, t) { k = split(s, v, ","); t = 0; for (i = 1; i <= k; i++) t += v[i]; return t / k }
    NR == FNR && FNR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; srand(42); next }
    NR == FNR {
      keep = (FNR <= 201) || ingene($c["GeneID"], "+")
      if ($c["FDR"] + 0 == 0 && z < 200) { z++; keep = 1 }
      if (keep) { K[FNR] = 1; next }
      ns = gm($c["IJC_SAMPLE_1"]) + gm($c["SJC_SAMPLE_1"]); n2 = gm($c["IJC_SAMPLE_2"]) + gm($c["SJC_SAMPLE_2"]); if (n2 < ns) ns = n2
      if (ns < 10) { m++; if (m <= 100) R[m] = FNR; else { j = int(rand() * m) + 1; if (j <= 100) R[j] = FNR } }
      next }
    FNR == 1 { for (j in R) K[R[j]] = 1; print; next }
    FNR in K' "$f" "$f" > "$OUT/star/rmats/KO_vs_WT/$t.MATS.JC.txt"
done
cp "$SRC/star/rmats/bamlist/KO_bamlist.txt" "$SRC/star/rmats/bamlist/WT_bamlist.txt" "$OUT/star/rmats/bamlist/"
# SUPPA2 dPSI files: gene = text before the first ;, loci split on _and_; p-val (column 3 of the rows) exactly 0
for k in per_local_event/local_KO-WT.dpsi per_isoform/transcript_KO-WT.dpsi; do
  awk -F'\t' -v genes="$GENES" "$AWKLIB"'
    FNR == 1 { print; next }
    { gid = $1; sub(/;.*$/, "", gid); keep = (FNR <= 201) || ingene(gid, "_and_")
      if ($3 != "nan" && $3 + 0 == 0 && z < 200) { z++; keep = 1 }
      if (keep) print }' "$SRC/salmon/suppa/diffsplice/$k" > "$OUT/salmon/suppa/diffsplice/$k"
done
```

Second part of the same file (append it to the first part; it uses `SRC`, `OUT`, `SCR`, `GENES` and `AWKLIB` from the first part):

```bash
# SUPPA2 PSI per sample: the rows of the kept local events, in the order of the PSI file
awk -F'\t' 'NR == FNR { if (FNR > 1) K[$1] = 1; next } FNR == 1 || ($1 in K)' "$OUT/salmon/suppa/diffsplice/per_local_event/local_KO-WT.dpsi" \
  "$SRC/salmon/suppa/psi_per_local_event/suppa_tpm.psi" > "$OUT/salmon/suppa/psi_per_local_event/suppa_tpm.psi"
# DEXSeq exon usage bins (column 7 padj): first 200, named genes (groups split on +), padj exactly 0 (none in this run)
awk -F, -v genes="$GENES" "$AWKLIB"'
  FNR == 1 { print; next }
  { keep = (FNR <= 201) || ingene($1, "+"); if ($7 != "NA" && $7 + 0 == 0 && z < 200) { z++; keep = 1 }; if (keep) print }' \
  "$SRC/star/dexseq_exon/results/DEXSeqResults.KO_vs_WT.csv" > "$OUT/star/dexseq_exon/results/DEXSeqResults.KO_vs_WT.csv"
# per-gene q-values of the groups of the kept bins
awk -F, 'NR == FNR { if (FNR > 1) K[$1] = 1; next } FNR == 1 || ($1 in K)' "$OUT/star/dexseq_exon/results/DEXSeqResults.KO_vs_WT.csv" \
  "$SRC/star/dexseq_exon/results/perGeneQValue.KO_vs_WT.csv" > "$OUT/star/dexseq_exon/results/perGeneQValue.KO_vs_WT.csv"
# edgeR exons (column 10 FDR): first 200, named genes, FDR exactly 0 (none in this run); gene tables: the genes of the kept exons
awk -F, -v genes="$GENES" "$AWKLIB"'
  FNR == 1 { print; next }
  { keep = (FNR <= 201) || ingene($1, "+"); if ($10 + 0 == 0 && z < 200) { z++; keep = 1 }; if (keep) print }' \
  "$SRC/star/edger/contrast_KO_vs_WT.usage.exon.csv" > "$OUT/star/edger/contrast_KO_vs_WT.usage.exon.csv"
for k in gene simes; do
  awk -F, 'NR == FNR { if (FNR > 1) K[$1] = 1; next } FNR == 1 || ($1 in K)' "$OUT/star/edger/contrast_KO_vs_WT.usage.exon.csv" \
    "$SRC/star/edger/contrast_KO_vs_WT.usage.$k.csv" > "$OUT/star/edger/contrast_KO_vs_WT.usage.$k.csv"
done
# DTU: whole genes. Genes of the first 200 rows of the DEXSeq table, the named genes, and the genes of the first 200 stageR rows
# with a transcript value exactly 0; then every row of those genes in the DEXSeq table, the stageR table and the per-gene q-values
awk -F'\t' -v genes="$GENES" "$AWKLIB"'
  NR == FNR { if (FNR > 1 && $4 != "NA" && $4 + 0 == 0 && z < 200) { z++; G[$1] = 1 }; next }
  FNR == 1 { print; next }
  FNR <= 201 || ingene($1, "+") || ($1 in G) { print $1 }' "$SRC/salmon/dexseq_dtu/results/stager/getAdjustedPValues.KO-WT.tsv" \
  "$SRC/salmon/dexseq_dtu/results/dexseq/DEXSeqResults.KO-WT.tsv" | sed 1d | sort -u > "$SCR/dtu_genes"
for k in dexseq/DEXSeqResults.KO-WT.tsv dexseq/perGeneQValue.KO-WT.tsv stager/getAdjustedPValues.KO-WT.tsv; do
  awk -F'\t' 'NR == FNR { K[$1] = 1; next } FNR == 1 || ($1 in K)' "$SCR/dtu_genes" "$SRC/salmon/dexseq_dtu/results/$k" > "$OUT/salmon/dexseq_dtu/results/$k"
done
cp "$SRC/salmon/dexseq_dtu/results/stager/stageRTx.KO-WT.rds" "$OUT/salmon/dexseq_dtu/results/stager/"
cp "$RUN/261006_project_samplesheet.csv" "$RUN/261006_project_contrasts.csv" "$RUN/261006_project_params.yaml" "$OUT/"
```

## sign_rows.tsv

The header lines of the file state how its rows were chosen. The column `fixture` names the directory of the row.

`real_testdata` (15 rows): 3 real rows per tool (rMATS SE JC, SUPPA2 local events, DEXSeq exon usage, DEXSeq transcript usage,
edgeR exon usage) taken from the files above, with the harmonised treatment-minus-control sign expected from the rule named in the
column `rule`. rMATS, DEXSeq exon and DTU rows have FDR/padj < 0.05, the SUPPA2 rows p-val < 0.05. edgeR has no exon row at
FDR < 0.05 in this run: two rows at FDR < 0.5 are used, plus one row with a negative logFC that is not significant.

`real_hipscm` (34 rows, treatment KO, control WT): 10 rMATS SE events, 6 SUPPA2 local events, 6 DEXSeq exon bins, 6 DTU transcripts
and 6 edgeR exons, every one a row of the trimmed tables, all significant. The column `verification` says how the sign of each row
was checked against the data of the run: rMATS by the sample lists of the run (group 1 is KO) and IncLevel1 and IncLevel2; SUPPA2
by the mean PSI per group in `suppa_tpm.psi`; DEXSeq exon and DTU by the fitted KO and WT columns of the table; edgeR by the DEXSeq
bin with the same coordinates (the `DEXSeq.gff` of the run, in its work directory). The rows include RBPMS2 (rMATS -0.89 and -0.854,
SUPPA2 raw +0.454 harmonised -0.454, DEXSeq E023 and E024, edgeR -1.588) and one DTU transcript beyond 10 log2 units
(ENST00000958519, raw -21.54).
