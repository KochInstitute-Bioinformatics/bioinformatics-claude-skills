# Test fixtures

All files in `real_testdata/` are unmodified copies of the output of one nf-core/rnasplice run on the nf-core test data
(4 paired-end human chrX samples; contrast `GBR_vs_YRI`, treatment GBR, control YRI), taken on 2026-10-03.

Run directory: `/net/bmc-lab3/data/bcc/yannvrb/rnasplice_accept4/a/`. The results directory of the run is `results/2026-10-03_a/`;
every path below is relative to it unless it says "run directory". `real_testdata/` mirrors the same relative layout.

Size rule: no file over 2 MB is stored; a larger file would be cut to its header, its first 300 rows and every row with
FDR/padj < 0.5 (SUPPA2: first 300 rows plus rows with p-val < 0.01). All files below are smaller than 2 MB, so none was trimmed.

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

## sign_rows.tsv

The header lines of the file state how its rows were chosen. In short: 3 real rows per tool (rMATS SE JC, SUPPA2 local events,
DEXSeq exon usage, DEXSeq transcript usage, edgeR exon usage) taken from the files above, with the harmonised treatment-minus-control
sign expected from the rule named in the last column. rMATS, DEXSeq exon and DTU rows have FDR/padj < 0.05, the SUPPA2 rows
p-val < 0.05. edgeR has no exon row at FDR < 0.05 in this run: two rows at FDR < 0.5 are used, plus one row with a negative logFC
that is not significant.
