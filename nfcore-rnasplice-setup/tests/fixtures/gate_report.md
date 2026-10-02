# nf-core/rnasplice verification gate (Task 0)

## GATE VERDICT

GATE: B revision 1b447239488097651d8eac44bca2c1556865eb0f Nextflow 26.04.6

Release 1.0.4 does not launch under the cluster's Nextflow 26.04.6: config parse error on `def check_max` (G1).
Both fallbacks completed the nf-core test profile (G2a dev commit + 26.04.6; G2b 1.0.4 + Nextflow 24.04.4).
At the mid-gate checkpoint the user chose B (2026-10-01): dev commit `1b447239488097651d8eac44bca2c1556865eb0f`
(1.1.0dev, 2026-09-24) with Nextflow 26.04.6 in conda env `nf-env`. G3-G8 ran only for B and all passed.
- The typed params file (all modules, local data) completed: 80/80 tasks.
- `aligner: "star"` + `pseudo_aligner: "salmon"` runs DTU and SUPPA2 once.
- Explicit `false` switches modules off.
- `reverse` reaches rMATS, DEXSeq and featureCounts, from FASTQ and from genome BAM.
- Genome-BAM input works with and without the optional `strandedness`/`single_end` columns.
- rMATS bamlists follow samplesheet order; b1 = treatment.
- Every container URL answers HTTP 200.

Evidence directory: `/net/bmc-lab3/data/bcc/yannvrb/rnasplice_gate/` (`src/`, `data/`, `jobs/*.sh|*.log`, `g*.yaml`, `g*_samplesheet.csv`, `results/2026-10-01_*`, `work/`).

## Environment

- Nextflow 26.04.6 build 12646 (env `nf-env`); plugin nf-schema 2.7.2 (dev manifest `nextflowVersion = '!>=26.04.0'`, `version = '1.1.0dev'`; the summary prints `nf-core/rnasplice: v1.1.0dev-g1b44723`).
- Module lines in every job, each guarded with `|| { echo "ERROR: ..."; exit 1; }`:
  - `module add miniconda3/v4`
  - `source /home/software/conda/miniconda3/bin/condainit`
  - `conda activate nf-env`
  - `module add singularity/3.10.4`
- `NXF_*` seen in the jobs: `NXF_OFFLINE=TRUE NXF_SINGULARITY_CACHEDIR=/home/yannvrb/.singularity/cache`, inherited from the login shell's `~/.bashrc`.
- Compute-node internet: yes.
  - G0b `wget` from raw.githubusercontent.com in job 11379423.
  - The runs pulled the pipeline from GitHub.
  - The runs pulled images from depot.galaxyproject.org and community-cr-prod.seqera.io.
- Pull behaviour: `NXF_OFFLINE=TRUE` did NOT block the first download of a revision. G1: `Pulling nf-core/rnasplice:1.0.4 ... downloaded from https://github.com/nf-core/rnasplice.git`. G2a pulled the dev commit the same way. OFFLINE_PULL_NEEDED = no. Clones live in `~/.nextflow/assets/.repos/nf-core/rnasplice/clones/<commit>/`.
- Release 1.0.4 tag commit: `1d0494ae3402d1a46e0adadad24f81a0ff855c77`.
- The gate's `nextflow.config` sits in the launch dir and is also passed with `-c`, so it is loaded twice (`configFiles` lists it twice). This was harmless here.

## Runs

| Step | SLURM job | Command (all `nextflow run nf-core/rnasplice ... -c nextflow.config`) | Log | Result dir | Verdict |
|---|---|---|---|---|---|
| G0b | 11379423 | wget test data | jobs/g0_fetch.log | data/ | PASS: `G0B DONE`; 8 FASTQ of 3.66-3.83 MB, X.fa.gz 45,816,680 B, genes_chrX.gtf 25,734,901 B |
| G1 | 11379424 | `-r 1.0.4 -profile test,slurm,singularity -params-file g1.yaml` (nf-env 26.04.6) | jobs/g1.log | - | LAUNCH FAILURE (`G1A EXIT 1`, `PULL EXIT 0`, `G1B EXIT 1`) |
| G2a | 11379425 | `-r 1b44723... -profile test,slurm,singularity` (nf-env) | jobs/g2a.log | results/2026-10-01_g2a_dev | PASS: 130/130 COMPLETED, 16 min 15 s |
| G2b env | 11379426 | `conda create -n nf-24.04.4 nextflow=24.04.4` | jobs/g2b_env.log | conda env | PASS (24.04.4 build 5917) |
| G2b | 11379427 | `-r 1.0.4 -profile test,slurm,singularity` (nf-24.04.4) | jobs/g2b.log | results/2026-10-01_g2b_nf24 | PASS by rule: 128 COMPLETED + 1 FAILED with a COMPLETED retry, 56 min 45 s, after a manual kill (incident below). Not chosen. |
| G3 | 11379745 | `-r 1b44723... -profile slurm,singularity -params-file g3.yaml` | jobs/g3.log | results/2026-10-01_g3_typed | PASS: 80/80 COMPLETED, 12 min 9 s (reference run) |
| G3neg | 11379746 | g3neg.yaml (quoted numbers, bare multiqc_title, unknown `max_cpus`) | jobs/g3neg.log | - | rejected at launch (expected) |
| G3neg2 | 11379759 | g3neg2.yaml (same without `max_cpus`) | jobs/g3neg2.log | - | rejected at launch (expected) |
| G3neg3 | 11379771 | g3neg3.yaml (only quoted numbers) | jobs/g3neg3.log | - | accepted; cancelled by me after launch (own job) |
| G4 | 11379747 | g4.yaml (rMATS only, sheet `reverse`) | jobs/g4.log | results/2026-10-01_g4_rmats_reverse | PASS: 51/51 COMPLETED, 8 min 7 s |
| G4o | 11379894 | g4o.yaml (rMATS only, rows reversed within groups, `rmats_paired_stats: true`) | jobs/g4o.log | results/2026-10-01_g4o_order_paired | PASS: 51/51 COMPLETED, 9 min 9 s |
| G5 | 11379886 | g5.yaml (`source: "genome_bam"`, sheet with `strandedness` reverse, `single_end` false) | jobs/g5.log | results/2026-10-01_g5_bam | PASS: 43/43 COMPLETED, 6 min 8 s |
| G5b | 11379887 | g5b.yaml (genome_bam, sheet `sample,condition,genome_bam` only) | jobs/g5b.log | results/2026-10-01_g5b_bam_nocols | PASS: 43/43 COMPLETED, 6 min 38 s |
| G3s | 11379904 | g3s.yaml (= g3 + `save_reference: true`) `-resume d92451ab-1ce8-4159-8854-e448eff6ff8f` | jobs/g3s.log | results/2026-10-01_g3s_save_reference | PASS: 79 CACHED + 1 COMPLETED, 1 min 40 s |

G1 evidence (jobs/g1.log):
```
WARN: Cannot read project manifest -- Cause:  Config parsing failed
Error nextflow.config:368:14: Unexpected input: '('
│ 368 | def check_max(obj, type) {
ERROR ~ Config parsing failed
G1A EXIT 1
 Already-up-to-date - revision: 1d0494ae34 [1.0.4]
PULL EXIT 0
G1B EXIT 1
```

## Typing test

Full record: `typing_test.txt` (the yaml, the validator lines, the parameter summary, the negative tests).

- G3 `g3.yaml` was accepted with no WARN, validator or error line. It is the brief's G3 list with these branch-B changes:
  - added `isoformswitchanalyzer: false` and `leafcutter: false`;
  - dropped `max_cpus`, `max_memory` and `max_time` (not in the dev schema).

  The summary prints `multiqc_title : 261001`, `rmats_read_len : 75`, `aligner : star`, `sashimi_plot : false`, `isoformswitchanalyzer : false`, `diffsplice_paired : false`.
- `rmats_splice_diff_cutoff: 0.0001` is printed as `1.0E-4` and reaches rMATS as `--cstat 1.0E-4`. This is valid for rMATS; the post step completed.
- Comparison with the dev `config_params.txt`:
  - Every DEXSeq, edgeR, DTU and SUPPA value equals the config default except `diffsplice_paired`: yaml `false`, config `true`. This is the plan's deliberate deviation (planner decision 2).
  - The other differences are intended: `aligner star`, `sashimi_plot false`, `isoformswitchanalyzer false`, `rmats_read_len 75`.
  - `fasta` and `gtf` are not in the config `params {}` block. They are schema parameters and were accepted.
- Negative tests (Nextflow 26.04.6 itself checks params before nf-schema):
  - Unknown key: ``[ERROR] Parameter `max_cpus` was specified on the command line or params file but is not declared in the script or config``. This is fatal (exit 1, no task).
  - Bare number for a string: ``[ERROR] Parameter `multiqc_title` with type String cannot be assigned to 261001 [Integer]``. This is fatal.
  - Quoted number for an integer/number (`"75"`, `"0.1"`): accepted without message; the run started and was then cancelled by me.
- VALIDATOR_REJECTED = none.

## Process names, labels and requested resources

Reference run G3 trace; columns cpus / memory / time. The gate config sets the selectors, process defaults (2 CPU / 8 GB / 4 h) and `resourceLimits` 16 / 64 GB / 24 h. Values above the process default come from the pipeline's labels.

| Process (prefix `NFCORE_RNASPLICE:`) | cpus | mem | time | note |
|---|---|---|---|---|
| PREPARE_GENOME:STAR_GENOMEGENERATE | 8 | 64 GB | 8h | selector applied |
| RNASPLICE:ALIGN_STAR:STAR_ALIGN | 8 | 48 GB | 8h | selector applied |
| RNASPLICE:RMATS:RMATS_PREP | 4 | 16 GB | 8h | selector applied; dev runs it once per SAMPLE (label process_single) |
| RNASPLICE:RMATS:RMATS_POST | 8 | 32 GB | 16h | selector applied; once per contrast |
| RNASPLICE:RMATS:CREATE_BAMLIST | 1 | 6 GB | 4h | process_single |
| RNASPLICE:DEXSEQ_DEU:DEXSEQ_COUNT | 2 | 8 GB | 8h | selector applied |
| RNASPLICE:DEXSEQ_DEU:DEXSEQ_EXON | 8 | 32 GB | 8h | selector applied |
| RNASPLICE:DEXSEQ_DEU:DEXSEQ_ANNOTATION | 6 | 36 GB | 8h | label |
| RNASPLICE:DRIMSEQ_DEXSEQ_DTU_SALMON:DEXSEQ_DTU | 8 | 32 GB | 8h | selector applied |
| RNASPLICE:DRIMSEQ_DEXSEQ_DTU_SALMON:DRIMSEQ_DMFILTER / STAGER | 6 | 36 GB | 8h | label (dev name DRIMSEQ_DMFILTER; 1.0.4 had DRIMSEQ_FILTER) |
| RNASPLICE:EDGER_DEU:SUBREAD_FEATURECOUNTS | 6 | 36 GB | 8h | label |
| RNASPLICE:EDGER_DEU:EDGER_EXON | 1 | 6 GB | 4h | label |
| RNASPLICE:SALMON_QUANT_SALMON | 6 | 36 GB | 8h | `.*:SALMON_QUANT` selector does NOT match; needs `.*:SALMON_QUANT.*` |
| PREPARE_GENOME:SALMON_INDEX | 6 | 36 GB | 8h | label |
| PREPARE_GENOME:MAKE_TRANSCRIPTS_FASTA | 12 | 64 GB | 16h | label |
| RNASPLICE:SUPPA_SALMON:DIFFSPLICE_IOE/IOI, CLUSTEREVENTS_IOE/IOI | 12 | 64 GB | 16h | label |
| RNASPLICE:SUPPA_SALMON:SUPPA_PSIPEREVENT / SUPPA_PSIPERISOFORM | 6 | 36 GB | 8h | label |
| RNASPLICE:TX2GENE_TXIMPORT_SALMON:TXIMPORT | 6 | 36 GB | 8h | label |
| RNASPLICE:ALIGN_STAR:BAM_SORT_STATS_SAMTOOLS:SAMTOOLS_SORT | 6 | 36 GB | 8h | label |
| RNASPLICE:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:TRIMGALORE | 6 | 1 GB | 8h | as set by the module |

- All process names of every B run (G2a, G3, G4, G4o, G5, G5b) are in `trace_process_names.txt` (75 names).
- With the dev defaults, the test profile (G2a) additionally shows:
  - `*_STAR_SALMON` subworkflows and `SALMON_QUANT_STAR`;
  - `ISOFORMSWITCHANALYZER`;
  - `VISUALISE_MISO:MISOPY_INDEX|MISOPY_RUN|MISOPY_SETTINGS|MISOPY_SASHIMIPLOT` (renamed from 1.0.4's `MISO_*`).
- Genome-BAM runs add `BAM_SORT_STATS_SAMTOOLS:*`.
- Selector change for Task 6: `.*:SALMON_QUANT.*` instead of `.*:SALMON_QUANT`. The other seven selectors apply as written.

## Output tree

Reference run G3; `output_tree.txt` has 473 entries. Top level: `fastqc/ genome/ mergeevents/ multiqc/ salmon/ star/ trimgalore/ pipeline_info/`.

- `star/` holds the genome branch, named after `aligner`.
  - `star/<sample>_sorted.bam(.bai)`, `star/log/`, `star/samtools_stats/`.
  - rMATS:
    - `star/rmats/<contrast>/{SE,MXE,A3SS,A5SS,RI}.MATS.{JC,JCEC}.txt`;
    - `fromGTF.*.txt`, `JC(EC).raw.input.*.txt`, `summary.txt`;
    - `star/rmats/<contrast>.log`;
    - `star/rmats/bamlist/{<treatment>,<control>}_bamlist.txt`;
    - `star/rmats/prep/` (per-sample `.rmats`).
  - With `rmats_paired_stats: true` the directory is `star/rmats/<contrast>_paired/` (G2a, G4o).
  - Columns of `SE.MATS.JC.txt`: `ID GeneID geneSymbol chr strand exonStart_0base exonEnd upstreamES upstreamEE downstreamES downstreamEE ID IJC_SAMPLE_1 SJC_SAMPLE_1 IJC_SAMPLE_2 SJC_SAMPLE_2 IncFormLen SkipFormLen PValue FDR IncLevel1 IncLevel2 IncLevelDifference`.
  - DEXSeq DEU: `star/dexseq_exon/counts/<sample>.clean.count.txt`, plus `star/dexseq_exon/results/{DEXSeqDataSet,DEXSeqResults,perGeneQValue}.<contrast>.{rds,csv}` and `plotDEXSeq.<contrast>.pdf`.
  - edgeR DEU: `star/edger/contrast_<contrast>.{exprs,usage.exon,usage.gene,usage.simes}.csv` (+ pdf), `DGEList.rds`, `DGEGLM.rds`, `DGELRT.{exprs,usage}.rds`; `star/featurecounts/<sample>.featureCounts.tsv(.summary)`.
- `salmon/` holds the pseudo-alignment branch.
  - `salmon/<sample>/quant.sf`, `salmon/tximport/salmon.merged.*`.
  - DTU: `salmon/dexseq_dtu/results/dexseq/{DEXSeqResults,perGeneQValue}.<treatment>-<control>.tsv` and `salmon/dexseq_dtu/results/stager/getAdjustedPValues.<treatment>-<control>.tsv`.
  - SUPPA2: `salmon/suppa/diffsplice/per_local_event/local_<treatment>-<control>.{dpsi,psivec}`, `.../per_isoform/transcript_<treatment>-<control>.{dpsi,psivec}`, plus `clusterevents/`, `clustergroups/`, `generate_events/`, `psi_per_*`.
- NAMING FACT: DTU and SUPPA2 files are named `<treatment>-<control>` (here `GBR-YRI`), NOT the `contrast` column (`GBR_vs_YRI`). rMATS, DEXSeq DEU and edgeR use the contrast name.
- Sign conventions read from the outputs:
  - rMATS `IncLevelDifference = mean(IncLevel1) - mean(IncLevel2)` with b1 = treatment. Row 1 of SE.MATS.JC: mean(0.771, 1.0) - 1.0 = -0.115. Positive = more inclusion in treatment.
  - SUPPA dPSI: CORRECTED AFTER THE TASK 7 REVIEW (2026-10-02). The header `local_GBR-local_YRI_dPSI` (isoforms: `transcript_GBR-transcript_YRI_dPSI`) reads treatment - control, but the values are mean(PSI control) - mean(PSI treatment). Evidence (G3, results/2026-10-01_g3_typed/salmon/suppa/diffsplice/, dpsi vs psivec, awk):
    - event `ENSG00000001497;SE:X:64744930-64748140:64748249-64749092:-`: GBR (treatment) 0.977, 0.668, mean 0.823; YRI (control) 0.968, 0.986, mean 0.977; dPSI +0.154 = YRI - GBR.
    - every non-NaN event agrees: per_local_event 1098/1098 and per_isoform 1854/1854 have dPSI = mean(YRI) - mean(GBR), and none has a nonzero dPSI = mean(GBR) - mean(YRI).
    - code: subworkflows/local/suppa/main.nf:176-182 at 1b44723 passes psi1/tpm1 = treatment and psi2/tpm2 = control; SUPPA diffSplice reports cond2 - cond1 under a `cond1-cond2_dPSI` header.
    - So SUPPA2 dPSI has the same direction as the DTU fold change (positive = more in the control) and the opposite direction to rMATS IncLevelDifference. The original line ("header ..., i.e. treatment - control") inferred the sign from the header alone. Rule: the gate infers no sign from a header alone; every sign is checked on the values or the code.
  - DTU DEXSeq column: `log2fold_YRI_GBR`, i.e. control over treatment. This is the reverse direction; the downstream note must say so.
  - DEXSeq exon usage (added after the Task 7 review, header + code): star/dexseq_exon/results/DEXSeqResults.GBR_vs_YRI.csv has the column `log2fold_YRI_GBR`, and bin/run_dexseq_exon.R:34-41 at 1b44723 sets `factor(condition, levels = contrast)` with contrast = (treatment, control), so the treatment is the reference level: control over treatment, the same direction as DTU.
  - edgeR exon usage (added after the Task 7 review, code only, not checked numerically): the header of star/edger/contrast_GBR_vs_YRI.usage.exon.csv is `Geneid,Chr,Start,End,Strand,Length,logFC,exon.F,P.Value,FDR` (no direction in the name), and modules/local/edger/exon/templates/run_edger_exon.R:228-230 at 1b44723 builds `makeContrasts(paste(treatment, control, sep = "-"))` for diffSpliceDGE: logFC = treatment - control (exon usage relative to its gene).
- `genome/index/` exists but is empty with `save_reference: false`. With `save_reference: true` (G3s) it holds `genome/index/star/` (1.7 GB on chrX), `genome/index/salmon/` (760 MB), `genome/genome.transcripts.fa`, `genome/genes_chrX.tx2gene.tsv` and `genome/rsem/`.

## Strandedness and read type

Recipe: `command_lines.txt`, 96 task blocks for G3, G4, G4o, G5 and G5b.

| Run | sheet | RMATS_PREP | DEXSEQ_COUNT | SUBREAD_FEATURECOUNTS |
|---|---|---|---|---|
| G3 FASTQ | `unstranded` | `-t paired --libType fr-unstranded` | `-p yes ... -s no` | `-p ... -s 0` |
| G4 FASTQ | `reverse` | `-t paired --libType fr-firststrand` (all 4) | off | off |
| G5 genome_bam | `strandedness` reverse, `single_end` false | `-t paired --libType fr-firststrand` (4) | `-p yes ... -s reverse` (4) | `-p ... -s 2` (4) |
| G5b genome_bam | no optional columns | `-t paired --libType fr-unstranded` (4) | `-p yes ... -s no` (4) | `-p ... -s 0` (4) |

- `--readLength 75` and `--variable-read-length --allow-clipping` are on every RMATS_PREP. RMATS_POST gets `--cstat 1.0E-4` and, with `rmats_paired_stats: true`, `--paired-stats` (G4o).
- G4 module switches: no DEXSEQ, EDGER, DRIMSEQ, STAGER, SUPPA or MISO task (count 0). The Salmon pseudo branch (SALMON_INDEX, SALMON_QUANT_SALMON, TXIMPORT) still runs in an rMATS-only run because `pseudo_aligner` is `salmon`; there is no `false` value in its enum.
- Dev code paths:
  - `utils_nfcore_rnasplice_pipeline/main.nf:156` puts `strandedness` and `single_end` from the genome-BAM sheet into meta.
  - `conf/modules.config` `RMATS_PREP` ext.args maps forward->fr-secondstrand, reverse->fr-firststrand, else fr-unstranded; `-t single|paired` comes from meta.single_end.
  - `modules/local/dexseq/count/main.nf` maps forward `-s yes`, reverse `-s reverse`, unstranded `-s no`.
  - nf-core featurecounts maps 0/1/2.
  - Schema defaults (`unstranded`, `false`) are applied when the columns are absent (G5b).
- The 1.0.4 inconsistency (BAM rows get rMATS unstranded but DEXSeq count without `-s`, i.e. htseq default stranded) is NOT APPLICABLE to dev. All three tools receive the same value.
- The pipeline re-sorts and indexes input BAMs (`samtools cat | samtools sort`, then SAMTOOLS_INDEX). The sheet references only the `.bam`, so BAM_NEEDS_BAI = no. The sorted copies are published to `star/` (a disk-use note).
- Single-end BAM (`single_end: true`) was not run; it is read from the code only.

## rMATS b1 and bamlist order

- G3: `CREATE_BAMLIST` wrote `GBR_bamlist.txt: ERR188383_sorted.bam,ERR188428_sorted.bam` and `YRI_bamlist.txt: ERR188454_sorted.bam,ERR204916_sorted.bam`. RMATS_POST used `--b1 GBR_bamlist.txt --b2 YRI_bamlist.txt`. GBR = treatment of contrast `GBR_vs_YRI`, so RMATS_B1_GROUP = treatment.
- G4o: the sheet rows were reversed within each group (ERR188428, ERR188383, ERR204916, ERR188454). The bamlists were `echo ERR188428_sorted.bam ERR188383_sorted.bam > GBR_bamlist.txt` and `echo ERR204916_sorted.bam ERR188454_sorted.bam > YRI_bamlist.txt`. This is samplesheet order, not name order, so RMATS_BAMLIST_ORDER = sheet.
- The code agrees: `subworkflows/local/rmats/main.nf` indexes the samples per condition in sheet order and sorts by (condition number, index).
- Paired stats: `--paired-stats` is passed when `rmats_paired_stats: true`. The output directory becomes `<contrast>_paired`. Unequal group sizes raise an error (code).

## Salmon route

- G3 (`aligner: "star"`, `pseudo_aligner: "salmon"`): SALMON_INDEX 1, SALMON_QUANT_SALMON 4, TXIMPORT_SALMON 1, DRIMSEQ_DEXSEQ_DTU_SALMON (DRIMSEQ_DMFILTER, DEXSEQ_DTU, STAGER) 1 each, SUPPA_SALMON:* 1 each. There is no `*STAR_SALMON*` task.
- SALMON_ROUTE = pseudo_only, PSEUDO_OFF_LINE = none, so G3 is the reference run.
- With the dev defaults (`star_salmon` + `salmon`, G2a), DTU and SUPPA run twice: `*_SALMON` and `*_STAR_SALMON`. The pipeline warns "Both --aligner=star_salmon and --pseudo_aligner=salmon specified. Downstream analyses will be performed on both salmon output files."
- DTU runs one DEXSEQ_DTU task per branch, covering all contrasts. SUPPA DIFFSPLICE runs once per contrast.
- DTU_FILTER_SCOPE = all_samples: one DRIMSEQ_DMFILTER task. Its template reads the whole samplesheet (`samps[, c("sample","condition")]`) and calls `dmFilter` once.
- EDGER_DEU_FUNCTION: the dev template `modules/local/edger/exon/templates/run_edger_exon.R` uses `diffSpliceDGE`, `glmQLFit` and `glmQLFTest`.

## Containers

- `verified_urls.txt` (2026-10-01) lists 38 URLs (16 depot.galaxyproject.org + 22 community-cr-prod.seqera.io) from the dev clone's `modules/`. All answer HTTP 200.
- Splicing-module images:

| Module | Image |
|---|---|
| RMATS_PREP | `depot .../rmats:4.3.0--py311hf2f0b74_5` |
| RMATS_POST | seqera `sha256/46/4697fb...` (r-pairadise + rmats) |
| DEXSEQ_COUNT/ANNOTATION | seqera `f8/f894c1...` |
| DEXSEQ_EXON/DTU | seqera `67/673dc7...` |
| DRIMSEQ_DMFILTER | `depot .../bioconductor-drimseq:1.18.0--r40_0` |
| STAGER | seqera `41/41a73b...` |
| EDGER_EXON | `depot .../mulled-v2-419bd7f...:709335c...-0` |
| SUPPA_* | seqera `d8/d887a6...` |
| MISOPY_* | seqera `0a/0a6a58...`, `52/52399c...` |
| ISOFORMSWITCHANALYZER | `depot .../bioconductor-isoformswitchanalyzer:2.12.0--r45hd2fad28_0` |

- All of them are present in `NXF_SINGULARITY_CACHEDIR` after the runs.
- The skill must not assume galaxy-depot-only URLs: dev also pulls from community-cr-prod.seqera.io.
- Tool versions (G3 `pipeline_info/nf_core_rnasplice_software_mqc_versions.yml`): rMATS 4.3.0 (PAIRADISE 1.0), STAR 2.7.11b, Salmon 1.10.3, DEXSeq 1.56.0, DRIMSeq 1.18.0, edgeR 3.36.0, stageR 1.32.0, SUPPA 2.4, htseq 2.1.2, subread 2.1.1 (flattenGTF 2.0.1).

## Indexes

- STAR_GENOMEGENERATE (STAR 2.7.11b) wrote `versionGenome 2.7.4a`, `sjdbOverhang 100` (STAR default; the pipeline passes no `--sjdbOverhang`) and `genomeSAindexNbases 13`. The last is computed in the module as `min(14, log2(genome length)/2 - 1)`.
- SALMON_INDEX (Salmon 1.10.3) built a decoy-aware gentrome (`-d decoys.txt`) with `"indexVersion": 5`.
- Dev includes `STAR_GENOMEPARAMS_UPGRADE`, which rewrites a given `star_index`'s `versionGenome 20201` to `2.7.4a`. A pre-2.7.4a index given via `star_index` is therefore upgraded, not rejected. This is from the code, not run.
- `save_reference: true` publishes the indexes to `{outdir}/genome/index/{star,salmon}` (G3s).

## miso_genes

- The default `miso_genes` in the test GTF: ENSG00000004961 57 lines, ENSG00000005302 268, ENSG00000147403 177, so 3/3.
- In G2a (dev test profile) every MISOPY task completed, including 3 MISOPY_SASHIMIPLOT tasks.
- The skill writes `sashimi_plot: false` anyway (MISO dropped).

## Facts that change the plan

1. Pin: `-r 1b447239488097651d8eac44bca2c1556865eb0f`, VERSION_TAG `dev-1b44723`, Nextflow >= 26.04.0, tested 26.04.6, env `nf-env`, no upper bound. Affects T1, T6.
2. Global Constraints change for B:
   - `isoformswitchanalyzer` and `leafcutter` ARE parameters and must be written: `isoformswitchanalyzer: false` is mandatory, since its config default is `true`.
   - `ignore_tx_version` is a dev parameter (config default `true`); it is left at default (not emitted) in G3.
   - No `max_cpus`, `max_memory` or `max_time`: an unknown key is FATAL under 26.04.6. HAS_MAX_PARAMS = no; caps only via `process.resourceLimits`.

   Affects T5, T6 and the checker.
3. `rmats_paired_stats` config default is `false` in dev. Affects T5.
4. `multiqc_title` must be quoted: a bare number is fatal. Affects T6.
5. Selector `.*:SALMON_QUANT` must become `.*:SALMON_QUANT.*`. Affects T6.
6. dev RMATS_PREP runs per sample, with label process_single, so the 4 CPU / 16 GB selector is generous. Affects T6 (judgement).
7. BAM sheet: `sample,condition,genome_bam,strandedness,single_end` (the last two are optional in the pipeline). The skill can always write them, so all four BAM values are `from_sheet`. BAM input needs no refusal by library type. Affects T2.
8. Sample names must be valid R identifiers in dev: a letter or `.` first; letters, digits, `.`, `_`; no R keywords. This matches planner decision 8. Affects T2, T3.
9. RMATS_BAMLIST_ORDER = sheet: paired rMATS can be offered, with pairs given by row order within each condition. Affects T3, T5.
10. DTU/SUPPA output names use `<treatment>-<control>`, not the contrast name. The DTU log2fold column is control/treatment. Affects T7 (hand-off note).
11. Index reuse: STAR `versionGenome 2.7.4a` (STAR 2.7.11b), Salmon `indexVersion 5` (Salmon 1.10.3, decoy-aware). `save_reference: true` puts them in `{outdir}/genome/index/`. Affects T4.
12. `NXF_OFFLINE=TRUE` does not block the first pull. Planner decision 14 (explicit pull) is harmless but not required for a first download. Affects T6.
13. The `fasta`/`gtf` keys are schema parameters outside the config `params {}` block. A checker comparing emitted keys against `config_params.txt` alone would flag them; it must use the schema. Affects T1/T5 checker.
14. An rMATS-only run still runs the Salmon pseudo branch (index + quant), because `pseudo_aligner` has no off value. Affects T5/T7 (runtime note).

## Decisions for the user

1. **Pin**: confirm dev commit `1b447239488097651d8eac44bca2c1556865eb0f` with Nextflow 26.04.6 (env `nf-env`). It is unreleased: the skill will say so and pin the 40-character commit.
2. **BAM input**: `from_sheet`, so no library-type restriction is needed. Confirm that the skill always writes `strandedness` (the asked value) and `single_end` columns in genome-BAM sheets.
3. **Salmon route**: `pseudo_only`. With `aligner: "star"` + `pseudo_aligner: "salmon"`, DTU/SUPPA2 run once (under `salmon/`). Confirm planner decision 10. Accept that an rMATS-only run still runs Salmon quantification (it cannot be switched off).
4. **Paired rMATS**: bamlist order = samplesheet order (verified with reversed rows). Offer `rmats_paired_stats: true` for a confirmed subject pairing with equal group sizes (pairs = row order within each condition), or keep it always `false`?
5. **IsoformSwitchAnalyzeR / LeafCutter** (dev-only modules): keep them off and out of the menu (`isoformswitchanalyzer: false`, `leafcutter: false`), or offer them (spec change)?
6. **diffsplice_paired**: the plan writes `false` (config default `true`). Confirm the deviation.
7. The plan's "Planner decisions (need user review)" list. Decisions 2, 3, 10, 12 and 15 are now informed by the facts above. 15 is moot, because BAM strandedness comes from the sheet.

## Incidents

- **G2b (candidate C, not chosen).** `RMATS_PREP (GBR-YRI)` (SLURM job 11379635, work/01/8de381265523f0517bf76c60985956, rMATS 4.1.2) hung:
  - It ran 43 min at 100% CPU of one thread, wrote no output, and had `ERR188454_sorted.bam` open; the other threads waited on a futex.
  - Its twin on the same BAMs took 2 s.
  - I sent SIGTERM to the `rmats.py` process (`srun --jobid=11379635 kill -TERM 104559`, 2026-10-01 20:57). It exited 143, 1.0.4 retried it (exit 130-145, maxRetries 1), and the retry finished in 2 s.
  - Root cause not verified. Under B (rMATS 4.3.0, per-sample prep) no hang occurred in 6 runs with rMATS.
- **G3neg3.** I cancelled my own head job with `scancel 11379771` (2026-10-01 21:28:31) once the typing result was visible, to avoid a duplicate full run. No orphan task jobs remained in `squeue`; the only other task job then running (11379774) belonged to G4 (checked in `.nextflow.log.2`).
- Nothing heavy ran on the login node. It was used only for `curl` of small text/JSON, HEAD requests of container URLs, grep/awk/sed, `zcat | head` of one FASTQ, `bash -n`, `squeue`, `scontrol` and `srun --jobid` into my own task job. Nothing was deleted.

## NOT verified

- Single-end FASTQ or BAM input (`single_end: true` -> `-t single`, no `-p`).
- `forward` strandedness end to end (code only: fr-secondstrand / `-s yes` / `-s 1`).
- More than one contrast in a B run with the skill's settings (the test profile ran 2 contrasts with the defaults).
- Whether quoted numbers stay Strings inside the pipeline (run cancelled before use).
- Reuse of an existing STAR/Salmon index via `star_index`/`salmon_index`, and `STAR_GENOMEPARAMS_UPGRADE` on a real old index.
- Behaviour with a real-size genome or real data; resources are judgement.
- `gff` input, `gencode: true`, `transcriptome_bam`/`salmon_results` sources.
- Whether the rMATS 4.1.2 hang is reproducible (C only, not chosen).
- The meaning of every column of the DEXSeq/edgeR/SUPPA tables beyond the headers quoted.
- Long-term availability of the seqera community images (HTTP 200 on 2026-10-01 only).

## gate_values.tsv (as committed)

```
GATE_OUTCOME	B
PIPELINE_REVISION	1b447239488097651d8eac44bca2c1556865eb0f
VERSION_TAG	dev-1b44723
NEXTFLOW_TESTED	26.04.6
NEXTFLOW_MIN	26.04.0
NEXTFLOW_MAX_EXCL	none
CONDA_ENV_TESTED	nf-env
HAS_MAX_PARAMS	no
HAS_RESOURCE_LIMITS	yes
SALMON_ROUTE	pseudo_only
PSEUDO_OFF_LINE	none
BAM_SHEET_HEADER	sample,condition,genome_bam,strandedness,single_end
BAM_RMATS_LIBTYPE	from_sheet
BAM_RMATS_READTYPE	from_sheet
BAM_DEXSEQ_STRAND	from_sheet
BAM_FC_STRAND	from_sheet
BAM_NEEDS_BAI	no
RMATS_BAMLIST_ORDER	sheet
RMATS_B1_GROUP	treatment
STAR_VERSION_GENOME	2.7.4a
SALMON_INDEX_VERSION	5
DTU_FILTER_SCOPE	all_samples
TEST_CONTRAST	GBR_vs_YRI
TEST_READ_LENGTH	75
EDGER_DEU_FUNCTION	diffSpliceDGE glmQLFit glmQLFTest
GATE_DATE	2026-10-01
STAR_VERSION	2.7.11b
SALMON_VERSION	1.10.3
RMATS_VERSION	4.3.0
MISO_GENES_IN_TEST_GTF	3/3 (MISOPY_* all COMPLETED in G2a test profile)
OFFLINE_PULL_NEEDED	no
COMPUTE_INTERNET	yes
VALIDATOR_REJECTED	none
RUN_DIRS	/net/bmc-lab3/data/bcc/yannvrb/rnasplice_gate/results/{2026-10-01_g2a_dev,2026-10-01_g3_typed,2026-10-01_g4_rmats_reverse,2026-10-01_g4o_order_paired,2026-10-01_g5_bam,2026-10-01_g5b_bam_nocols,2026-10-01_g3s_save_reference}
```
