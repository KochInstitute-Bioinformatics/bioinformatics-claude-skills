# nf-core/rnasplice verification gate (Task 0) — INTERIM, mid-gate checkpoint

## GATE VERDICT

GATE: PENDING (A failed at launch; B and C both ran the test profile to completion) — awaiting the user's choice at the mid-gate checkpoint

Release 1.0.4 does not launch under the cluster's Nextflow 26.04.6: the strict config parser rejects the
`def check_max(obj, type)` function in the 1.0.4 `nextflow.config` (G1, both before and after an explicit pull).
Both fallback candidates of Step G2 ran the nf-core `test` profile (all modules on, including MISO) to
`Pipeline completed successfully`: B = dev commit `1b447239488097651d8eac44bca2c1556865eb0f` with Nextflow 26.04.6
(130/130 tasks COMPLETED), C = release 1.0.4 with Nextflow 24.04.4 in a new conda env `nf-24.04.4`
(128 COMPLETED + 1 FAILED whose retry completed; the failure was an rMATS 4.1.2 hang that I killed by hand,
see "Incident / rMATS hang"). Steps G3-G9 are not run yet: per the plan they use only the candidate the user chooses.
Gate evidence directory: `/net/bmc-lab3/data/bcc/yannvrb/rnasplice_gate/` (src/, data/, jobs/, results/, work/).

## Environment

- Login shell: `NXF_OFFLINE=TRUE`, `NXF_SINGULARITY_CACHEDIR=/home/yannvrb/.singularity/cache` (from `~/.bashrc`); both inherited by the jobs (G1 log line `NXF_OFFLINE=TRUE NXF_SINGULARITY_CACHEDIR=/home/yannvrb/.singularity/cache`).
- Module lines used (all guarded): `module add miniconda3/v4`, `source /home/software/conda/miniconda3/bin/condainit`, `conda activate <env>`, `module add singularity/3.10.4`.
- Nextflow: `nf-env` = 26.04.6 build 12646; `nf-24.04.4` (created by job 11379426 from conda-forge/bioconda `nextflow=24.04.4`) = 24.04.4 build 5917.
- Compute-node internet: yes (G0b `wget` from raw.githubusercontent.com on a compute node, job 11379423; the runs pulled the pipeline from GitHub and images from depot.galaxyproject.org and community-cr-prod.seqera.io).
- Pull behaviour: `NXF_OFFLINE=TRUE` did NOT block the first download. G1 log: `Pulling nf-core/rnasplice:1.0.4 ... downloaded from https://github.com/nf-core/rnasplice.git`. So OFFLINE_PULL_NEEDED = no (the G1 failure that followed is a parse error, not an offline error). Under 26.04.6 the clone sits in `~/.nextflow/assets/.repos/nf-core/rnasplice/clones/<commit>/`.
- 1.0.4 tag commit: `1d0494ae3402d1a46e0adadad24f81a0ff855c77` (`src/1.0.4/tag_ref.json`). Manifest `nextflowVersion = '!>=23.04.0'`, plugin `nf-validation@1.1.3`.
- dev commit tested: `1b447239488097651d8eac44bca2c1556865eb0f` (2026-09-24, "Merge pull request #291 from piplus2/rmats-subworkflow-nfcore"); manifest `version = '1.1.0dev'`, `nextflowVersion = '!>=26.04.0'`, plugin `nf-schema@2.7.2`.
- Note: the launch directory holds the gate's `nextflow.config`, so Nextflow loads it twice (auto-loaded from the launch dir and via `-c`); the run header lists it twice under `configFiles`. Harmless here.

## Runs

| Step | Job id | Command | Log | Result dir | Verdict |
|---|---|---|---|---|---|
| G0b | 11379423 | wget test data on a compute node | jobs/g0_fetch.log | data/ | PASS (`G0B DONE`; 8 FASTQ 3.66-3.83 MB, X.fa.gz 45,816,680 B, genes_chrX.gtf 25,734,901 B) |
| G1 | 11379424 | `nextflow run nf-core/rnasplice -r 1.0.4 -profile test,slurm,singularity -c nextflow.config -params-file g1.yaml` (nf-env, 26.04.6) | jobs/g1.log | none | LAUNCH FAILURE (`G1A EXIT 1`, `PULL EXIT 0`, `G1B EXIT 1`) |
| G2a | 11379425 | same with `-r 1b447239488097651d8eac44bca2c1556865eb0f` (nf-env, 26.04.6) | jobs/g2a.log | results/2026-10-01_g2a_dev | PASS: `G2A-A EXIT 0`, `Pipeline completed successfully`, 130 COMPLETED, 0 FAILED, 57 min |
| G2b env | 11379426 | `conda create -n nf-24.04.4 nextflow=24.04.4` | jobs/g2b_env.log | ~/.conda env | PASS (`version 24.04.4 build 5917`) |
| G2b | 11379427 | `-r 1.0.4` in env nf-24.04.4 | jobs/g2b.log | results/2026-10-01_g2b_nf24 | PASS by the brief's rule (128 COMPLETED, 1 FAILED + retry COMPLETED), `G2B-A EXIT 0`, 56 min 44 s; only after a manual kill of a hung task |

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
This is a launch failure (config compile error before the first task): 1.0.4 is incompatible with Nextflow 26.04.6.

## Pin candidates

### B — dev commit `1b44723` with Nextflow 26.04.6 (env nf-env)
- Test profile: 130/130 COMPLETED (job 11379425). All modules incl. IsoformSwitchAnalyzeR (on by default in dev) and MISO (MISOPY_* all COMPLETED, 3 sashimi plots).
- Tool versions (pipeline_info/nf_core_rnasplice_software_mqc_versions.yml): rMATS 4.3.0, STAR 2.7.11b (+2.7.10a), Salmon 1.10.3, DEXSeq 1.56.0, stageR 1.32.0, SUPPA 2.4, htseq 2.1.2, misopy 0.5.4.
- Schema differences that change the plan: no `max_cpus`/`max_memory`/`max_time` (resourceLimits only); new params `isoformswitchanalyzer` (config default TRUE: must be written `false` explicitly), `leafcutter` (false), `ignore_tx_version` (config default true); `rmats_paired_stats` default false; FASTQ sheet header unchanged `sample,fastq_1,fastq_2,strandedness,condition` (sample names must now be valid R identifiers); genome-BAM sheet (`assets/schema_input_genome_bam.json`) = `sample,condition,genome_bam` required + optional `strandedness` (enum forward/reverse/unstranded, default unstranded) and `single_end` (default false).
- rMATS code: per-sample `RMATS_PREP` (nf-core module), then `CREATE_BAMLIST` + `RMATS_POST` per contrast; bamlists built in explicit samplesheet order (`ids1.withIndex()` ... sort by condition number and sheet index), b1 = treatment. So B gives RMATS_BAMLIST_ORDER = sheet (paired rMATS could be offered).
- Risks: unreleased moving branch (the rMATS subworkflow was rewritten one week before this gate, PR #291, 2026-09-24); no release tag; the skill pins a 40-character commit.
- Output layout seen: `star_salmon/rmats/<contrast>_paired/{A3SS,...}.MATS.JC(EC).txt`, `star_salmon/rmats/bamlist/`, plus `isoformswitchanalyzer/`, `mergeevents/`, `salmon/`, `star_salmon/{dexseq_dtu,dexseq_exon,edger,featurecounts,suppa,tximport,...}`.

### C — release 1.0.4 with Nextflow 24.04.4 (new env nf-24.04.4)
- Test profile: 128 COMPLETED + 1 FAILED (`RMATS_PREP (GBR-YRI)`, exit 143) whose retry completed in 2 s (job 11379427). MISO_* all COMPLETED.
- Tool versions (pipeline_info/software_versions.yml): rMATS 4.1.2, STAR 2.7.9a (+2.7.10a), Salmon 1.10.1, DEXSeq 1.36.0, htseq 2.0.2, misopy 0.5.4.
- Needs a second conda env (the skill would have to create/check `nf-24.04.4` and refuse Nextflow >= 26.04: NEXTFLOW_MAX_EXCL = 26.04.0). 24.04.4 also supports `process.resourceLimits`.
- rMATS code (1.0.4 `subworkflows/local/rmats.nf`, `workflows/rnasplice.nf` lines 397-413): BAMs grouped per condition with `groupTuple(by:0)` (channel arrival order, no sort) -> RMATS_BAMLIST_ORDER = unordered (paired rMATS must not be offered); `--b1` = `${cond1}_bamlist.txt` with cond1 = contrast `treatment` -> b1 = treatment.
- Genome-BAM input in 1.0.4 (from code, not yet run): `create_genome_bam_channel` sets no `strandedness`/`single_end` -> RMATS_PREP `-t paired --libType fr-unstranded`; DEXSEQ_COUNT passes no `-s` (htseq dexseq_count.py default is stranded `yes`) and `-p yes`; featureCounts strand from the nf-core module default (not read). To be confirmed by G5 if C is chosen.

### Incident / rMATS hang (candidate C)
`RMATS_PREP (GBR-YRI)` (SLURM job 11379635, work/01/8de381265523f0517bf76c60985956) ran 43 min at 100% CPU of one thread with 40 MB RSS and wrote nothing (`rmats_temp/` empty, `.command.trace` empty); the other threads waited on a futex; fd 5 = `ERR188454_sorted.bam`. Its twin `RMATS_PREP (YRI-GBR)` on the same four BAMs finished in 2 s. I sent SIGTERM to the `rmats.py` process (pid 104559 via `srun --jobid=11379635 kill -TERM`, 2026-10-01 20:57), the task exited 143, 1.0.4's retry rule (exit 130-145, maxRetries 1) resubmitted it and the retry completed in 2 s. Without the kill the task would have held its slot until its 8 h limit (whether that would then have been retried is NOT VERIFIED). Interpretation: a non-deterministic hang in rMATS 4.1.2 prep with `--nthread 4` (not reproduced, root cause not verified). Candidate B (rMATS 4.3.0, per-sample prep) showed no hang in this single run.

### Other facts already measured (both candidates)
- Selector check (G2a trace, requested cpus/mem/time): STAR_GENOMEGENERATE 8/64 GB/8h, STAR_ALIGN 8/48 GB/8h, RMATS_PREP 4/16 GB/8h, RMATS_POST 8/32 GB/16h, DEXSEQ_COUNT 2/8 GB/8h, DEXSEQ_EXON 8/32 GB/8h, DEXSEQ_DTU 8/32 GB/8h — applied. The `.*:SALMON_QUANT` selector does NOT match: the real names are `SALMON_QUANT_SALMON` and `SALMON_QUANT_STAR` (they ran at the label default 6/36 GB/8h). Task 6 must use `.*:SALMON_QUANT.*` (or both names).
- With the defaults (`aligner star_salmon`, `pseudo_aligner salmon`) both Salmon branches run DTU and SUPPA twice (`DRIMSEQ_DEXSEQ_DTU_SALMON` and `..._STAR_SALMON`, `SUPPA_SALMON` and `SUPPA_STAR_SALMON`), with the WARN "Downstream analyses will be performed on both salmon output files." DEXSEQ_DTU runs once per branch (not once per contrast).
- Default `miso_genes` in the test GTF: ENSG00000004961 57 lines, ENSG00000005302 268, ENSG00000147403 177 -> 3/3; MISO tasks completed in both candidates.
- Test FASTQ read length: 75 (1000/1000 reads of `ERR188383_chrX_1.fastq.gz`; the brief's one-liner names `ERR188383_1.fastq.gz`, the real file is `ERR188383_chrX_1.fastq.gz`).
- 1.0.4 static code (src/1.0.4): DRIMSeq filter runs once on all samples of the samplesheet (`run_drimseq_filter.R`, one DRIMSEQ_FILTER call) -> DTU_FILTER_SCOPE all_samples; edgeR exon script uses `glmQLFit glmQLFTest diffSpliceDGE`.
- Test samplesheet header `sample,fastq_1,fastq_2,strandedness,condition`; contrasts `GBR-YRI,GBR,YRI` and `YRI-GBR,YRI,GBR`.

## Decisions for the user (mid-gate checkpoint)

1. Pinned dev commit `1b447239488097651d8eac44bca2c1556865eb0f` with Nextflow 26.04.6 in `nf-env` (GATE_OUTCOME B).
2. Release 1.0.4 with Nextflow 24.04.4 in the new env `nf-24.04.4` (GATE_OUTCOME C); note the rMATS 4.1.2 hang seen once.
3. Stop the project.

G3-G9 then run with the chosen candidate only.

## NOT verified yet

Everything of G3-G8 for the chosen candidate: params-file typing, Salmon route with `aligner: "star"`, G4 switch-off and `reverse` mapping, G5 BAM input, G6 container HTTP checks, STAR/Salmon index versions, command lines, output tree fixture, gate_values.tsv. Whether the rMATS 4.1.2 hang recurs. Whether an unpatched 1.0.4 could run under 26.04.6 with a legacy-parser setting (not part of the plan's decision tree; not tried).
