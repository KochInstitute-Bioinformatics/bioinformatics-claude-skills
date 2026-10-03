# `/nfcore-rnavar-setup` — nf-core/rnavar Variant-Calling Setup Skill

An interactive Claude Code skill that walks you through setting up and submitting an [nf-core/rnavar](https://nf-co.re/rnavar) RNA-seq short-variant calling pipeline on an HPC cluster running SLURM and Singularity. rnavar implements the GATK best-practices RNA-seq workflow: 2-pass STAR, MarkDuplicates, SplitNCigarReads, base recalibration (BQSR), HaplotypeCaller, soft variant filtering, and optional SnpEff / VEP annotation. The skill follows the same wrap-an-nf-core-pipeline style as `/nfcore-rnaseq-setup`.

---

## Installation

```bash
cp nfcore-rnavar-setup.md ~/.claude/commands/
```

Then invoke it in Claude Code:

```
/nfcore-rnavar-setup
```

---

## Prerequisites

| Requirement | Notes |
|-------------|-------|
| SLURM scheduler | Script targets the bcc queue |
| Singularity ≥ 3.10 | Loaded via `module add singularity/3.10.4` |
| bcftools biocontainer | `depot.galaxyproject.org-singularity-bcftools-1.20--h8b25389_0.img` in `$NXF_SINGULARITY_CACHEDIR` (or `~/.singularity/cache`); the `prepare_known_sites` helper runs `bcftools`, `tabix` and `bgzip` from it via Singularity (there is no htslib/tabix module), downloading it once if absent |
| Nextflow ≥ 24.04 | Available in a conda environment; needed for the `resourceLimits` directive in the generated config (validated with 26.04.6) |
| Internet access | Required from the login node (GitHub API, Ensembl, known-sites downloads); compute nodes may not have it |

---

## What the skill does

The skill is a guided wizard that collects your settings one step at a time and generates a ready-to-submit SLURM script, plus helper scripts for any missing resources. The wizard asks only what it cannot detect.

| Step | Topic | Asked or auto-detected |
|------|-------|------------------------|
| 0 | Working directory | Auto (`pwd`); all outputs are written here |
| 1 | Email for SLURM notifications | Asked |
| 2 | Conda environment with Nextflow | Asked |
| 3 | nf-core/rnavar version | Auto (latest release via the GitHub API); if newer than 1.3.0 the parameter names are re-checked against that release's schema |
| 4 | Raw data and samplesheet | FASTQ/BAM/CRAM found automatically; paired-end vs single-end and sequencing date auto-detected; sample names sanitised; you review the sheet, choose naming and filename |
| 5 | Read length | Auto-detected from the reads (asked only for BAM/CRAM input) |
| 6 | Organism and genome files | Organism and genome base directory asked; choose an Ensembl release (version, existing FASTA/GTF and STAR index auto-detected) or a custom reference (you give the FASTA, GTF and index directory) |
| 7 | Known sites for base recalibration | Asked: supply known-sites VCFs or skip recalibration |
| 8 | Variant calling options | Asked: duplicates, two-pass, thresholds, gVCF, large chromosomes, save intermediates |
| 9 | Optional annotation | Asked: none, SnpEff, VEP or both; cache locations and identifiers asked |
| 10 | MultiQC title, output directory, `nextflow.config` | Titles asked; config written only if absent (literal `resourceLimits`, report overwrite enabled) |
| 11 | Params file and submission script | Generated |
| 12 | Helper scripts | Generated only for missing resources |
| 13 | Hand-off note | Prints where results will be, confirmed against the docs for the chosen version |

All pipeline parameters are written to a generated params YAML and the launch line carries only `-params-file`, so the exact settings of a run are recorded in one reviewable file and are not scattered over the command line. Parameters that equal the pipeline default are omitted.

### Key design points

- **Read length is always detected and written as `read_length` in the params file.** rnavar's default is 150, which is wrong for most other libraries, and it sets STAR's `sjdbOverhang` (`read_length − 1`).
- **A STAR index is built per read length**, in `{genome_base}/{organism}/{assembly}_ens{version}/index/star_rnavar_sjdb{N-1}/`. An index made for another read length (for example by another pipeline) is never reused. `star_index` is always set so rnavar does not rebuild it inside the workflow. FASTA/GTF are shared with the other nf-core skills' folder convention and reused, not re-downloaded.
- **Known sites are required, or base recalibration is explicitly skipped.** rnavar does not skip BQSR automatically when known sites are missing and errors late, after alignment. The wizard therefore asks you to supply dbSNP and known-indels VCFs (with `.tbi` indexes) or sets `skip_baserecalibration: true`. Download URLs are resolved at run time and shown to you for confirmation, never typed from memory. For GRCh38 the files are `Homo_sapiens_assembly38.dbsnp138.vcf.gz` and `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` from the GATK resource bundle; when the GATK page cannot be read (it returned HTTP 403), they are resolved from the listing of the public Google Cloud Storage bucket that hosts the bundle.
- **The pipeline head job is small but long:** `-n 2 --mem=8G -t 2-00:00:00`. It only coordinates the tasks, but it must outlive every one of them, so it asks for more than the usual 4 h on purpose.
- **Known-sites contigs must match the FASTA.** The wizard requires bgzipped `.vcf.gz` + `.tbi`. Contig names (Ensembl `1` vs GATK `chr1`) are checked against the FASTA in the wizard when the VCFs already exist; otherwise a contig guard in the `prepare_known_sites` helper script checks after download and either renames the contigs (`bcftools annotate --rename-chrs`) or fails fast with a non-zero exit. The pipeline is submitted with `--dependency=afterok`, so it cannot start after a failed guard.
- **Annotation caches must be pre-downloaded.** The skill does not silently set `download_cache`, because compute nodes may lack internet and the job then stalls. It offers a helper script to fetch the cache where internet is available.

---

## Output files

After running the skill you will have:

| File | Description |
|------|-------------|
| `{date}_{project}_samplesheet.csv` | Input samplesheet for nf-core/rnavar (FASTQ, BAM or CRAM form) |
| `nextflow.config` | SLURM + Singularity resource profiles (written only if none exists; an existing one is never overwritten) |
| `{date}_{project}_params.yaml` | All rnavar parameters for the run, passed with `-params-file` |
| `nf-core_rnavar_{version}.sh` | Pipeline SLURM submission script |
| `build_star_index_rnavar_{ref_tag}.sh` | STAR index build for this read length (if missing) |
| `prepare_known_sites_{ref_tag}.sh` | Known-sites download, `bgzip`, `tabix` indexing (if needed) and an always-run contig guard (rename with `bcftools annotate --rename-chrs`, or fail fast with non-zero exit), and a final check that every file named in the params file exists |
| `prepare_annotation_cache_{tool}.sh` | Annotation-cache pre-download, to run where internet is available (if no cache) |

---

## Typical workflow

```bash
# 1. Run any helper scripts that were generated (each is an sbatch script)
sbatch build_star_index_rnavar_GRCh38_ens115.sh

#    If the STAR helper downloads the FASTA, the known-sites helper (which reads it) must wait for it
sbatch --dependency=afterok:<star_jobid> prepare_known_sites_GRCh38_ens115.sh

# 2. Submit the pipeline, listing every generated helper
sbatch --dependency=afterok:<star_jobid>:<known_sites_jobid> nf-core_rnavar_1.3.0.sh
```

---

## Validation status

- The skill text is checked by `nfcore-rnavar-setup/tests/check_skill.sh`, which schema-validates every `--parameter` against nf-core/rnavar 1.3.0, checks that the launch line carries no `--flag` (including continuation lines) and that every key in the params template is a schema parameter, and asserts that required text is present and forbidden text is absent. It also checks the resource values of every `withName` selector and that each selector matches a process name of the real-data run (`tests/fixtures/realdata_gm12878_processes.txt`). `tests/test_star_tmp.sh` runs the STAR index block with stub `STAR` and `module` executables (temp directory placement and removal, exit codes, module check).
- **One end-to-end run on the nf-core rnavar test data was completed on the cluster** (Nextflow 26.04.6): paired-end reads, custom reference, known sites from local VCFs, no annotation. The known-sites helper logic (container wrappers for `bcftools`/`tabix`/`bgzip`, `tabix -f` re-run, contig-guard rename) was exercised in standalone batch tests on earlier helper text; the params-file check exercised then was the older `grep` check. The new params-file existence loop and the download/verify/move rule were exercised in standalone batch tests (job 11375528: container download and bad-URL failure, `--bind` de-duplication, truncated/valid/absent VCF download, params-file existence loop). The real-data test below then ran both generated helpers end to end (STAR index build; known-sites download from the GATK hg38 bucket, contig rename and params-file check). The Ensembl FASTA/GTF download path, the Mouse Genomes Project known-sites download, annotation, BAM/CRAM input and gVCF output remain untested.
- The `withName` selectors in the generated `nextflow.config` (`STAR_ALIGN`, `GATK4_SPLITNCIGARREADS`, `GATK4_BASERECALIBRATOR`, `GATK4_HAPLOTYPECALLER`) were verified to match real tasks in that run and again in the real-data test below. The `PICARD_MARKDUPLICATES` selector was added after the real-data test; it matches the process name in that run's trace (`tests/fixtures/realdata_gm12878_processes.txt`, which the checker uses), but its values have not yet been applied in a run.
- The 1.3.0 output layout in the hand-off note (`variant_calling/`, `preprocessing/`, `reports/`, `pipeline_info/`) was confirmed from that run; `annotation/` was not observed.

---

## Real-data test (one human dataset)

One end-to-end run on real data, with the skill at master 77bfb5d followed literally, on 2026-10-02. It is ONE dataset, one sample and one run; the numbers below hold for that dataset only.

- **Dataset:** SRR5665260 (BioProject PRJNA389940), GM12878 = GIAB HG001, whole-transcriptome RNA-seq, paired-end, 39,598,362 read pairs, NextSeq 500, reads pre-trimmed (151 bp mode). Library selection and strandedness are not stated by ENA.
- **What was run:** nf-core/rnavar 1.3.0, Nextflow 26.04.6, SLURM partition `bcc`. Ensembl 116 GRCh38 primary assembly and GTF (already on disk), a STAR index built by the skill's helper for sjdbOverhang 150, GATK hg38 bundle dbSNP138 and Mills/1000G indels (contigs renamed `chrN`→`N` by the helper's guard), default calling and filtering, no annotation.
- **Outcome:** 82/82 tasks completed, 0 failed, 0 retried. Pipeline wall time **3 h 10 min 42 s**; **22.9 CPU-hours** used by all tasks; 3 h 54 min end to end including the 43-minute STAR index build. STAR uniquely mapped 89.95%; Picard duplication 40.9%.

### Accuracy against the GIAB HG001 v4.2.1 truth (PASS calls)

Matching with `rtg vcfeval` (haplotype-aware), autosomes only. "genotype" requires the right genotype; "allele" (`--squash-ploidy`) ignores zygosity. Regions: GIAB confident regions ∩ Ensembl 116 exons (all biotypes, or protein-coding transcripts only) ∩ RNA depth ≥ 10 or ≥ 20. Indel means non-SNP (indels plus complex).

| Region | Mode | SNV P | SNV R | SNV F1 | Indel P | Indel R | Indel F1 |
|---|---|---|---|---|---|---|---|
| exons_all_dp10 | genotype | 0.885 | 0.917 | 0.901 | 0.674 | 0.815 | 0.738 |
| exons_all_dp10 | allele | 0.891 | 0.923 | 0.907 | 0.733 | 0.886 | 0.802 |
| exons_all_dp20 | genotype | 0.924 | 0.922 | 0.923 | 0.743 | 0.856 | 0.795 |
| exons_all_dp20 | allele | 0.927 | 0.925 | 0.926 | 0.800 | 0.922 | 0.857 |
| exons_pc_dp10 | genotype | 0.905 | 0.935 | 0.919 | 0.675 | 0.821 | 0.741 |
| exons_pc_dp10 | allele | 0.909 | 0.940 | 0.924 | 0.738 | 0.897 | 0.810 |
| exons_pc_dp20 | genotype | 0.937 | 0.940 | 0.939 | 0.738 | 0.866 | 0.797 |
| exons_pc_dp20 | allele | 0.940 | 0.943 | 0.942 | 0.796 | 0.935 | 0.860 |

PASS vs all records (PASS plus soft-filtered), exons_all_dp10, genotype mode:

| Set | SNV P | SNV R | SNV F1 | Indel P | Indel R | Indel F1 |
|---|---|---|---|---|---|---|
| PASS | 0.885 | 0.917 | 0.901 | 0.674 | 0.815 | 0.738 |
| all records | 0.766 | 0.957 | 0.851 | 0.534 | 0.823 | 0.648 |

**Evaluation limits.** The main region (exons_all_dp10) is 32.9 Mb, about 1.3% of the GIAB confident genome. No accuracy claim is made outside GIAB confident regions ∩ exons ∩ depth ≥ 10 (or ≥ 20), for other samples, libraries, tissues or read lengths, or for X/Y (the truth is autosomal).

**Call-set statistics:** 85,799 calls, 84.6% PASS; ts/tv 2.92 (PASS SNVs, autosomes); het/hom-alt 1.15; 72.9% of PASS calls in dbSNP138 (exact position, REF and ALT).

**What the errors were (measured; the editing-mask gain is inferred):**
- 78% of the false-positive SNVs were A>G/T>C (2,300 of 2,932), and 88.5% of those (2,036) sit on known REDIportal editing sites, against 1.37% of true-positive SNVs: they are RNA editing, real at the RNA level but absent from the DNA truth. Masking known editing sites would raise SNV precision to about 0.96 (an estimate, not re-run through vcfeval). The skill's hand-off note now recommends masking known editing sites (for example REDIportal).
- 76.5% of the missed truth variants were heterozygous. Of the 1,706 missed heterozygous SNVs, 841 (49%) showed strong allelic imbalance (no alt read or alt fraction < 0.2 at a median depth of 32×) and 777 (46%) were called but removed by the soft filters, 91% of the filtered ones by SnpCluster; the groups overlap partly. Depth was not the main cause: only 17 missed SNVs had depth < 10 at the site.

### Resources: requested vs observed (that run)

Requested values are those the tasks received in that run (the skill at 77bfb5d); the last column is this skill's request after the fix round. Observed values are from Nextflow's `execution_trace.txt` (peak_rss, realtime, %CPU); the helper memory is a lower bound from 2-minute `sstat` samples.

| Process | Requested then (cpu / mem / time) | Observed peak memory | Observed time | CPU used | Request now |
|---|---|---|---|---|---|
| STAR_ALIGN | 8 / 64 GB / 8 h | 43.0 GB (67%) | 49 min 09 s | 704% | 8 / 64 GB / 8 h (unchanged) |
| PICARD_MARKDUPLICATES | 6 / 36 GB / 8 h (rnavar label) | 27.8 GB (**77%**, the closest to its limit) | 20 min 47 s | 158% | 2 / 48 GB / 8 h (new selector) |
| GATK4_SPLITNCIGARREADS (25 scatters) | 2 / 16 GB / 8 h | 12.2 GB (**76%**) | at most 26 min 48 s | 335 to 746% | 4 / 24 GB / 4 h |
| GATK4_BASERECALIBRATOR | 2 / 16 GB / 8 h | 4.4 GB (27%) | 20 min 07 s | 110% | 2 / 8 GB / 4 h |
| GATK4_HAPLOTYPECALLER (25 scatters) | 2 / 16 GB / 8 h | 2.2 GB (14%) | at most 8 min 48 s | 212 to 249% | 2 / 8 GB / 4 h |
| RECALIBRATE:GATK4_APPLYBQSR | 2 / 12 GB / 4 h (rnavar label) | 3.5 GB (29%) | 26 min 56 s | 152% | unchanged (no selector) |
| Head job (Nextflow) | 32 cores, no memory or time | about 1 core | 3 h 11 min | about 1 core | `-n 2 --mem=8G -t 2-00:00:00` |
| STAR index helper (GRCh38) | 8 / 64G / 4 h | at least 43.1 GB (67%) | 42 min 27 s | n/a | unchanged |
| Known-sites helper | 2 / 8G / 4 h | about 45 MB | about 7 min | n/a | unchanged |

No task exceeded 80% of its memory or time. The head job's former 32-core request held about 101 core-hours against the 22.9 CPU-hours used by all tasks. The new 2-day head-job limit is a deliberate exception to the usual 4 h default: the job only coordinates, but it must outlive every task. Memory is not a consumable resource on the `bcc` partition, so memory requests and the SplitNCigarReads CPU oversubscription were advisory there; on a cluster that enforces memory they would be enforced.

**Disk use:** about 67 GB under the test folder (FASTQ, truth and REDIportal 7.6 GB; project with `work/` and outputs 58 GB; evaluation 1.7 GB), plus 33 GB written to the shared genome folder: STAR index 30 GB and `known_sites/` 2.9 GB.

### Defects found and fixed in this round

- D1, D2: the GATK resource-bundle page returned HTTP 403, and the skill named no files. Step 7 now falls back to the public bucket listing and names the two `.vcf.gz` objects to use (not the 11 GB plain `.vcf`).
- D8: the head job requested 32 cores with no memory or time; now `-n 2 --mem=8G -t 2-00:00:00`.
- D3: the sample-name rule now covers `_1.fastq.gz` (SRA/ENA names). D4: when no sequencing date is found, the identical date options are offered once. D5: every `module add` line in the helpers and the pipeline script is checked. D6: `wget -nv` keeps helper logs short (the known-sites log had 30,922 lines). D9: STAR's temporary directory goes next to the index (`--outTmpDir`) instead of the working directory, and is removed after the build. D7: the generic `cpus`/`memory`/`time` lines in the config are kept as a fallback for unlabelled processes, and the comment now says so (no task in the run used them; every rnavar process carries a label).
- Resource tiers updated as in the table above; hand-off note gains the RNA-editing and missed-variant guidance.

### Not verified

- The Ensembl FASTA/GTF download branch of the STAR helper and the bcftools container download branch (both files already existed).
- Annotation (SnpEff / VEP), BAM/CRAM input, gVCF, `remove_duplicates`, custom thresholds, multi-sample or multi-lane merging, single-end data.
- A cluster that enforces memory limits.
- Accuracy outside GIAB confident regions ∩ exons ∩ depth ≥ 10; X/Y; other samples, libraries or tissues.
- A second run for reproducibility (one run only).
- The fixes of this round (new tiers, head-job request, `--outTmpDir`, the bucket fallback) have not been run on real data; they were checked statically and with stubs.

---

## Known limitations

- Annotation caches (SnpEff / VEP) must be pre-downloaded; the skill does not rely on `--download_cache`
- Known sites are required, otherwise base recalibration is skipped, which is slightly less accurate
- A separate STAR index is built for each read length
- Parameter names were verified against rnavar 1.3.0 only; newer versions are re-checked at run time against the release schema
- The resource tiers are measured on one human dataset (one 39.6 M-pair library, 151 bp); much deeper libraries may need more memory or time. Compare the selectors with `pipeline_info/execution_trace.txt` after the first run, and after the first run on a new version
- RNA editing sites appear as A>G/T>C false positives against a DNA truth; mask known editing sites before treating them as genomic variants
- The head job requests 2 days; a longer run must be resumed with `-resume`
- Organisms other than mouse and human require manual FASTA/GTF paths and known-sites files
