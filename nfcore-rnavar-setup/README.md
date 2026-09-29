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
- **Known sites are required, or base recalibration is explicitly skipped.** rnavar does not skip BQSR automatically when known sites are missing and errors late, after alignment. The wizard therefore asks you to supply dbSNP and known-indels VCFs (with `.tbi` indexes) or sets `skip_baserecalibration: true`. Download URLs are resolved at run time and shown to you for confirmation, never typed from memory.
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

- The skill text is checked by `nfcore-rnavar-setup/tests/check_skill.sh`, which schema-validates every `--parameter` against nf-core/rnavar 1.3.0, checks that the launch line carries no `--flag` (including continuation lines) and that every key in the params template is a schema parameter, and asserts that required text is present and forbidden text is absent.
- **One end-to-end run on the nf-core rnavar test data was completed on the cluster** (Nextflow 26.04.6): paired-end reads, custom reference, known sites from local VCFs, no annotation. The known-sites helper logic (container wrappers for `bcftools`/`tabix`/`bgzip`, `tabix -f` re-run, contig-guard rename) was exercised in standalone batch tests on earlier helper text; the params-file check exercised then was the older `grep` check. The new params-file existence loop and the download/verify/move rule were exercised in standalone batch tests (job 11375528: container download and bad-URL failure, `--bind` de-duplication, truncated/valid/absent VCF download, params-file existence loop). The Ensembl FASTA/GTF download path, the known-sites download from real GATK/Mouse Genomes Project URLs, annotation, BAM/CRAM input, gVCF output and the generated helper end to end remain untested.
- The `withName` selectors in the generated `nextflow.config` (`STAR_ALIGN`, `GATK4_SPLITNCIGARREADS`, `GATK4_BASERECALIBRATOR`, `GATK4_HAPLOTYPECALLER`) were verified to match real tasks in that run.
- The 1.3.0 output layout in the hand-off note (`variant_calling/`, `preprocessing/`, `reports/`, `pipeline_info/`) was confirmed from that run; `annotation/` was not observed.

---

## Known limitations

- Annotation caches (SnpEff / VEP) must be pre-downloaded; the skill does not rely on `--download_cache`
- Known sites are required, otherwise base recalibration is skipped, which is slightly less accurate
- A separate STAR index is built for each read length
- Parameter names were verified against rnavar 1.3.0 only; newer versions are re-checked at run time against the release schema
- Process-name selectors were verified against real tasks only for the four listed in the validation status; compare them with `pipeline_info/execution_trace.txt` after the first run on a new version
- Organisms other than mouse and human require manual FASTA/GTF paths and known-sites files
