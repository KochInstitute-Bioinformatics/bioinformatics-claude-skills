# bioinformatics-claude-skills

[Claude Code](https://claude.ai/code) skills that help you **set up and run RNA-seq analyses on an HPC cluster (SLURM, PBS, SGE or LSF) or a single machine, with Singularity/Apptainer containers**. They cover:

- **Pipeline setup**: interactive wizards that configure and submit nf-core pipelines: rnaseq (bulk), scrnaseq (single-cell), rnavar (variant calling) and rnasplice (alternative splicing).
- **Downstream analysis**: skills that generate R Markdown reports and job scripts for bulk RNA-seq (DESeq2/edgeR, GSEA, optional DTU), alternative splicing (nf-core/rnasplice outputs) and single-cell (Seurat).
- **Allele-specific expression**: a wizard that runs its own alignment and counting jobs (no Nextflow) for F1 crosses and outbred/human samples.
- **Session handoff**: a structured end-of-session report so the next Claude Code session can resume.

> **Please read first**
> - These are **AI-agent instructions** (prompts for Claude Code). They are **unofficial and not endorsed by nf-core** or its maintainers.
> - They generate samplesheets, `nextflow.config` files and job scripts for your scheduler. **Review everything the agent writes before you submit it.** You are responsible for what you run on your cluster and your data.
> - They were developed and tested on one cluster (Koch Institute, MIT). Cluster-specific settings live in one config file (see [Configure your site](#configure-your-site)), but expect to adapt things to your environment. Issues and pull requests are welcome.

## Available skills

| Skill | Command | Description |
|-------|---------|-------------|
| nf-core/rnaseq setup | `/nfcore-rnaseq-setup` | Interactive setup wizard for nf-core/rnaseq bulk RNA-seq on your scheduler + Singularity |
| nf-core/scrnaseq setup | `/nfcore-scrnaseq-setup` | Interactive setup wizard for nf-core/scrnaseq single-cell RNA-seq (cellranger / star / simpleaf / kallisto, CellBender, 10x v2–v4) on your scheduler + Singularity |
| nf-core/rnavar setup | `/nfcore-rnavar-setup` | Interactive setup wizard for nf-core/rnavar RNA-seq variant calling (GATK best practices: 2-pass STAR, SplitNCigarReads, BQSR, HaplotypeCaller, optional SnpEff/VEP) on your scheduler + Singularity |
| nf-core/rnasplice setup | `/nfcore-rnasplice-setup` | Interactive setup wizard for nf-core/rnasplice differential alternative splicing (rMATS, SUPPA2, DEXSeq/edgeR exon usage, DEXSeq DTU; every module switch explicit) on your scheduler + Singularity |
| nf-core/rnasplice report | `/nfcore-rnasplice-report` | Interactive wizard that turns an nf-core/rnasplice output directory into R Markdown reports and job scripts: rMATS, SUPPA2, DEXSeq/edgeR exon usage and DTU tables and volcano plots, cross-tool comparison, GO enrichment and a summary page; every effect is treatment minus control. Needs no internet |
| Bulk RNA-seq pipeline | `/bulk-rnaseq-pipeline` | Generates a full bulk RNA-seq downstream analysis pipeline (tximport → DESeq2/edgeR → GSEA, optional DTU) from nf-core/rnaseq star_salmon output as R Markdown + job scripts |
| Seurat scRNA-seq pipeline | `/seurat-scrna-pipeline` | Generates a full single-cell RNA-seq downstream analysis pipeline (QC → Harmony integration → annotation → DEG/GSEA/LIANA) as R Markdown + job scripts |
| Allele-specific expression pipeline | `/ase-pipeline` | Interactive wizard for allele-specific expression from bulk RNA-seq (F1 cross with a third-allele masked reference, or outbred/human with STAR WASP): STAR alignment, GATK ASEReadCounter, per-sample beta-binomial tests, reference-bias diagnostic, R Markdown reports and a summary page on SLURM + Singularity. Analyses: per-sample imbalance, reciprocal F1 (strain and parent-of-origin effects) and differential ASE between conditions; in outbred mode, phASER haplotype counts reported next to the unphased results |
| Session handoff | `/handoff` | Writes a structured, validated end-of-session report under `.claude/reports/` and indexes it in the project `CLAUDE.md` so the next session can resume. Folder skill: install to `~/.claude/skills/` |

Each skill folder has its own README with the full step list and known limitations. The `nfcore-*-setup` skills fetch the **latest** nf-core release each time they run, so generated commands follow current pipeline parameters; still check them against the pipeline docs for the version you run.

## Quick start

```bash
# 1. Get the skills and install the ones you want
git clone https://github.com/KochInstitute-Bioinformatics/bioinformatics-claude-skills.git
cd bioinformatics-claude-skills
mkdir -p ~/.claude/commands
cp nfcore-rnaseq-setup/nfcore-rnaseq-setup.md ~/.claude/commands/

# 2. Configure your site (partition, modules, container paths) - see below
cp bioinfo-skills.conf.example ~/.claude/bioinfo-skills.conf
$EDITOR ~/.claude/bioinfo-skills.conf

# 3. In your project directory, start Claude Code and run the skill
cd /path/to/my/project && claude
> /nfcore-rnaseq-setup
```

Skills are plain Markdown files; no restart is needed after copying. Folder skills that bundle scripts (currently `handoff`) are copied as a whole folder into `~/.claude/skills/` instead:

```bash
cp -r bioinformatics-claude-skills/handoff ~/.claude/skills/
```

### Configure your site

Everything cluster-specific (scheduler and queue, resource caps, how to load conda/Singularity, where container images live, optional report URL prefix) is read from `~/.claude/bioinfo-skills.conf`. Any key you leave out is asked interactively, and the skill offers to save your answer. [`bioinfo-skills.conf.example`](bioinfo-skills.conf.example) documents every key; [`examples/koch-institute.conf`](examples/koch-institute.conf) is a filled-in example.

## Schedulers

Set `SCHEDULER` in `~/.claude/bioinfo-skills.conf` (the skill asks if it is missing). The skills write their scripts in SLURM form and adapt them to the scheduler with a table that is embedded in every skill (`docs/scheduler-table.md`); `local` writes plain `bash` scripts without job headers and tells you the order to run them in. `MAX_TIME`, `MAX_MEM_GB` and `MAX_CPUS` cap what the skills request (for example 12 h on a public partition).

| Scheduler | Status |
|-----------|--------|
| `slurm` | tested |
| `local` | tested |
| `pbs`, `sge`, `lsf` | translated from the schedulers' documentation, not yet run on a real cluster: please report problems |

`ase-pipeline` builds its own job chain and supports SLURM only for now. Examples of configs: `examples/koch-institute.conf` (SLURM) and `examples/mit-orcd-engaging.conf` (MIT ORCD Engaging, SLURM).

## Example

Prompt (in a directory containing paired-end `*.fastq.gz` files):

```
/nfcore-rnaseq-setup
```

The skill asks one question at a time (illustrative, abbreviated):

```
Working directory: /path/to/Project1
Which email should the scheduler notify?            > me@example.org
Which conda environment has Nextflow?       > nf-env
Latest nf-core/rnaseq release: 3.x.y
Found 12 paired-end samples, sample sheet:  (table shown for review)
Trimming?  1) already trimmed  2) fastp (recommended)  3) trimgalore   > 2
Organism?  1) mouse  2) human  3) other                                > 1
Aligner?   1) star_salmon (recommended) ...                            > 1
```

Result: `samplesheet.csv`, `nextflow.config`, and a launch script. **Open and read them**, then submit:

```
Script written: nf-core_rnaseq_3.x.y.sh
To submit:  sbatch nf-core_rnaseq_3.x.y.sh
```

## Requirements

- [Claude Code](https://claude.ai/code) CLI installed and authenticated
- A scheduler (SLURM, PBS, SGE or LSF) or a single machine (`SCHEDULER=local`), with Singularity/Apptainer available
- Nextflow in a conda environment (for the `nfcore-*` setup skills; `/ase-pipeline` does not use Nextflow)
- Internet access from the login node (Ensembl files, container images); `/ase-pipeline` downloads any missing biocontainer image in its prep job on a compute node; `/nfcore-rnasplice-setup` needs internet on the compute nodes: its job downloads the pinned nf-core/rnasplice revision and the containers
- R is not needed on the host: the analysis skills run R inside a Singularity container, never on a login node. The images used by default and their Dockerfiles are in [`containers/`](containers)

## Contributing and feedback

Bug reports, missing options, and requests for other nf-core pipelines are welcome: please [open an issue](https://github.com/KochInstitute-Bioinformatics/bioinformatics-claude-skills/issues). See [CONTRIBUTING.md](CONTRIBUTING.md).

## License

[MIT](LICENSE)
