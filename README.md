# bioinformatics-claude-skills

Claude Code custom skills (slash commands) for bioinformatics pipelines at the Koch Institute.

## What are skills?

Skills are Markdown files placed in `~/.claude/commands/` that define reusable, interactive workflows for Claude Code. They are invoked with a `/` prefix in the Claude Code CLI (e.g. `/nfcore-rnaseq-setup`).

## Installation

Clone this repository and copy the desired skill file(s) to your Claude commands directory:

```bash
git clone https://github.com/KochInstitute-Bioinformatics/bioinformatics-claude-skills.git
cp bioinformatics-claude-skills/nfcore-rnaseq-setup/nfcore-rnaseq-setup.md ~/.claude/commands/
```

The skill is immediately available — no restart required.

## Available skills

| Skill | Command | Description |
|-------|---------|-------------|
| nf-core/rnaseq setup | `/nfcore-rnaseq-setup` | Interactive setup wizard for nf-core/rnaseq bulk RNA-seq on SLURM + Singularity |
| nf-core/scrnaseq setup | `/nfcore-scrnaseq-setup` | Interactive setup wizard for nf-core/scrnaseq single-cell RNA-seq (cellranger / star / simpleaf / kallisto, CellBender, 10x v2–v4) on SLURM + Singularity |
| Bulk RNA-seq pipeline | `/bulk-rnaseq-pipeline` | Generates a full bulk RNA-seq downstream analysis pipeline (tximport → DESeq2/edgeR → GSEA) from nf-core/rnaseq star_salmon output as R Markdown + SLURM scripts |
| Seurat scRNA-seq pipeline | `/seurat-scrna-pipeline` | Generates a full single-cell RNA-seq downstream analysis pipeline (QC → Harmony integration → annotation → DEG/GSEA/LIANA) as R Markdown + SLURM scripts |

## Requirements

- [Claude Code](https://claude.ai/code) CLI installed and authenticated
- HPC cluster with SLURM scheduler
- Singularity/Apptainer available as a module
- Nextflow available via a conda environment
- Internet access from login node (for fetching Ensembl files and container images)

## Contributing

Add new skills as folders named after the slash command, each containing:
- `<skill-name>.md` — the skill definition
- `README.md` — documentation
