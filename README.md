# bioinformatics-claude-skills

Claude Code custom skills (slash commands) for RNA-seq pipelines at the Koch Institute — covering both **pipeline setup** (nf-core/rnaseq, nf-core/scrnaseq) and **downstream analysis** (bulk RNA-seq with DESeq2/edgeR, single-cell with Seurat), all wired for SLURM + Singularity on HPC.

## What are skills?

Skills are Markdown files placed in `~/.claude/commands/` that define reusable, interactive workflows for Claude Code. They are invoked with a `/` prefix in the Claude Code CLI (e.g. `/nfcore-rnaseq-setup`).

## Installation

Clone this repository and copy the desired skill file(s) to your Claude commands directory:

```bash
git clone https://github.com/KochInstitute-Bioinformatics/bioinformatics-claude-skills.git
cp bioinformatics-claude-skills/nfcore-rnaseq-setup/nfcore-rnaseq-setup.md ~/.claude/commands/
```

The skill is immediately available — no restart required.

Folder skills that bundle scripts (currently `handoff`) are copied as a whole folder into `~/.claude/skills/` instead:

```bash
cp -r bioinformatics-claude-skills/handoff ~/.claude/skills/
```

## Available skills

| Skill | Command | Description |
|-------|---------|-------------|
| nf-core/rnaseq setup | `/nfcore-rnaseq-setup` | Interactive setup wizard for nf-core/rnaseq bulk RNA-seq on SLURM + Singularity |
| nf-core/scrnaseq setup | `/nfcore-scrnaseq-setup` | Interactive setup wizard for nf-core/scrnaseq single-cell RNA-seq (cellranger / star / simpleaf / kallisto, CellBender, 10x v2–v4) on SLURM + Singularity |
| Bulk RNA-seq pipeline | `/bulk-rnaseq-pipeline` | Generates a full bulk RNA-seq downstream analysis pipeline (tximport → DESeq2/edgeR → GSEA) from nf-core/rnaseq star_salmon output as R Markdown + SLURM scripts |
| Seurat scRNA-seq pipeline | `/seurat-scrna-pipeline` | Generates a full single-cell RNA-seq downstream analysis pipeline (QC → Harmony integration → annotation → DEG/GSEA/LIANA) as R Markdown + SLURM scripts |
| Session handoff | `/handoff` | Writes a structured, validated end-of-session report under `.claude/reports/` and indexes it in the project `CLAUDE.md` so the next session can resume. Folder skill: install to `~/.claude/skills/` (see [handoff/README.md](handoff/README.md)) |

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

A skill that needs bundled scripts or reference files can instead be a folder skill: `SKILL.md` plus its supporting files and a `README.md`, installed to `~/.claude/skills/<skill-name>/`.
