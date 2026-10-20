# Contributing

Thanks for helping improve these skills. Issues and pull requests are both welcome.

## Reporting problems

Open an issue with:
- the skill name and the pipeline version it targeted (e.g. nf-core/rnaseq 3.x.y),
- what the generated file or command was, and what was wrong or missing,
- your scheduler/container setup if relevant (SLURM version, Singularity vs Apptainer).

**Never paste credentials, internal hostnames, or unpublished sample/project names.**

## Changing or adding a skill

Each skill is a folder named after its slash command:

- `<skill-name>.md`: the skill definition (the prompt Claude Code follows)
- `README.md`: documentation (what it does, prerequisites, steps, known limitations)

Guidelines:
- **No site-specific values in skills.** Partition names, module lines, paths and image locations come from `~/.claude/bioinfo-skills.conf` (add a key to `bioinfo-skills.conf.example` if you need a new one).
- Keep generated output reviewable: the skill should show files it writes and never submit jobs itself.
- Prefer numbered-option questions over open-ended ones when the answer set is finite.
- Say which pipeline versions you tested against.

## Running a skill's own tests

The `nfcore-rnavar-setup`, `nfcore-rnasplice-setup`, `nfcore-rnasplice-report` and `ase-pipeline` folders ship a `tests/` directory. Most of it is static (bash, awk and sed against the skill text, with stub `sbatch`, `module`, `singularity` and `nextflow`) and can be run from the skill folder:

```bash
cd nfcore-rnasplice-setup && bash tests/run_all_tests.sh        # about an hour: submit it with sbatch, not on a login node
cd nfcore-rnavar-setup && bash tests/check_skill.sh nfcore-rnavar-setup.md
```

A few tests run R or tools in a container on a compute node (`nfcore-rnasplice-report/tests/test_r.sh`, `ase-pipeline/tests/r/run_stats_tests.sh`, the `ase-pipeline/tests/synthetic/` scripts). They read your site settings from `~/.claude/bioinfo-skills.conf` (or the file named in `BIOINFO_SKILLS_CONF`) and need a scratch directory that the compute nodes can see:

| Variable | Used by | Meaning |
|----------|---------|---------|
| `BIOINFO_SKILLS_CONF` | all of the above | alternative path of the site config |
| `RSR_SHARED_TMP` | `nfcore-rnasplice-report` R tests and mutation proof | scratch directory visible from compute nodes (required) |
| `RSR_PARTITION`, `ASE_PARTITION` | the same | override `QUEUE` for the test jobs |
| `ASE_SCRATCH_ROOT` | `ase-pipeline/tests/chain/dry_run_chain.sh` | the only tree under which that test deletes and re-creates its scratch directory |

When you change text that a test asserts (for example a generated script line), change the test in the same commit. The tests compare against the placeholders (`{QUEUE}`, `{SINGULARITY_LOAD}`, `{CONDA_SETUP}`, `{PCLI_FLAG}`), never against site values.
