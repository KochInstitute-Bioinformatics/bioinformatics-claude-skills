# Scheduler support for the bioinformatics skills — design

Status: ready for approval. Items marked **(my call)** are decisions I made, not ones you made. Decided by you: ask when `SCHEDULER` is missing; rename `SLURM_PARTITION` to `QUEUE`; round 1 excludes ASE; Engaging is a SLURM configuration and you have access to it.

## Intent

The skills write job scripts and submit them. Today every one of them assumes SLURM. The goal is that a user on another scheduler, or on a workstation with no scheduler, can use the skills without editing them. "As easy and usable by everyone as possible" is the stated success measure; there is no specific cluster in mind. MIT ORCD Engaging was named as a target; it runs SLURM, so it is a configuration, not a scheduler.

## Goals

1. A new site-config key `SCHEDULER` with values `slurm`, `pbs`, `sge`, `lsf`, `local`. When it is missing the skill asks (decided), so nobody silently gets SLURM scripts.
2. Every generated job script uses the right header lines, submit command, dependency syntax and log options for that scheduler. With `SCHEDULER=slurm` the generated output stays byte-identical to today (regression guard: all existing tests keep passing unchanged except for renamed keys).
3. `local`: no job headers; the skills write plain `bash` scripts, state the order to run them in, and for the Nextflow head job advise `nohup`/`tmux`. Nextflow gets `executor = 'local'`.
4. Resource caps: new keys `MAX_TIME`, `MAX_MEM_GB`, `MAX_CPUS` (empty = no cap). Generated requests are clamped to them. When a request that must be long (the Nextflow head job: 48 h in rnasplice, 2 days in rnavar and scrnaseq) exceeds `MAX_TIME`, the skill requests `MAX_TIME`, and tells the user to resubmit with `-resume` (rnavar already documents this) or to use a longer queue. This matters on Engaging, whose public partitions allow 12 h (`mit_normal`) and 48 h (`mit_preemptable`).
5. Scheduler-neutral key names: `QUEUE` replaces `SLURM_PARTITION`, `QUEUE_ON_CLI` replaces `SLURM_PARTITION_ON_CLI`. The old names are still read as aliases (decided). Placeholders in the skills and tests become `{QUEUE}`.
6. Ready-made example configs: `examples/mit-orcd-engaging.conf` next to `koch-institute.conf`.
7. The README states per scheduler what has been tested.
8. Verification on Engaging (you have access): a short real run of the generated scripts there, and confirmation of the Apptainer module name for the example config. These runs happen after the implementation; the example config ships with the verified module line or with the blank line and a note if verification is not done.

## Non-goals (round 1)

- `ase-pipeline` (round 2): it builds its own job chain with arrays, `--dependency`, `squeue` and `scancel`, which each scheduler spells differently. It stays SLURM-only and the README says so.
- Container engines other than Singularity/Apptainer (Docker, Podman, conda-only), cloud batch services, Kubernetes, GPU requests.
- Choosing nf-core institutional profiles (`-profile engaging`). The skills keep generating their own `nextflow.config`.

## Anti-goals

- Claiming support that was not exercised. PBS, SGE and LSF can only be tested here through stubs and header translation, not on a real cluster. The README marks them "translated from scheduler documentation, not yet run on a real cluster" until someone confirms; `slurm` is "tested" and `local` is "tested" (it can run here).
- Silently wrong scripts: a SLURM option with no mapping must stop the skill with a message, not be dropped.

## Constraints

- A skill is one Markdown file copied to `~/.claude/commands/`; it cannot read other files of the repo at run time. Whatever table it needs is embedded in it.
- Existing tests assert exact script text and have mutation suites; they must keep passing, and new behavior gets its own tests and mutation rows.
- No heavy work on login nodes (user rule), so test suites run through `sbatch`.

## How it works

**One source of truth.** `docs/scheduler-table.md` holds, per scheduler: header line forms for job name, queue, nodes/cpus, memory, walltime, log file, mail, array; the submit command; how to read the job id; dependency flag; and notes on unit differences (SLURM `D-HH:MM:SS` vs hours; SGE memory is per slot; LSF memory in MB and walltime `HH:MM`). A generator script embeds the same table block, verbatim, into every skill's "Site configuration" section, and a test fails if any skill's block differs from the source.

**Skills whose scripts the model writes** (rnaseq, scrnaseq, rnavar, rnasplice-setup, bulk, Seurat): the templates stay in SLURM form and the skill says: for any other `SCHEDULER`, rewrite each `#SBATCH` line and `sbatch` command through the table, convert units as the table says, and stop on an option with no row. **(my call: lowest churn; the SLURM text stays the tested reference.)** A coverage test extracts every `#SBATCH` option and `sbatch` flag used in the templates of these skills and checks that each has a row in the table.

**Skills that render scripts with their own code** (rnasplice-report: `sed` templates and `run_all.tpl`): the render block gets a deterministic translation step for the header and submit lines, tested with stub schedulers for all five values.

**Nextflow config.** `executor`, `queue` and the `queueSize`/`submitRateLimit` settings come from `SCHEDULER`/`QUEUE`: `slurm`, `pbspro` or `pbs`, `sge`, `lsf`, `local`. **(my call: PBS maps to `pbspro` by default; a `PBS_FLAVOR=torque` key switches to `pbs`.)**

**Docs.** README gets a scheduler section and a support-status table. `examples/mit-orcd-engaging.conf`: `SCHEDULER=slurm`, `QUEUE=mit_normal`, `MAX_TIME=12:00:00`. The Apptainer module name on Engaging could not be confirmed from public documentation; the file leaves `SINGULARITY_LOAD` with a comment to check `module avail apptainer`.

## Left to the builder

Exact wording of table rows and of the per-skill instruction; the generator script's language (bash/awk, no Python on the login nodes); test design for stubs; naming of test files.

## Plan of rounds

Round 1 (this design): table, generator, the 8 non-ASE skills, rnasplice-report's deterministic translation, tests, mutation rows, README, Engaging example. Round 2: `ase-pipeline`.

## Decided questions

1. Missing `SCHEDULER`: the skill asks.
2. `SLURM_PARTITION` is renamed `QUEUE` (`SLURM_PARTITION_ON_CLI` to `QUEUE_ON_CLI`), old names kept as aliases.
3. Engaging can be tested for real by the user; PBS, SGE and LSF have no test cluster and ship labeled untested.
