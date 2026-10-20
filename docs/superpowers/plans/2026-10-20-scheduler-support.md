# Scheduler Support Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let the seven job-writing skills run on SLURM, PBS, SGE, LSF or no scheduler (`local`), chosen by `SCHEDULER` in `~/.claude/bioinfo-skills.conf`, with resource caps and a MIT ORCD Engaging example.

**Architecture:** One source-of-truth Markdown table (`docs/scheduler-table.md`) is embedded verbatim by a generator into each skill between marker comments. Skills keep their SLURM templates as the reference and tell the model to rewrite them through the table (stopping on an unmapped option); `nfcore-rnasplice-report`, whose scripts are rendered by Bash blocks, gets a tested deterministic translator. `ase-pipeline` is untouched (round 2).

**Tech Stack:** Bash and awk only (the login nodes have no usable Python; tests under `env -i` with stubs), Markdown, the existing `check_skill.sh` / mutation-proof suites.

**Spec:** `docs/superpowers/specs/2026-10-20-scheduler-support-design.md`

## Global Constraints

- Seven skills are in scope: `nfcore-rnaseq-setup`, `nfcore-scrnaseq-setup`, `nfcore-rnavar-setup`, `nfcore-rnasplice-setup`, `nfcore-rnasplice-report`, `bulk-rnaseq-pipeline`, `seurat-scrna-pipeline`. `ase-pipeline` and `handoff` are not touched.
- `SCHEDULER` values: `slurm`, `pbs`, `sge`, `lsf`, `local`. A missing key: the skill asks (numbered list). Never assume `slurm`.
- `QUEUE` replaces `SLURM_PARTITION`; `QUEUE_ON_CLI` replaces `SLURM_PARTITION_ON_CLI`; the old names are read as aliases when the new ones are absent. Placeholders `{SLURM_PARTITION}` become `{QUEUE}`; `{PCLI_FLAG}` keeps its name.
- New keys: `MAX_TIME` (`[D-]HH:MM:SS`), `MAX_MEM_GB`, `MAX_CPUS`, `PBS_FLAVOR` (`pbspro` default, or `torque`). Empty = no cap.
- With `SCHEDULER=slurm` the generated scripts equal today's. Existing `check_skill.sh` assertions on SLURM text keep passing; only renamed keys change.
- A SLURM option with no row in the table stops the skill with a message; it is never dropped.
- README status labels: `slurm` tested, `local` tested, `pbs`/`sge`/`lsf` "translated from scheduler documentation, not yet run on a real cluster".
- No heavy work on login nodes: long suites run through `sbatch` (default requests at most 64G and 4h). No push without the user's explicit approval; commit messages end with `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`.
- Results of Engaging checks go to the user to run; do not invent module names.

## Review Focus

1. `SCHEDULER` written as `SLURM`, ` slurm ` or `torque`-like typos: trim and lowercase it; anything not in the five values is asked again, never guessed.
2. An old config that has only `SLURM_PARTITION` / `SLURM_PARTITION_ON_CLI`: behaves exactly as `QUEUE` / `QUEUE_ON_CLI`.
3. `MAX_TIME` lower than a fixed request (2-day or 48 h head job): the request becomes `MAX_TIME`, `-resume` guidance is printed; times compared after converting `D-HH:MM:SS` to seconds.
4. SGE memory is per slot: `--mem=64G` with `-n 8` becomes `h_vmem=8G`.
5. `local` never prints or writes `sbatch`, `qsub`, `bsub` or `#SBATCH`/`#PBS`/`#$`/`#BSUB` lines, and runs stages strictly in order, stopping at the first failure.
6. Empty `QUEUE`: no queue line or `-q`/`-p` option on any scheduler.

## Lean rulings (user chose the lean version on 2026-10-20)

- Ruling: no deterministic translator and no `test_sched_translate.sh` (old Task 6) — the model adapts all seven skills, `nfcore-rnasplice-report` included, from the embedded table; its Step 10 gets a `**Scheduler adaptation.**` paragraph (Task 4) telling the model to rewrite the header lines of `run_job.sbatch.tpl` and the `JOB` lines of `run_all.tpl` through the table before the render block runs — cost if wrong: non-slurm report scripts vary between runs and are untested.
- Ruling: `tests/check_scheduler.sh` has parts (a) and (c) only; the construct-coverage part (b) and the extra `sched-*` mutation rows of Tasks 3 and 4 are not written. The simple `need` assertions of Task 3 and Task 4 stay.
- Ruling: Task 4's list adds `nfcore-rnasplice-report` (paragraph in Step 10).
- Work happens inline on branch `scheduler-support`; the final review is a self-review unless a fresh reviewer is requested.

---

### Task 1: Rename `SLURM_PARTITION` to `QUEUE`

**Files:**
- Modify: all seven in-scope skill `.md` files, `ase-pipeline/ase-pipeline.md` and its tests, every `*/README.md`, `README.md`, `CONTRIBUTING.md`, `bioinfo-skills.conf.example`, `examples/koch-institute.conf`, every `tests/` file that mentions the old names (`grep -rl SLURM_PARTITION`).

**Interfaces:**
- Produces: placeholders `{QUEUE}`, keys `QUEUE`, `QUEUE_ON_CLI`, template marker `__QUEUE__` (was `__PARTITION__`) and shell variable `QUEUE` (was `PARTITION`) in `nfcore-rnasplice-report`; the alias sentence in the shared rules block (Task 3).

- [ ] **Step 1:** In a scratch copy count `grep -rIc 'SLURM_PARTITION'` per file; this list is the checklist.
- [ ] **Step 2:** Replace `SLURM_PARTITION_ON_CLI` by `QUEUE_ON_CLI`, then `SLURM_PARTITION` by `QUEUE`, then in `nfcore-rnasplice-report` only `__PARTITION__` by `__QUEUE__` and its `PARTITION=` variable by `QUEUE=` (skill, `tests/test_render_slurm.sh`, `tests/check_skill.sh`, `tests/mutations.tsv`). Keep `ase-pipeline` consistent (it uses `{SLURM_PARTITION}` and `{PCLI_FLAG}`): rename there too, and `assemble_outbred_from_skill.sh` header comment.
- [ ] **Step 3:** In `bioinfo-skills.conf.example` document the aliases: a comment `# SLURM_PARTITION / SLURM_PARTITION_ON_CLI are still read when QUEUE / QUEUE_ON_CLI are absent`.
- [ ] **Step 4:** Run the four static checkers (`check_skill.sh`) and `test_render_slurm.sh`, `dry_run_submit.sh`, `dry_run_chain.sh` (with `ASE_SCRATCH_ROOT` set). Expected: PASS lines as before.
- [ ] **Step 5:** `grep -rIn 'SLURM_PARTITION' . --exclude-dir=.git` must print only the alias comment lines. Commit `Rename SLURM_PARTITION to QUEUE`.

### Task 2: Scheduler table, generator and root tests

**Files:**
- Create: `docs/scheduler-table.md`, `tools/embed-scheduler.sh`, `tests/check_scheduler.sh`, `tests/test_embed_scheduler.sh`

**Interfaces:**
- Produces: in `docs/scheduler-table.md`, two regions delimited by the lines `<!-- SCHEDULER:BEGIN table -->` / `<!-- SCHEDULER:END table -->` and `<!-- SCHEDULER:BEGIN rules -->` / `<!-- SCHEDULER:END rules -->`. The table's header row is `| SLURM construct | pbs | sge | lsf | local |`; the first column holds the SLURM construct in backticks, exactly as it appears in skills (`` `#SBATCH -p` ``).
- Produces: `tools/embed-scheduler.sh [--check] <skill.md>...` — for each file, replaces the content between the same marker pairs with the source content; if a pair is absent, inserts both blocks immediately before the first line starting `**Safety notice` (the skill's "Site configuration" section). With `--check` it changes nothing and exits 1 if any file differs. Exit 2 on a file with no `**Safety notice` line and no markers.
- Produces: `tests/check_scheduler.sh` — exits 0 and prints `SCHEDULER CHECK PASS` when (a) every in-scope skill's two blocks equal the source (calls the generator with `--check`), (c) `bioinfo-skills.conf.example` documents `SCHEDULER`, `QUEUE`, `QUEUE_ON_CLI`, `MAX_TIME`, `MAX_MEM_GB`, `MAX_CPUS`, `PBS_FLAVOR`.
- Table rows (the SLURM constructs the seven skills use; no automated coverage test): `#SBATCH -J`, `-p`, `-N`, `-n`, `-c`, `--mem`, `-t`, `-o`, `--mail-type`, `--mail-user`; `sbatch`, `sbatch --parsable`, `sbatch --dependency`, `squeue -u`, `scancel`, `sacct`, `srun`/`salloc`; Nextflow `executor`.

- [ ] **Step 1: Write `tests/test_embed_scheduler.sh`** (under `env -i`, in a `mktemp -d` under `/tmp/tmp.XXXXXXXXXX`): a fake skill with a `**Safety notice` line and no markers → generator inserts both blocks before it (names `inserted`); running again changes nothing (`idempotent`); editing text inside a block then `--check` exits 1 (`drift detected`) and a plain run restores it (`repaired`); a file without the anchor exits 2 (`no anchor`). Prints `EMBED TEST PASS (4)`.
- [ ] **Step 2:** Run it. Expected: FAIL (generator missing).
- [ ] **Step 3: Implement `tools/embed-scheduler.sh`** in bash and awk; the source is read from `$(dirname "$0")/../docs/scheduler-table.md`.
- [ ] **Step 4:** Write `docs/scheduler-table.md`. Table rows (SLURM construct → pbs / sge / lsf / local), values exactly:
  - `#SBATCH -J NAME` → `#PBS -N NAME` / `#$ -N NAME` / `#BSUB -J NAME` / omit
  - `#SBATCH -p Q` → `#PBS -q Q` / `#$ -q Q` / `#BSUB -q Q` / omit
  - `#SBATCH -N 1` → merged into the select line / omit / omit / omit
  - `#SBATCH -n N`, `#SBATCH -c N` → `#PBS -l select=1:ncpus=N` (memory joined as `:mem=Mgb`) / `#$ -pe smp N` / `#BSUB -n N` plus `#BSUB -R "span[hosts=1]"` / omit
  - `#SBATCH --mem=MG` → `mem=Mgb` in the select line / `#$ -l h_vmem=(M/N)G` (per slot) / `#BSUB -M (M*1024)` and `#BSUB -R "rusage[mem=(M*1024)]"` (MB; site may use another unit — say so) / omit
  - `#SBATCH -t D-HH:MM:SS` → `#PBS -l walltime=HH:MM:SS` (hours = 24*D+HH) / `#$ -l h_rt=HH:MM:SS` / `#BSUB -W HH:MM` / omit
  - `#SBATCH -o FILE` → `#PBS -o FILE` + `#PBS -j oe` / `#$ -o FILE` + `#$ -j y` / `#BSUB -o FILE` / run as `bash script.sh > FILE 2>&1`
  - `#SBATCH --mail-type=END,FAIL` + `--mail-user=A` → `#PBS -m ae` + `#PBS -M A` / `#$ -m ea` + `#$ -M A` / `#BSUB -N` + `#BSUB -u A` / omit (`END` only: `e`/`e`/`-N`; `FAIL` only: `a`/`a`/omit)
  - `sbatch SCRIPT` → `qsub SCRIPT` / `qsub SCRIPT` / `bsub < SCRIPT` / `bash SCRIPT`
  - `sbatch --parsable SCRIPT` (prints the job id) → `qsub SCRIPT` (prints the id) / `qsub -terse SCRIPT` / `bsub < SCRIPT` then `sed -n 's/^Job <\([0-9]*\)>.*/\1/p'` / not applicable
  - `sbatch --dependency=afterok:ID SCRIPT` → `qsub -W depend=afterok:ID SCRIPT` / `qsub -hold_jid ID SCRIPT` / `bsub -w "done(ID)" < SCRIPT` / run after the previous script succeeded (`&&`)
  - `squeue -u $USER` → `qstat -u $USER` / `qstat -u $USER` / `bjobs` / `jobs`
  - `scancel ID` → `qdel ID` / `qdel ID` / `bkill ID` / `kill`
  - `sacct` → `qstat -x` / `qacct -j ID` / `bhist -l ID` / not applicable
  - `srun` / `salloc` (prose, debugging) → `qsub -I` / `qrsh` / `bsub -Is bash` / run directly
  - Nextflow `executor` → `pbspro` (`pbs` when `PBS_FLAVOR=torque`) / `sge` / `lsf` / `local`; the profile name in `-profile NAME,singularity` and in `profiles { NAME {` equals the scheduler value (`slurm`, `pbs`, `sge`, `lsf`, `local`)
  Rules block text (one paragraph each, imperative): keys and aliases (`SCHEDULER`, `QUEUE`, `QUEUE_ON_CLI`, `MAX_TIME`, `MAX_MEM_GB`, `MAX_CPUS`, `PBS_FLAVOR`; `SLURM_PARTITION[_ON_CLI]` read when the new name is absent); trim and lowercase `SCHEDULER`, ask (numbered) when missing or not one of the five; rewrite every SLURM construct through the table, stop with `No row for <construct> in the scheduler table` on an unmapped one; clamp requests to the caps (compare times in seconds), and for a request that must be longer than `MAX_TIME` write `MAX_TIME`, tell the user to resubmit with `-resume` or choose a longer queue; `local` mode writes no header lines and no scheduler commands, lists the scripts in the order to run them (`bash a.sh && bash b.sh`), and for the Nextflow head job recommends `nohup bash script.sh > script.log 2>&1 &` or `tmux`; `{PCLI_FLAG}` is the queue option of the submit command (`-p Q ` for `sbatch`, `-q Q ` for `qsub` and `bsub`), empty when `QUEUE_ON_CLI` is not `yes` or `QUEUE` is empty.
- [ ] **Step 5:** Run `tests/test_embed_scheduler.sh`. Expected: `EMBED TEST PASS (4)`.
- [ ] **Step 6: Write `tests/check_scheduler.sh`** per the interface; run it. Expected: FAIL at (a) (skills not yet embedded) and PASS on (b), (c) once the conf example is updated in Task 7 (until then expect FAIL at (c) too).
- [ ] **Step 7:** Commit `Add scheduler table, generator and checks`.

### Task 3: Embed the table and rules in the seven skills; pin them in the skill tests

**Files:**
- Modify: the seven skill `.md` files (via the generator); `nfcore-rnavar-setup/tests/check_skill.sh`, `nfcore-rnasplice-setup/tests/check_skill.sh`, `nfcore-rnasplice-report/tests/check_skill.sh`; `*/tests/mutations.tsv` of rnasplice-setup and rnasplice-report.

**Interfaces:**
- Consumes: `tools/embed-scheduler.sh` (Task 2).
- Produces: each skill contains the marker pair `<!-- SCHEDULER:BEGIN rules -->` and `<!-- SCHEDULER:BEGIN table -->` in its Site configuration section.

- [ ] **Step 1: Failing tests.** In the three `check_skill.sh` files add `need '<!-- SCHEDULER:BEGIN table -->'`, `need '<!-- SCHEDULER:BEGIN rules -->'`, `need 'No row for <construct> in the scheduler table'`, `need 'resubmit with \`-resume\`'`. Run each checker. Expected: FAIL missing required text.
- [ ] **Step 2:** `tools/embed-scheduler.sh` on the seven skills; `tests/check_scheduler.sh` part (a) now passes.
- [ ] **Step 3:** Re-run the three checkers. Expected: PASS.
- [ ] **Step 4:** Commit `Embed the scheduler table and rules in the skills`.

### Task 4: Nextflow setup skills — executor, profile, head-job cap, local mode

**Files:**
- Modify: `nfcore-rnaseq-setup/nfcore-rnaseq-setup.md`, `nfcore-scrnaseq-setup/nfcore-scrnaseq-setup.md`, `nfcore-rnavar-setup/nfcore-rnavar-setup.md`, `nfcore-rnasplice-setup/nfcore-rnasplice-setup.md`, `nfcore-rnasplice-report/nfcore-rnasplice-report.md` (Step 10, before the render block); checkers of rnavar and rnasplice-setup.

**Interfaces:**
- Consumes: the rules block (Task 3).
- Produces: in each skill, directly after the launch-script code block, a paragraph starting `**Scheduler adaptation.**` stating: for `SCHEDULER` other than `slurm`, the `nextflow.config` profile block, its `executor` line, the `-profile` name and the script's header lines follow the table; for `local` the profile uses `executor = 'local'`, the script has no header and the final message says `bash <script>` (or the `nohup` form) instead of the submit command; the head job's walltime (`-t 2-00:00:00` rnavar and scrnaseq, `-t 48:00:00` rnasplice-setup, none in rnaseq) is clamped to `MAX_TIME` as the rules say.

- [ ] **Step 1: Failing tests.** In the rnavar and rnasplice-setup checkers add `need '**Scheduler adaptation.**'` and `need "executor = 'local'"`. Run. Expected: FAIL.
- [ ] **Step 2:** Add the paragraph to the four skills (the same wording in all four, with the skill's own head-job time).
- [ ] **Step 3:** Run the checkers and `dry_run_submit.sh` of rnasplice-setup. Expected: PASS (the slurm templates are unchanged).
- [ ] **Step 4:** Commit `Nextflow setup skills: scheduler adaptation`.

### Task 5: Downstream skills (bulk, Seurat) — local run order

**Files:**
- Modify: `bulk-rnaseq-pipeline/bulk-rnaseq-pipeline.md`, `seurat-scrna-pipeline/seurat-scrna-pipeline.md`

**Interfaces:**
- Consumes: the rules block. Produces: in each skill, after the "Final submission instructions" list, a paragraph `**Scheduler adaptation.**`: other schedulers rewrite each stage script through the table and print the submit commands with the table's submit command; `local` writes scripts without headers, runs nothing itself, and prints `bash run_01...sh && bash run_02...sh ...` in stage order, with a note that the `singularity exec` lines stay.

- [ ] **Step 1:** These skills have no test suites. Add the paragraph, then run `tests/check_scheduler.sh` part (a) (still PASS: the paragraph is outside the markers) and `grep -c 'Scheduler adaptation' ` on both files (expected 1 each).
- [ ] **Step 2:** Commit `Downstream skills: scheduler adaptation`.

### Task 6: (dropped — see Lean rulings)

### Task 7: Docs and example configs

**Files:**
- Modify: `bioinfo-skills.conf.example`, `README.md`, `CONTRIBUTING.md`, `examples/koch-institute.conf`
- Create: `examples/mit-orcd-engaging.conf`, `docs/engaging-verification.md`

**Interfaces:**
- Produces: documented keys `SCHEDULER`, `QUEUE`, `QUEUE_ON_CLI`, `MAX_TIME`, `MAX_MEM_GB`, `MAX_CPUS`, `PBS_FLAVOR` in the conf example; `examples/koch-institute.conf` gains `SCHEDULER=slurm`, `QUEUE=bcc`, `QUEUE_ON_CLI=yes`.
- `examples/mit-orcd-engaging.conf`: `SCHEDULER=slurm`, `QUEUE=mit_normal`, `MAX_TIME=12:00:00`, `MAX_CPUS=96`, and `SINGULARITY_LOAD=` empty with the comment `# check with: module avail apptainer (name not confirmed)`.

- [ ] **Step 1:** Run `tests/check_scheduler.sh`. Expected: FAIL at (c).
- [ ] **Step 2:** Update the conf example and the two example configs.
- [ ] **Step 3:** README: replace "HPC with SLURM" requirement lines with a Schedulers section and the status table (labels from Global Constraints; `ase-pipeline` listed as SLURM-only for now); rewrite the "Please read first" bullet to say the skills generate job scripts for the configured scheduler; CONTRIBUTING: how to change `docs/scheduler-table.md` and re-run `tools/embed-scheduler.sh`.
- [ ] **Step 4:** Write `docs/engaging-verification.md`: the commands the user runs on Engaging (module check; copy `examples/mit-orcd-engaging.conf` to a scratch conf; run `/nfcore-rnaseq-setup` or `/bulk-rnaseq-pipeline` on a tiny input; submit with the printed command; what to paste back).
- [ ] **Step 5:** Run `tests/check_scheduler.sh` and `tests/test_embed_scheduler.sh`. Expected: both PASS.
- [ ] **Step 6:** Commit `Document schedulers; add Engaging example`.

### Task 8: Full regression

- [ ] **Step 1:** Run `check_skill.sh` of rnavar, rnasplice-setup, rnasplice-report, ase-pipeline (`ASE_SCRATCH_ROOT` set) on the login node. Expected: four `PASS`.
- [ ] **Step 2:** Through `sbatch` (scratchpad job scripts, `-p` from the test conf): `tests/run_all_tests.sh` of rnasplice-setup and rnasplice-report, and `prove_mutations.sh` of both. Expected: `ALL TESTS PASS`, `MUTATIONS PASS`.
- [ ] **Step 3:** Report to the user: pass counts, the 97 `rtest` rows and the PBS/SGE/LSF real-cluster runs not performed. Wait for the user before pushing.

### Task 9: Engaging verification (user-run)

- [ ] **Step 1:** Give the user `docs/engaging-verification.md`; they run it on Engaging.
- [ ] **Step 2:** From the pasted output, fill `SINGULARITY_LOAD` in `examples/mit-orcd-engaging.conf` (or keep it empty with the note), fix anything the run exposed in the same way as Tasks 4–6 (test first), re-run Task 8 step 1, commit `Engaging: verified example config`.
