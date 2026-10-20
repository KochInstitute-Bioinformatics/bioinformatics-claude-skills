# Scheduler table (source of truth)

`tools/embed-scheduler.sh` copies the two marked regions below, verbatim, into the "Site configuration" section of every job-writing skill. Edit the regions here, then run `tools/embed-scheduler.sh */*.md` (the skills that have no markers are skipped by the list you pass). `tests/check_scheduler.sh` fails when a skill's copy differs.

<!-- SCHEDULER:BEGIN rules -->
**Scheduler.** Besides the keys above, read `SCHEDULER` (`slurm`, `pbs`, `sge`, `lsf` or `local`), `QUEUE`, `QUEUE_ON_CLI`, `MAX_TIME`, `MAX_MEM_GB`, `MAX_CPUS` and `PBS_FLAVOR` (`pbspro`, the default, or `torque`). `SLURM_PARTITION` and `SLURM_PARTITION_ON_CLI` are older names: read them as `QUEUE` and `QUEUE_ON_CLI` when the new names are absent.
- Trim and lowercase `SCHEDULER`. When it is missing, or is not one of the five values, ask (numbered: 1. slurm · 2. pbs · 3. sge · 4. lsf · 5. local). Never assume `slurm`.
- The scripts in this skill are written in SLURM form. For any other scheduler, rewrite every `#SBATCH` line and every `sbatch`, `squeue`, `scancel`, `sacct`, `srun` or `salloc` command through the table below, converting units as the table says. When a construct has no row, stop and say `No row for <construct> in the scheduler table`; never drop an option silently.
- Clamp every request to `MAX_TIME`, `MAX_MEM_GB` and `MAX_CPUS` when they are set (compare times after converting `D-HH:MM:SS` to seconds). A job that must run longer than `MAX_TIME` (a Nextflow head job) is requested with `MAX_TIME`, and you tell the user to resubmit with `-resume` or to pick a longer queue.
- `local`: write no header lines and no scheduler commands, list the scripts in the order to run them (`bash a.sh && bash b.sh`), and for a Nextflow head job recommend `nohup bash script.sh > script.log 2>&1 &` or `tmux`. The Nextflow executor is `local`.
- `{PCLI_FLAG}` is the queue option of the submit command (`-p Q ` for `sbatch`, `-q Q ` for `qsub` and `bsub`); it is empty unless `QUEUE_ON_CLI` is `yes` and `QUEUE` is not empty. An empty `QUEUE` means no queue line and no queue option on any scheduler.
<!-- SCHEDULER:END rules -->

<!-- SCHEDULER:BEGIN table -->
| SLURM construct | pbs | sge | lsf | local |
|---|---|---|---|---|
| `#SBATCH -J NAME` | `#PBS -N NAME` | `#$ -N NAME` | `#BSUB -J NAME` | omit |
| `#SBATCH -p Q` | `#PBS -q Q` | `#$ -q Q` | `#BSUB -q Q` | omit |
| `#SBATCH -N 1` | part of the `select` line | omit | omit | omit |
| `#SBATCH -n N` or `-c N` | `#PBS -l select=1:ncpus=N:mem=Mgb` (one line with the memory) | `#$ -pe smp N` | `#BSUB -n N` and `#BSUB -R "span[hosts=1]"` | omit |
| `#SBATCH --mem=MG` | `mem=Mgb` in the `select` line | `#$ -l h_vmem=XG` where X = M / N (per slot) | `#BSUB -M MB` and `#BSUB -R "rusage[mem=MB]"` where MB = M * 1024 (the unit is a site setting; say so) | omit |
| `#SBATCH -t D-HH:MM:SS` | `#PBS -l walltime=HH:MM:SS` (hours = 24 * D + HH) | `#$ -l h_rt=HH:MM:SS` | `#BSUB -W HH:MM` | omit |
| `#SBATCH -o FILE` (`%j` in FILE is the job id) | `#PBS -o FILE` and `#PBS -j oe` (drop `%j`) | `#$ -o FILE` and `#$ -j y` (`%j` becomes `$JOB_ID`) | `#BSUB -o FILE` (`%j` becomes `%J`) | run as `bash script.sh > FILE 2>&1` (drop `%j`) |
| `#SBATCH --mail-type=END,FAIL` with `--mail-user=A` | `#PBS -m ae` and `#PBS -M A` | `#$ -m ea` and `#$ -M A` | `#BSUB -N` and `#BSUB -u A` | omit |
| `sbatch SCRIPT` | `qsub SCRIPT` | `qsub SCRIPT` | `bsub < SCRIPT` | `bash SCRIPT` |
| `sbatch --parsable SCRIPT` (prints the job id) | `qsub SCRIPT` (prints the id) | `qsub -terse SCRIPT` | `bsub < SCRIPT`, id from `sed -n 's/^Job <\([0-9]*\)>.*/\1/p'` | not applicable |
| `sbatch --dependency=afterok:ID SCRIPT` | `qsub -W depend=afterok:ID SCRIPT` | `qsub -hold_jid ID SCRIPT` | `bsub -w "done(ID)" < SCRIPT` | run after the previous script succeeded (`&&`) |
| `sbatch --dependency=afterany:ID1:ID2 SCRIPT` (runs when the jobs ended, whatever their status) | `qsub -W depend=afterany:ID1:ID2 SCRIPT` | `qsub -hold_jid ID1,ID2 SCRIPT` | `bsub -w "ended(ID1) && ended(ID2)" < SCRIPT` | run after the previous scripts, whatever their status (`;`) |
| `--kill-on-invalid-dep=yes` | omit | omit | omit | omit |
| `squeue -u $USER` | `qstat -u $USER` | `qstat -u $USER` | `bjobs` | `jobs` |
| `scancel ID` | `qdel ID` | `qdel ID` | `bkill ID` | `kill` |
| `sacct` | `qstat -x` | `qacct -j ID` | `bhist -l ID` | not applicable |
| `srun` or `salloc` (debugging) | `qsub -I` | `qrsh` | `bsub -Is bash` | run directly |
| Nextflow `executor = 'slurm'` and `-profile slurm` | `pbspro` (`pbs` when `PBS_FLAVOR=torque`) | `sge` | `lsf` | `local` |
<!-- SCHEDULER:END table -->
