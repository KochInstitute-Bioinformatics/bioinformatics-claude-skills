# Checking the skills on MIT ORCD Engaging

Run these on an Engaging login node (light commands only), then paste the output back.

1. Find the Apptainer/Singularity module name:
   ```
   module avail apptainer
   module avail singularity
   ```
2. Copy the example config and put the module line from step 1 into `SINGULARITY_LOAD`:
   ```
   cp examples/mit-orcd-engaging.conf ~/.claude/bioinfo-skills.conf
   ```
3. In a small project directory run `/bulk-rnaseq-pipeline` (needs nf-core/rnaseq `star_salmon` output) or `/nfcore-rnaseq-setup` on tiny FASTQ files, and read every generated script. Check: `#SBATCH -p mit_normal`, no walltime above 12:00:00, the Singularity load line you gave.
4. Submit with the command the skill prints and wait for the job to finish: `squeue -u $USER`, then the log.
5. Paste back: the module output from step 1, the first 15 lines of one generated job script, and the last 20 lines of its log.
