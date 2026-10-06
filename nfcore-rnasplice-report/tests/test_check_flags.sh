#!/bin/bash
# Test of the --flag rule of check_skill.sh on small skills built in a temp dir: prose outside ``` fences must never contain a
# --flag; inside a fence a --flag is allowed only in a command segment (the line split on ; && || and |) that starts with (or continues, after a trailing backslash)
# singularity exec, sbatch, Rscript or srun, and on #SBATCH lines. Only the flag message of the checker is looked at.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
T=$(mktemp -d /tmp/tmp.XXXXXXXXXX) || exit 1
trap 'rm -f "$T/nfcore-rnasplice-report.md" "$T/README.md" "$T/out" "$T/root.md"; rmdir "$T"' EXIT
case "$T" in /tmp/tmp.??????????) ;; *) echo "unexpected temp dir: $T"; exit 1 ;; esac
: > "$T/README.md"
bad=0
# flagcase <name> <expected: red|ok> <body text after the title line>
flagcase() {
  { echo "# nf-core/rnasplice Report Skill"; printf '%s\n' "$3"; } > "$T/nfcore-rnasplice-report.md"
  bash "$HERE/check_skill.sh" "$T/nfcore-rnasplice-report.md" > "$T/out" 2>&1
  if grep -q '^FAIL: a --flag appears' "$T/out"; then got=red; else got=ok; fi
  if [ "$got" = "$2" ]; then echo "ok   $1"; else echo "FAIL $1: expected $2, got $got"; bad=1; fi
}
flagcase prose-flag red 'Run `report --now` here.'
flagcase prose-flag-plain red 'Call report --now here.'
flagcase fenced-singularity ok $'```bash\nsingularity exec --bind /x img true\n```'
flagcase fenced-sbatch ok $'```bash\nsbatch --parsable job.sh\n```'
flagcase fenced-rscript ok $'```bash\nRscript --vanilla x.R\n```'
flagcase fenced-srun ok $'```bash\nsrun --mem=4G true\n```'
flagcase fenced-sbatch-header ok $'```bash\n#SBATCH --mem=4G\n```'
flagcase fenced-continuation ok $'```bash\nsingularity exec \\\n  --bind /x img true\n```'
flagcase fenced-ls red $'```bash\nls --color\n```'
flagcase fenced-continuation-of-ls red $'```bash\nls \\\n  --color\n```'
flagcase fenced-assignment-red red $'```bash\nout=$(ls --color)\n```'
flagcase capture-sbatch ok $'```bash\njid=$(sbatch --parsable job.sh)\n```'
flagcase capture-sbatch-dep ok $'```bash\nJOB02=$(sbatch --parsable --dependency=afterok:$JOB01 run_02.sbatch) || exit 1\n```'
flagcase capture-placeholder ok $'```\nJOB__N__=$(sbatch --parsable --dependency=afterok:$JOB01 run___N__.sbatch) || exit 1\n```'
flagcase capture-ls red $'```bash\njid=$(ls --color)\n```'
flagcase segment-semicolon red $'```bash\nsbatch x; ls --color\n```'
flagcase segment-and red $'```bash\nsbatch x && rm --force y\n```'
flagcase segment-or red $'```bash\nsbatch x || ls --color\n```'
flagcase segment-pipe red $'```bash\nsbatch x | tee --append y\n```'
flagcase segment-in-capture red $'```bash\njid=$(sbatch x; ls --color)\n```'
flagcase segment-continuation-red red $'```bash\nsbatch x; ls \\\n  --color\n```'
flagcase segment-all-allowed ok $'```bash\nsbatch --parsable a.sh && sbatch --parsable b.sh\n```'
flagcase flag-after-fence-closed red $'```bash\ntrue\n```\nRun report --now here.'
# Test hooks CHECK_ROOT_README and CHECK_RUN_ALL: honoured only for a file under a scratch directory made by mktemp (/tmp/tmp.XXXXXXXXXX,
# as prove_mutations.sh makes it); any other path fails the run. A hook that is honoured prints a NOTE line.
hookcase() { # name var path expect(red|note)
  { echo "# nf-core/rnasplice Report Skill"; } > "$T/nfcore-rnasplice-report.md"; echo row > "$T/root.md"
  env -i PATH="$PATH" "$2=$3" bash "$HERE/check_skill.sh" "$T/nfcore-rnasplice-report.md" > "$T/out" 2>&1
  if [ "$4" = red ]; then grep -qF "FAIL: test hook $2=$3 is not under a scratch directory" "$T/out" && r=ok || r=bad
  else grep -qF "NOTE: test hook $2 is set ($3)" "$T/out" && ! grep -qF 'is not under a scratch directory' "$T/out" && r=ok || r=bad; fi
  if [ $r = ok ]; then echo "ok   $1"; else echo "FAIL $1"; bad=1; fi
}
hookcase hook-root-stray CHECK_ROOT_README "$HERE/fixtures/sign_rows.tsv" red
hookcase hook-runall-stray CHECK_RUN_ALL "$HERE/run_all_tests.sh" red
hookcase hook-root-scratch CHECK_ROOT_README "$T/root.md" note
hookcase hook-runall-scratch CHECK_RUN_ALL "$T/root.md" note
[ $bad -eq 0 ] && echo "FLAG RULE PASS" || exit 1
