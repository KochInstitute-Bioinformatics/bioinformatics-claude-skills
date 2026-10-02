#!/bin/bash
# Static checks for nfcore-rnasplice-setup.md. Usage: check_skill.sh <skill.md>
# Uses only the recorded fixtures of the gated pipeline revision (tests/fixtures/); no network, no python.
# The README checked is the README.md next to the skill file. CHECK_LIST_NEEDS=1 also prints "NEED <line> <text>" per need.
# Exit codes: 0 PASS, 1 a check failed, 2 a fixture is missing or invalid (gate values are validated against their allowed sets).
# Call need/needr at top level (also inside a loop, `case` or `||`), never from inside a helper function: BASH_LINENO[0]
# would then be the line inside the helper, and prove_red.sh could not map the need to an added checker line.
set -u
SKILL=${1:?usage: check_skill.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd); FIX="$HERE/fixtures"; CUT="$HERE/cut_block.sh"
README="$(dirname "$SKILL")/README.md"
fail=0
[ -s "$SKILL" ] || { echo "FAIL: skill file missing or empty: $SKILL"; exit 1; }
for f in rnasplice_schema.json config_params.txt gate_values.tsv trace_process_names.txt output_tree.txt test_samplesheet.csv test_contrastsheet.csv schema_input.json schema_input_genome_bam.json; do
  [ -s "$FIX/$f" ] || { echo "cannot read fixture $FIX/$f"; exit 2; }
done
GATE_KEYS="GATE_OUTCOME PIPELINE_REVISION VERSION_TAG NEXTFLOW_TESTED NEXTFLOW_MIN NEXTFLOW_MAX_EXCL CONDA_ENV_TESTED HAS_MAX_PARAMS HAS_RESOURCE_LIMITS SALMON_ROUTE PSEUDO_OFF_LINE BAM_SHEET_HEADER BAM_RMATS_LIBTYPE BAM_RMATS_READTYPE BAM_DEXSEQ_STRAND BAM_FC_STRAND BAM_NEEDS_BAI RMATS_BAMLIST_ORDER RMATS_B1_GROUP STAR_VERSION_GENOME SALMON_INDEX_VERSION DTU_FILTER_SCOPE TEST_CONTRAST TEST_READ_LENGTH EDGER_DEU_FUNCTION"
for k in $GATE_KEYS; do
  [ "$(awk -F'\t' -v k="$k" '$1 == k && $2 != ""' "$FIX/gate_values.tsv" | wc -l)" -eq 1 ] \
    || { echo "gate value $k must appear exactly once, non-empty, in fixtures/gate_values.tsv"; exit 2; }
done
gv() { awk -F'\t' -v k="$1" '$1 == k {print $2; exit}' "$FIX/gate_values.tsv"; }
# Allowed values of every gate key (Task 0 brief table); anything else is a broken fixture (exit 2).
gate_ok() { local v; v=$(gv "$1"); [[ $v =~ $2 ]] || { echo "gate value $1='$v' is not allowed (expected $2)"; exit 2; }; }
SEMVER='^[0-9]+\.[0-9]+\.[0-9]+$'
gate_ok GATE_OUTCOME '^(A|B|C)$'
gate_ok PIPELINE_REVISION '^([0-9a-f]{40}|[0-9]+\.[0-9]+\.[0-9]+)$'
gate_ok NEXTFLOW_TESTED "$SEMVER"
gate_ok NEXTFLOW_MIN "$SEMVER"
gate_ok NEXTFLOW_MAX_EXCL '^(none|[0-9]+\.[0-9]+\.[0-9]+)$'
gate_ok CONDA_ENV_TESTED '^[A-Za-z0-9][A-Za-z0-9._-]*$'
gate_ok HAS_MAX_PARAMS '^(yes|no)$'
gate_ok HAS_RESOURCE_LIMITS '^(yes|no)$'
gate_ok SALMON_ROUTE '^(pseudo_only|star_salmon_only|star_salmon_both)$'
gate_ok PSEUDO_OFF_LINE '^(none|[a-z_0-9]+: .+)$'
gate_ok BAM_SHEET_HEADER '^sample,condition,genome_bam(,[a-z_0-9]+)*$'
gate_ok BAM_RMATS_LIBTYPE '^(fr-unstranded|fr-firststrand|fr-secondstrand|from_sheet|none)$'
gate_ok BAM_RMATS_READTYPE '^(paired|single|from_sheet|none)$'
gate_ok BAM_DEXSEQ_STRAND '^(no|yes|reverse|from_sheet|none)$'
gate_ok BAM_FC_STRAND '^(0|1|2|from_sheet|none)$'
gate_ok BAM_NEEDS_BAI '^(yes|no|n/a)$'
gate_ok RMATS_BAMLIST_ORDER '^(sheet|sorted_by_name|unordered)$'
gate_ok RMATS_B1_GROUP '^(treatment|control)$'
gate_ok STAR_VERSION_GENOME '^[0-9]+\.[0-9]+\.[0-9]+[a-z]?$'
gate_ok SALMON_INDEX_VERSION '^[0-9]+$'
gate_ok DTU_FILTER_SCOPE '^(all_samples|per_contrast)$'
gate_ok TEST_CONTRAST '^[A-Za-z][A-Za-z0-9_]*_vs_[A-Za-z][A-Za-z0-9_]*$'
gate_ok TEST_READ_LENGTH '^[1-9][0-9]*$'
gate_ok EDGER_DEU_FUNCTION '^(diffSpliceDGE|diffSplice|exactTest|glmQLFTest|glmQLFit|glmLRT|glmFit)( (diffSpliceDGE|diffSplice|exactTest|glmQLFTest|glmQLFit|glmLRT|glmFit))*$'
# Branch consistency: B = a pinned commit, tag dev-<7>, no upper Nextflow bound; A/C = a release, tag = release; C has an upper bound.
case "$(gv GATE_OUTCOME)" in
  B) gate_ok PIPELINE_REVISION '^[0-9a-f]{40}$'; gate_ok VERSION_TAG "^dev-$(gv PIPELINE_REVISION | cut -c1-7)\$"; gate_ok NEXTFLOW_MAX_EXCL '^none$' ;;
  A) gate_ok PIPELINE_REVISION "$SEMVER"; gate_ok VERSION_TAG "^$(gv PIPELINE_REVISION | sed 's/\./\\./g')\$"; gate_ok NEXTFLOW_MAX_EXCL '^none$' ;;
  C) gate_ok PIPELINE_REVISION "$SEMVER"; gate_ok VERSION_TAG "^$(gv PIPELINE_REVISION | sed 's/\./\\./g')\$"; gate_ok NEXTFLOW_MAX_EXCL "$SEMVER" ;;
esac
[ "$(gv SALMON_ROUTE)" = star_salmon_only ] || gate_ok PSEUDO_OFF_LINE '^none$'

need()   { [ -n "${CHECK_LIST_NEEDS:-}" ] && echo "NEED ${BASH_LINENO[0]} $1"; grep -qF -- "$1" "$SKILL" || { echo "FAIL: missing required text: $1"; fail=1; }; }
needr()  { [ -n "${CHECK_LIST_NEEDS:-}" ] && echo "NEED ${BASH_LINENO[0]} $1"; grep -qF -- "$1" "$README" 2>/dev/null || { echo "FAIL: README missing required text: $1"; fail=1; }; }
forbid() { ! grep -qF -- "$1" "$SKILL" || { echo "FAIL: forbidden text present: $1"; fail=1; }; }
anchor_once() { [ "$(grep -cF -- "$1" "$SKILL")" -eq 1 ] || { echo "FAIL: anchor must appear exactly once: $1"; fail=1; }; }
# forbid_re <extended regex> <label>: no line of the skill may match the regex.
forbid_re() { ! grep -qE -- "$1" "$SKILL" || { echo "FAIL: forbidden pattern present: $2"; fail=1; }; }

# Schema parameters = keys whose parent object key is "properties" (group names and structural keys are excluded).
schema_names=$(awk '
  { s = $0
    while (s != "") {
      if (match(s, /^[ \t]*"[^"]*"[ \t]*:[ \t]*\{/)) {
        k = s; sub(/^[ \t]*"/, "", k); sub(/".*$/, "", k)
        if (d > 0 && st[d] == "properties") print k
        st[++d] = k; s = substr(s, RSTART + RLENGTH); continue
      }
      c = substr(s, 1, 1); s = substr(s, 2)
      if (c == "\"") { while (s != "") { c2 = substr(s, 1, 1); s = substr(s, 2); if (c2 == "\\") { s = substr(s, 2); continue }; if (c2 == "\"") break } }
      else if (c == "{") st[++d] = ""
      else if (c == "}") d--
    } }' "$FIX/rnasplice_schema.json" | sort -u)
for p in input contrasts source outdir rmats rmats_read_len; do
  echo "$schema_names" | grep -qx "$p" || { echo "schema extraction broken: $p missing"; exit 2; }
done
for p in properties definitions input_output_options; do
  ! echo "$schema_names" | grep -qx "$p" || { echo "schema extraction broken: structural key $p listed"; exit 2; }
done
# Recorded pipeline config defaults: key<TAB>value with surrounding quotes and trailing // comments removed.
config_kv=$(awk -v q="'" '/^params[ \t]*\{/ {f = 1; next} f && /^\}/ {exit}
  f && /^[ \t]*[A-Za-z_0-9]+[ \t]*=/ { k = $0; sub(/^[ \t]*/, "", k); sub(/[ \t]*=.*$/, "", k)
    v = $0; sub(/^[^=]*=[ \t]*/, "", v); sub(/[ \t]+\/\/.*$/, "", v); sub(/[ \t]+$/, "", v)
    gsub("^[\"" q "]|[\"" q "]$", "", v); print k "\t" v }' "$FIX/config_params.txt")
echo "$config_kv" | grep -q "^rmats	" || { echo "config extraction broken: rmats missing"; exit 2; }

# (a) no pipeline parameter as a --flag anywhere; only allowlisted tool flags (never a schema name)
allow="mail-type mail-user mem dependency parsable variable-read-length allow-clipping"
for a in $allow; do ! echo "$schema_names" | grep -qx -- "$a" || { echo "allowlist contains a pipeline parameter: $a"; exit 2; }; done
for flag in $(grep -oE '(^|[^A-Za-z0-9_-])--[A-Za-z_][A-Za-z_0-9-]*' "$SKILL" | sed -E 's/^[^-]*--//' | sort -u); do
  echo "$allow" | tr ' ' '\n' | grep -qx -- "$flag" || { echo "FAIL: --$flag is not an allowlisted tool flag (pipeline parameters go in the params file, never as --flags)"; fail=1; }
done
# (b) every key in a fenced yaml block is a parameter of the recorded schema. Fences are 3 or more backticks; an info
# string yaml or yml (any case) marks a yaml block; a bare fence at least as long as the open one closes it; any other
# fence line inside a block opens a nested block (stricter than CommonMark on purpose, so nested yaml examples are checked).
while IFS= read -r k; do
  echo "$schema_names" | grep -qx -- "$k" || { echo "FAIL: params key not in the recorded schema: $k"; fail=1; }
done < <(awk '
  match($0, /^`+/) && RLENGTH >= 3 {
    n = RLENGTH; info = tolower(substr($0, n + 1)); gsub(/^[ \t]+|[ \t]+$/, "", info)
    if (d > 0 && info == "" && n >= len[d]) { d--; next }
    len[++d] = n; y[d] = (info == "yaml" || info == "yml"); next
  }
  d > 0 && y[d]' "$SKILL" | grep -oE '^[A-Za-z_0-9-]+:' | tr -d ':' | sort -u)
forbid "«"
forbid "»"
forbid "python3"
forbid "python -c"
forbid_re '(^|[^A-Za-z0-9_.-])gh[[:space:]]+[a-z]' "a gh command (gh is not installed; use WebFetch)"
forbid_re '(^|[[:space:]])\[(A/C|A/B|B|C)\]([[:space:]]|$)' "a plan variant marker [A/C], [A/B], [B] or [C]"
forbid "nextflow pull"

# (c) content checks, appended per task -------------------------------
# --- Task 1 (Steps 0-3)
need "# nf-core/rnasplice Pipeline Setup Skill"
for n in 0 1 2 3; do need "## Step $n — "; done
need "{TODAY_ISO}"
need "api.github.com/repos/nf-core/rnasplice/releases/latest"
need "Do NOT use the \`gh\` CLI"
forbid "gh api"
need "verified against nf-core/rnasplice $(gv PIPELINE_REVISION) with Nextflow $(gv NEXTFLOW_TESTED)"
need "Nextflow $(gv NEXTFLOW_MIN) or newer"
[ "$(gv NEXTFLOW_MAX_EXCL)" = none ] || need "and older than $(gv NEXTFLOW_MAX_EXCL)"
need "Set \`{VERSION_TAG}\`"
need "Do not run \`nextflow\` here"
case "$(gv GATE_OUTCOME)" in
  B) need "that is a commit of the development branch (\`$(gv VERSION_TAG)\`)" ;;
  C) need "create_nextflow_env_$(gv NEXTFLOW_TESTED).sh" ;;
esac
# Steps 0-2 rules (as in nfcore-rnavar-setup)
need "run \`pwd\` to record the current working directory (\`{CWD}\`)"
need "\`{WD_NAME}\` = basename of \`{CWD}\`"
need "\`{TODAY_YYMMDD}\` = today's date as \`YYMMDD\`"
need "\`{TODAY_ISO}\` = today's date as \`YYYY-MM-DD\`"
need "Do **not** pre-suggest or pre-fill any email address. Wait for the answer before Step 2."
need "Do **not** pre-suggest any environment name. Ask this as a separate question after Step 1"
need "never combine Steps 1 and 2 in one message"
# Step 3: the pinned revision is the default; a newer release only after its schema has every key this skill writes
need "raw.githubusercontent.com/nf-core/rnasplice/{TAG}/nextflow_schema.json"
need "(the Step 8 module block and the Step 11 template)"
case "$(gv GATE_OUTCOME)" in
  B) need "- Set \`{VERSION}\` = $(gv PIPELINE_REVISION) (\`nextflow run -r\` accepts a commit)."
     need "1. Use the pinned commit (verified; default)"
     need "- If the fetch fails (no network, rate limit, or no \`tag_name\` in the answer): tell the user that the latest release could not be checked and keep the pinned commit."
     need "- If \`{TAG}\` is 1.0.4: tell the user that the pinned commit is used"
     need "were verified for the pinned commit only, not for release {TAG}" ;;
esac
# --- end Task 1

# --- Task 2 (Step 4)
need "## Step 4 — Input files and samplesheet rows"
need "sample,fastq_1,fastq_2,strandedness,condition"
need "$(head -n1 "$FIX/test_samplesheet.csv" | tr -d '\r')"
need "Strandedness is asked, never guessed."
need "1. unstranded · 2. forward · 3. reverse · 4. I don't know"
need "rnasplice has no \`auto\` strandedness"
forbid "strandedness to \`auto\`"
forbid "strandedness: auto"
need "Never continue with \"I don't know\""
need "rMATS needs one strandedness and one read type for all samples"
need "if a name then does not start with a letter (a digit or \`_\`), prefix \`S\`"
need "if it is an R reserved word, append \`_S\`"
need "every other character outside \`A-Za-z0-9_\`"
need "1. the same sample (lanes or technical replicates; the pipeline merges their reads)"
need "2. different samples — rename them"
anchor_once "**Strandedness from an existing nf-core/rnaseq run.**"
anchor_once "**BAM input rule.**"
need "infer_strandedness() {"
need "bam_input_allowed() {"
need "RMATS_LIBTYPE=\"$(gv BAM_RMATS_LIBTYPE)\""
need "RMATS_READTYPE=\"$(gv BAM_RMATS_READTYPE)\""
need "DEXSEQ_STRAND=\"$(gv BAM_DEXSEQ_STRAND)\""
need "FC_STRAND=\"$(gv BAM_FC_STRAND)\""
need "1. Start from the FASTQ files instead (go to Step 4a)"
need "DTU and SUPPA2 need Salmon quantification from reads"
[ "$(gv BAM_RMATS_LIBTYPE)" = none ] || need "$(gv BAM_SHEET_HEADER)"
# Controller ruling (b): a BAM sheet with strandedness/single_end columns carries the asked values in every row.
case "$(gv BAM_SHEET_HEADER)" in *,strandedness,single_end*)
  need "write \`{STRANDEDNESS}\` and \`true\` (single-end) or \`false\` in every row" ;; esac
# Task 2 review fixes (I1-I3, minors)
need "*reverse* = read 1 comes from the strand opposite to the transcript"
need "*forward* = read 1 comes from the transcript strand"
need "do not continue until they answer 1, 2 or 3"
need "so Step 8 switches them off"
need "if the directory mixes paired-end and single-end files, stop"
need "\`_R1.\`/\`_R2.\`"
need "unstranded if the two fractions differ by at most 0.1"
need "No \`*.infer_experiment.txt\` file"
forbid_re 'infer_experiment\.txt.*\|[[:space:]]*head' "head on the RSeQC file list (every file must be checked)"
need "-path \"*/work\" -prune"
need "\`{BAM_DIR}\`"
forbid "{DIR}"
need "the pipeline re-sorts and indexes them itself"
forbid "be coordinate-sorted"
need "prefer \`.markdup.sorted.bam\`"
anchor_once "**Sample names and read pairs.**"
need "fastq_r2_name() {"
need "fastq_sample_name() {"
need "sanitize_sample_name() {"
need "sample_name_collisions() {"
need "input_path_ok() {"
# The sample rule of the pinned samplesheet schemas (FASTQ and genome BAM) is quoted in the skill.
samp_pat=$(awk '/"sample": \{/ {f = 1} f && /"pattern"/ {sub(/^[ \t]*"pattern": "/, ""); sub(/",[ \t]*$/, ""); print; exit}' "$FIX/schema_input.json" | sed 's/\\\\/\\/g')
[ -n "$samp_pat" ] || { echo "sample pattern extraction broken: fixtures/schema_input.json"; exit 2; }
need "\`$samp_pat\`"
# --- end Task 2

[ $fail -eq 0 ] && echo "PASS" || exit 1
