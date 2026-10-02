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
# SALMON_INDEX_VERSION stays a validated fixture value, but the skill no longer uses it: by controller ruling (Task 4 review I3)
# the skill never reuses a Salmon index (the pipeline always builds its own from the GTF-derived transcripts).
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

# --- Task 3 (Step 5)
need "## Step 5 — Conditions, contrasts and the two sheets"
need "contrast,treatment,control"
need "$(head -n1 "$FIX/test_contrastsheet.csv" | tr -d '\r')"
need "{treatment}_vs_{control}"
need "1. All pairwise comparisons · 2. Only the comparisons I list"
need "1. No, the samples are independent (default) · 2. Yes, every sample has a partner"
# The paired-test defaults come from the recorded config (at the gate rmats_paired_stats = false, diffsplice_paired = true).
cfg_default() { echo "$config_kv" | awk -F'\t' -v k="$1" '$1 == k {print $2; exit}'; }
[ -n "$(cfg_default rmats_paired_stats)" ] && [ -n "$(cfg_default diffsplice_paired)" ] || { echo "config extraction broken: rmats_paired_stats or diffsplice_paired missing"; exit 2; }
[ "$(cfg_default diffsplice_paired)" = true ] && need "the pipeline default \`diffsplice_paired: true\` assumes a pairing"
[ "$(cfg_default rmats_paired_stats)" = true ] && need "the pipeline default \`rmats_paired_stats: true\` assumes a pairing"
need "\`rmats_paired_stats\` defaults to \`$(cfg_default rmats_paired_stats)\` at this revision"
[ "$(cfg_default rmats_paired_stats)" = true ] || forbid "\`rmats_paired_stats: true\` assumes"
need "Every condition used in a contrast needs at least 2 samples"
anchor_once "**Sheet validation.**"
need "validate_rnasplice_sheets() {"
need "bamhdr=\"$(gv BAM_SHEET_HEADER)\""
# Step 5 repeats the BAM header (validator, procedure): pin the Step 4b sentence itself so T2-bamhdr still bites.
[ "$(gv BAM_RMATS_LIBTYPE)" = none ] || need "The samplesheet header is \`$(gv BAM_SHEET_HEADER)\`"
need "{SHEET_PREFIX}_samplesheet.csv"
need "{SHEET_PREFIX}_contrasts.csv"
need "never edit the validator"
need "1. overwrite · 2. choose another filename"
forbid "rm -rf"
case "$(gv RMATS_BAMLIST_ORDER)" in
  sheet) need "write the rows ordered by subject within each condition"
         need "validate_rnasplice_sheets \"\$T/samplesheet.csv\" \"\$T/contrasts.csv\" \"\$src\" 1 \"\$T/pairs.csv\" || rc=1"
         need "the rows of one sample (technical replicates) stay next to each other" ;;
  sorted_by_name) need "propose names \`{SUBJECT}_{CONDITION}\`" ;;
  unordered) need "a paired rMATS test cannot be set up safely" ;;
esac
need "Ask only when there are exactly two conditions with the same number of samples"
need "sets both \`rmats_paired_stats\` and SUPPA2's \`diffsplice_paired\` in Step 8"
need "treatment = the later and control = the earlier one in \`{CONDITIONS}\`"
# --- end Task 3

# --- Task 3 review fixes
# I1: the key Step 5 decisions
need "Answer 1: \`{PAIRED_DESIGN}\` = \`false\`"
need "this skill writes both \`false\` unless the user confirms pairing"
need "otherwise set \`{PAIRED_DESIGN}\` = \`false\` without asking"
need "validate_rnasplice_sheets \"\$T/samplesheet.csv\" \"\$T/contrasts.csv\" \"\$src\" 0 || rc=1"
need "ask for an order with the reference first"
need "validated there, and only then moved into \`{CWD}\`"
# I2: shell variables and functions do not persist between Bash tool calls
need "**Bash tool calls.** Shell variables and functions do not persist between Bash tool calls."
need "**Writing the sheets.** Shell variables and functions do not persist between Bash tool calls, so the whole procedure is one Bash call"
need "paste the whole block at the start of every Bash call that uses its functions"
need "In one Bash call, define the function below and run it on every one of those files"
need "In one Bash call, define this function and run \`bam_input_allowed {STRANDEDNESS} {LAYOUT}\`"
need "install_rnasplice_sheets() {"
need "install_rnasplice_sheets \"\$T\" {SOURCE} 0 {SAMPLESHEET_CSV} {CONTRASTS_CSV}"
need "install_rnasplice_sheets \"\$T\" {SOURCE} 1 {SAMPLESHEET_CSV} {CONTRASTS_CSV} {SHEET_PREFIX}_pairs.csv"
need "The scratch files and directory are removed on every path"
forbid_re 'in the Bash tool' "'in the Bash tool' (a function defined in one Bash call cannot be used in another)"
# $T (the scratch directory) may appear only inside fenced code blocks: a prose step that uses it runs in another call.
t_outside=$(awk '/^```/ { inb = !inb; next } !inb && /\$T([^A-Za-z0-9_]|$)/ { print FNR }' "$SKILL" | tr '\n' ' ')
[ -z "$t_outside" ] || { echo "FAIL: \$T used outside a code block (line $t_outside): shell variables do not persist between Bash tool calls"; fail=1; }
# Minors
need "start with a letter and are not R reserved words"
need "tell the user the direction"
[ "$(gv RMATS_B1_GROUP)" = treatment ] && need "a positive rMATS IncLevelDifference means more inclusion in the treatment"
need "For BAM input, \`{SEQ_DATE}\` ="
need "store the new name in \`{SAMPLESHEET_CSV}\` or \`{CONTRASTS_CSV}\`"
need "A missing or empty input file stops the call"
# --- end Task 3 review fixes

# --- Task 4 (Steps 6-7)
need "## Step 6 — Read length for rMATS"
anchor_once "**Read-length detection.**"
need "detect_read_length() {"
need "the pipeline default of 40 is wrong for almost all data"
need "\`--variable-read-length\` and \`--allow-clipping\`"
need "an integer between 20 and 1000"
need "## Step 7 — Organism and genome files"
need "{genome_base}/{organism}/{assembly}_ens{version}/"
need "Ask the base directory only for option 1"
need "custom_{WD_NAME}"
need "versionGenome"
need "$(gv STAR_VERSION_GENOME)"
need "1. Reuse (default) · 2. Let the pipeline build its own"
need "\`gencode: false\` is written either way"
need "download_genome_{REF_TAG}.sh"
need "Never type a URL from memory"
# Controller rulings: one Bash call per procedure; STAR index reuse only for the gate's index format; sjdbOverhang note;
# the wizard never builds an index; pipeline-built indexes are not kept.
need "only when \`versionGenome\` is \`$(gv STAR_VERSION_GENOME)\`"
need "grep -E '^(versionGenome|sjdbOverhang)[[:space:]]' \"{STAR_DIR}/genomeParameters.txt\""
[ "$(cfg_default rmats_read_len)" = 40 ] || { echo "config extraction broken: rmats_read_len default is not 40"; exit 2; }
[ "$(cfg_default save_reference)" = false ] || { echo "config extraction broken: save_reference default is not false"; exit 2; }
need "If its \`sjdbOverhang\` is not 100"
need "the junction database is tuned for reads of"
need "The wizard never builds an index itself"
need "\`save_reference: false\`"
need "check each with a HEAD request"
# --- end Task 4

# --- Task 4 review fixes
# I1: the Ensembl release lookup (current_README is a 404); fallback listing
need "the latest Ensembl release from \`https://ftp.ensembl.org/pub/current/README\`"
need "its line \"Ensembl Release N Databases.\" gives N"
need "use the highest \`release-N/\` directory in the listing of \`https://ftp.ensembl.org/pub/\`"
forbid "pub/current_README"
# I2: Step 6-7 decisions
need "with any other value, say why and let the pipeline build it"
need "Any \`MISSING\` line: there is no usable index; let the pipeline build it."
need "\`{READ_LENGTH}\` = the number on the \`READ_LENGTH\` line (the most common R1 length over all samples), never a single sample's value"
need "Every row of the samplesheet is checked: its R1 file and, for paired-end data, its R2 file"
forbid "up to 5 different samples"
need "Version: the highest existing \`{assembly}_ens{N}\` directory unless the user asks otherwise"
need "FASTA \`Mus_musculus.GRCm39.dna.primary_assembly.fa\`"
need "Human: GRCh38, FASTA \`Homo_sapiens.GRCh38.dna.primary_assembly.fa\`"
forbid "GRCh37"
need "Check that both files exist (\`test -s\`)."
need "for f in SA SAindex Genome sjdbList.out.tab genomeParameters.txt; do"
need "(\`{STAR_DIR}\` = \`{GENOME_DIR}/index/star\` for option 1"
need "whether a STAR index built from exactly this FASTA and GTF exists: 1. No (default) · 2. Yes"
need "STAR_GENOMEGENERATE needs about 32 GB of memory and an hour or more"
need "must show 200), show them to the user,"
need "2. (numbered): 1. Ensembl release in the standard folder (default) · 2. Custom reference"
need "gene IDs with a version suffix (\`ENSG00000000003.15\`) indicate GENCODE, without one Ensembl"
need "BAM input: ask \"What is the read length of the sequencing (for example 100 or 150)?\""
forbid_re 'Ensembl release [0-9]+' "a hard-coded Ensembl release (read it at run time)"
# I3 (controller ruling): no Salmon index reuse
need "\`{SALMON_INDEX}\` is always empty: the pipeline always builds its own Salmon index from the transcripts it extracts from the GTF"
need "lacks the GTF's non-coding transcripts, which would get no quantification (DTU, SUPPA2) without any error"
need "(BAM input uses no index: skip this part; \`{STAR_INDEX}\` is empty.)"
need "FASTA, GTF and the STAR index are reused"
need "not used: rnasplice always builds its own Salmon index"
need "SALMON_INDEX: 6 CPUs, 36 GB and 8 h"
need "likely takes more than an hour (not measured by this skill's verification)"
forbid "{SALMON_DIR}"
forbid "versionInfo.json"
forbid "indexVersion"
forbid "STAR or Salmon"
# I4: every sample, grouped by condition, confounding with condition
anchor_once "**Running the read-length check.**"
need "read_length_report() {"
need "In one Bash call, paste this block, then the lines under **Running the read-length check** below."
need "After the block above, in the same Bash call"
need "Show the \`SAMPLE\` lines as a table grouped by condition"
need "- \`CONFOUNDED\`: warn loudly"
need "their PSI values, and the differences between the conditions, are biased, and rMATS gives no error"
need "1. Stop here (default) — trim all reads to one common length"
need "- \`WARNING\`: warn loudly"
# Minors: errors, informational sjdbOverhang note, miso_read_len, STAR reuse copy cost
need "never continue with an empty or guessed read length"
need "A \`gzip: stdout: Broken pipe\` message is harmless"
need "add a note (information only; it is not a reason to rebuild"
forbid "warn before asking"
need "\`miso_read_len\` belongs to the MISO sashimi plots"
forbid "rnasplice has no other read-length setting"
need "copies a given STAR index into its \`work/\` directory before aligning (about 30 GB for a human index, for every run)"
need "did not exercise the reuse path (unverified)"
# --- end Task 4 review fixes

# --- Task 5 (Steps 8-9): module block
need "## Step 8 — Analyses (modules) and their settings"
need "every analysis module is switched on, so a module that is not mentioned in the params file runs anyway"
need "Default (empty answer): 1 only"
need "This skill always writes \`sashimi_plot: false\`"
need "it is the same analysis twice"
need "## Step 9 — Trimming and QC"
anchor_once "**Module and option keys (written to the params file).**"
MOD=$(bash "$CUT" "$SKILL" "**Module and option keys (written to the params file).**" 2>/dev/null)
[ -n "$MOD" ] || { echo "FAIL: module block not found"; fail=1; }
for sw in rmats dexseq_exon edger_exon dexseq_dtu suppa sashimi_plot; do
  [ "$(printf '%s\n' "$MOD" | grep -cE "^$sw: (true|false|\{RUN_[A-Z_]+\})$")" -eq 1 ] \
    && [ "$(printf '%s\n' "$MOD" | grep -c "^$sw:")" -eq 1 ] \
    || { echo "FAIL: module switch $sw must be written exactly once with an explicit value"; fail=1; }
done
printf '%s\n' "$MOD" | grep -qx 'sashimi_plot: false' || { echo "FAIL: sashimi_plot must be written as false"; fail=1; }
for kv in 'rmats: {RUN_RMATS}' 'dexseq_exon: {RUN_DEXSEQ_EXON}' 'edger_exon: {RUN_EDGER_EXON}' 'dexseq_dtu: {RUN_DEXSEQ_DTU}' 'suppa: {RUN_SUPPA}' \
          'rmats_read_len: {READ_LENGTH}' 'rmats_paired_stats: {PAIRED_DESIGN}' 'diffsplice_paired: {PAIRED_DESIGN}' 'rmats_novel_splice_site: {RMATS_NOVEL}' \
          'min_samps_gene_expr: {MIN_SAMPS_GENE_EXPR}' 'min_samps_feature_expr: {MIN_SAMPS_FEATURE}' 'min_samps_feature_prop: {MIN_SAMPS_FEATURE}' \
          'min_gene_expr: 10' 'min_feature_expr: 10' 'min_feature_prop: 0.1' 'dtu_txi: "dtuScaledTPM"' 'rmats_splice_diff_cutoff: 0.0001'; do
  printf '%s\n' "$MOD" | grep -qxF -- "$kv" || { echo "FAIL: module block must contain the line: $kv"; fail=1; }
done
case "$(gv SALMON_ROUTE)" in
  pseudo_only) r1='aligner: "star"'; r2='pseudo_aligner: "salmon"' ;;
  star_salmon_only) r1='aligner: "star_salmon"'; r2="$(gv PSEUDO_OFF_LINE)" ;;
  star_salmon_both) r1='aligner: "star_salmon"'; r2='pseudo_aligner: "salmon"' ;;
esac
for kv in "$r1" "$r2"; do printf '%s\n' "$MOD" | grep -qxF -- "$kv" || { echo "FAIL: Salmon route $(gv SALMON_ROUTE) needs the line: $kv"; fail=1; }; done
for sw in isoformswitchanalyzer leafcutter; do
  if echo "$schema_names" | grep -qx "$sw"; then printf '%s\n' "$MOD" | grep -qx "$sw: false" || { echo "FAIL: $sw is a parameter of this revision and must be written as false"; fail=1; }; fi
done
case "$(gv DTU_FILTER_SCOPE)" in
  all_samples) need "the number of samples in the samplesheet (distinct names), because the filter is applied once to all samples" ;;
  per_contrast) need "the smallest number of samples in a contrast (treatment plus control), because the filter is applied per contrast" ;;
esac
# literal values equal the recorded pipeline config, except the deliberate deviations in DEV
# (key pattern widened from [a-z_0-9] to [A-Za-z_0-9] so that isoformswitchanalyzer_dIF is compared too)
DEV=" aligner pseudo_aligner rmats dexseq_exon edger_exon dexseq_dtu suppa sashimi_plot isoformswitchanalyzer leafcutter max_cpus max_memory max_time "
defaults_equal() {
  local line k v cv
  while IFS= read -r line; do
    k=${line%%:*}; v=${line#*: }
    case "$v" in *"{"*) continue ;; esac
    case "$DEV" in *" $k "*) continue ;; esac
    v=${v#\"}; v=${v%\"}
    cv=$(printf '%s\n' "$config_kv" | awk -F'\t' -v k="$k" '$1 == k {print $2; f = 1} END {exit !f}') \
      || { echo "FAIL: key $k has no default in the recorded pipeline config"; fail=1; continue; }
    [ "$v" = "$cv" ] || { echo "FAIL: key $k = $v differs from the pipeline config default $cv (a deliberate deviation needs a written reason and an entry in DEV)"; fail=1; }
  done < <(printf '%s\n' "$1" | grep -E '^[A-Za-z_0-9]+: ')
}
defaults_equal "$MOD"
# --- end Task 5

# --- Task 5 additions (controller rulings b, c, e; fixtures of the pinned revision)
# Ruling (e): every option key of every module group of the recorded schema is written exactly once (also for modules that are
# switched off), except keys without a default (config value null), which are never written. No key appears twice.
mod_keys=$(awk '/^        "[A-Za-z_0-9]+": \{/ { g = $1; gsub(/[":{ ]/, "", g) }
  /^                "[A-Za-z_0-9]+": \{/ { k = $1; gsub(/[":{ ]/, "", k); print g "\t" k }' "$FIX/rnasplice_schema.json")
for gk in "rmats_options	rmats_read_len" "suppa_options	clusterevents_method" "dexseq_dtu_options	min_samps_gene_expr" "miso	sashimi_plot"; do
  printf '%s\n' "$mod_keys" | grep -qxF -- "$gk" || { echo "schema group extraction broken: $gk"; exit 2; }
done
MODULE_GROUPS=" rmats_options dexseq_deu_options edger_deu_options dexseq_dtu_options miso suppa_options isoformswitchanalyzer leafcutter_options "
while IFS=$'\t' read -r g k; do
  case "$MODULE_GROUPS" in *" $g "*) ;; *) continue ;; esac
  cv=$(cfg_default "$k")
  [ -n "$cv" ] || { echo "FAIL: module option $k ($g) has no value in the recorded pipeline config"; fail=1; continue; }
  nk=$(printf '%s\n' "$MOD" | grep -c "^$k:")
  if [ "$cv" = null ]; then
    [ "$nk" -eq 0 ] || { echo "FAIL: module option $k has no default (null) and must not be written"; fail=1; }
  else
    [ "$nk" -eq 1 ] || { echo "FAIL: module option key missing from the module block or written more than once: $k"; fail=1; }
  fi
done < <(printf '%s\n' "$mod_keys")
dups=$(printf '%s\n' "$MOD" | grep -oE '^[A-Za-z_0-9]+:' | sort | uniq -d | tr '\n' ' ')
[ -z "$dups" ] || { echo "FAIL: key written more than once in the module block: $dups"; fail=1; }
# Never written by this block: salmon_index (no Salmon index reuse), genome (no iGenomes), max_* (gate HAS_MAX_PARAMS=no)
for k in salmon_index genome max_cpus max_memory max_time; do
  ! printf '%s\n' "$MOD" | grep -q "^$k:" || { echo "FAIL: the module block must not write $k"; fail=1; }
done
[ "$(gv HAS_MAX_PARAMS)" = no ] && need "\`max_cpus\`, \`max_memory\` and \`max_time\` (not parameters of this revision"
need "\`salmon_index\` (Step 7: the pipeline always builds its own Salmon index)"
need "Options without a default value"
# Leafcutter is the one module that is off in the pipeline's own config at this revision.
[ "$(cfg_default leafcutter)" = false ] && need "the one exception at this revision is LeafCutter, which is off by default"
forbid "the pipeline default is true for all of them"
# Ruling (b): BAM input switches DTU and SUPPA2 off; no Salmon runs then (gate G5: zero Salmon tasks).
need "With BAM input, options 2 and 5 are not offered (they need Salmon quantification from reads), and \`{RUN_SUPPA}\` and \`{RUN_DEXSEQ_DTU}\` are \`false\`."
need "Set \`{RUN_RMATS}\`, \`{RUN_SUPPA}\`, \`{RUN_DEXSEQ_EXON}\`, \`{RUN_EDGER_EXON}\` and \`{RUN_DEXSEQ_DTU}\` to \`true\` for the chosen analyses and \`false\` for all others."
need "At least one analysis must be chosen."
need "With BAM input no Salmon step runs"
# rMATS options
need "1. Annotated splice sites only (default) · 2. Also detect unannotated splice sites"
need "\`{RMATS_NOVEL}\` = \`false\` (answer 1, and whenever rMATS is not chosen) or \`true\` (answer 2)"
need "it is the threshold of rMATS's null hypothesis"
need "\`rmats_paired_stats\` = \`{PAIRED_DESIGN}\` (Step 5)"
need "SUPPA2's paired test \`diffsplice_paired\` = \`{PAIRED_DESIGN}\`"
[ "$(cfg_default diffsplice_paired)" = true ] && need "(its pipeline default \`true\` assumes paired samples)"
# DTU filter values: the pipeline defaults are quoted from the recorded config (schema and config agree at this revision).
need "The pipeline defaults (\`min_samps_gene_expr\` $(cfg_default min_samps_gene_expr), \`min_samps_feature_expr\` $(cfg_default min_samps_feature_expr), \`min_samps_feature_prop\` $(cfg_default min_samps_feature_prop))"
forbid "schema (6/0/0)"
need "\`{MIN_SAMPS_FEATURE}\` = the size of the smallest condition used in a contrast"
need "Tell the user the six filter values"
# Salmon route (pseudo_only at this gate): cost stated honestly (Task 4: the Salmon index build likely takes more than an hour)
[ "$(gv SALMON_ROUTE)" = pseudo_only ] && need "the STAR alignments feed rMATS, DEXSeq and edgeR; Salmon, run on the reads, feeds DTU and SUPPA2 once"
need "because \`pseudo_aligner\` has no off value at this revision"
forbid "it is quick"
# MISO dropped: why, in one sentence; never written keys of this revision
need "**MISO is not used**: in this pipeline it only draws sashimi plots for a short gene list"
[ "$(gv GATE_OUTCOME)" = B ] && need "\`isoformswitchanalyzer\` and \`leafcutter\` are parameters of this development revision and are written \`false\` (this skill does not offer them)"
need "\`rmats_variable_read_len\` and \`local_events\`"
need "writes no trimming key"
need "(edgeR function in this revision: $(gv EDGER_DEU_FUNCTION))"
# --- end Task 5 additions

# --- Task 5 review fixes
# I1: Step 8 decisions pinned with the left-hand side of each rule
need "if it ran (or will run) there for these samples, do not choose it here"
case "$(gv DTU_FILTER_SCOPE)" in
  all_samples) need "- \`{MIN_SAMPS_GENE_EXPR}\` = the number of samples in the samplesheet" ;;
  per_contrast) need "- \`{MIN_SAMPS_GENE_EXPR}\` = the smallest number of samples in a contrast" ;;
esac
need "**rMATS settings** (asked only when rMATS is chosen; otherwise the values below are written unchanged)"
need "with differential splicing between conditions. FASTQ input only."
need "DEXSeq, stageR), from Salmon. FASTQ input only."
# I2: yaml types of the literal lines of the module block, from the recorded schema (types at the parameter depth, 16 spaces;
# the group named isoformswitchanalyzer is itself an object one level up). string => "quoted"; integer, number, boolean => bare.
schema_types=$(awk '/^                "[A-Za-z_0-9]+": \{/ { k = $1; gsub(/[":{ ]/, "", k); next }
  k != "" && /^                    "type": "/ { t = $2; gsub(/[",]/, "", t); print k "\t" t; k = "" }' "$FIX/rnasplice_schema.json")
for kt in "aggregation	boolean" "n_dexseq_plot	integer" "min_feature_prop	number" "dtu_txi	string" "isoformswitchanalyzer	boolean"; do
  printf '%s\n' "$schema_types" | grep -qxF -- "$kt" || { echo "schema type extraction broken: $kt"; exit 2; }
done
# types_ok "<block text>": every key line with a literal value (no placeholder) has the yaml form of its schema type.
types_ok() {
  local line k v t
  while IFS= read -r line; do
    k=${line%%:*}; v=${line#*: }
    case "$v" in *"{"*) continue ;; esac
    t=$(printf '%s\n' "$schema_types" | awk -F'\t' -v k="$k" '$1 == k {print $2; exit}')
    case "$t" in
      string)  [[ $v =~ ^\"[^\"]*\"$ ]] ;;
      integer) [[ $v =~ ^-?[0-9]+$ ]] ;;
      number)  [[ $v =~ ^-?[0-9]+(\.[0-9]+)?([eE]-?[0-9]+)?$ ]] ;;
      boolean) [[ $v =~ ^(true|false)$ ]] ;;
      *) false ;;
    esac || { echo "FAIL: key $k = $v is not a yaml ${t:-?} (string: double-quoted; integer, number, boolean: bare)"; fail=1; }
  done < <(printf '%s\n' "$1" | grep -E '^[A-Za-z_0-9]+: ')
}
types_ok "$MOD"
# M1: the Akerberg et al. 2022 cut-offs are an example of one study, not a standard
need "As an example (not a standard), Akerberg et al. 2022 kept rMATS events with 0 uncalled replicates, FDR < 0.1 (zebrafish) or FDR < 0.05 (human), and |IncLevelDifference| > 0.1"
need "(from the parts of the paper's Methods that were accessible)"
forbid "usual cut-offs are"
# M2: when to choose DTU here
need "Choose DTU here only to get it in the same pipeline run from the FASTQ files"
need "for event-level splicing questions (which exons or events change) use rMATS"
need "differential expression plus DTU on the Salmon output of an existing nf-core/rnaseq run stays in \`/bulk-rnaseq-pipeline\`"
# --- end Task 5 review fixes
[ $fail -eq 0 ] && echo "PASS" || exit 1
