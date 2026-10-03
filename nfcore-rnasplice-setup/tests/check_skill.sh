#!/bin/bash
# Static checks for nfcore-rnasplice-setup.md. Usage: check_skill.sh <skill.md>
# Uses only the recorded fixtures of the gated pipeline revision (tests/fixtures/); no network, no python.
# The README checked is the README.md next to the skill file. CHECK_LIST_NEEDS=1 also prints "NEED <line> <text>" per need.
# CHECK_RUN_ALL=<file> checks that file instead of tests/run_all_tests.sh: a test hook of prove_mutations.sh (target runall) only.
# Exit codes: 0 PASS, 1 a check failed, 2 a fixture is missing or invalid (gate values are validated against their allowed sets).
# Call need/needr at top level (also inside a loop, `case` or `||`), never from inside a helper function: BASH_LINENO[0]
# would then be the line inside the helper, and prove_red.sh could not map the need to an added checker line.
set -u
SKILL=${1:?usage: check_skill.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd); FIX="$HERE/fixtures"; CUT="$HERE/cut_block.sh"
README="$(dirname "$SKILL")/README.md"
fail=0
[ -s "$SKILL" ] || { echo "FAIL: skill file missing or empty: $SKILL"; exit 1; }
for f in rnasplice_schema.json config_params.txt gate_values.tsv trace_process_names.txt trace_resources_g3.tsv output_tree.txt test_samplesheet.csv test_contrastsheet.csv schema_input.json schema_input_genome_bam.json modules.config realdata_star_align_args.txt realdata_resource_table.txt; do
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
# SALMON_INDEX_VERSION stays a validated fixture value, but the skill no longer uses it: by design
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
# STAR options of the STAR_ALIGN ext.args override (Step 10) and its explanation (copied from the pinned conf/modules.config)
allow="$allow quantMode quantTranscriptomeSAMoutput twopassMode outSAMtype readFilesCommand runRNGseed outFilterMultimapNmax alignSJDBoverhangMin outSAMattributes outSAMattrRGline outReadsUnmapped"
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
# a BAM sheet with strandedness/single_end columns carries the asked values in every row.
case "$(gv BAM_SHEET_HEADER)" in *,strandedness,single_end*)
  need "write \`{STRANDEDNESS}\` and \`true\` (single-end) or \`false\` in every row" ;; esac
# Task 2 review fixes (I1-I3, minors)
need "*reverse* = read 1 comes from the strand opposite to the transcript"
need "*forward* = read 1 comes from the transcript strand"
need "ask the user to check the kit; do not continue until they answer 1, 2 or 3"  # RSeQC bullet (the Salmon rules repeat the second half)
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
need "install_rnasplice_sheets \"\$T\" {SOURCE} 0 \"{SAMPLESHEET_CSV}\" \"{CONTRASTS_CSV}\""
need "install_rnasplice_sheets \"\$T\" {SOURCE} 1 \"{SAMPLESHEET_CSV}\" \"{CONTRASTS_CSV}\" \"{SHEET_PREFIX}_pairs.csv\""
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
# Design: one Bash call per procedure; STAR index reuse only for the gate's index format; sjdbOverhang note;
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
# no Salmon index reuse
need "\`{SALMON_INDEX}\` is always empty: the pipeline always builds its own Salmon index from the transcripts it extracts from the GTF"
need "lacks the GTF's non-coding transcripts, which would get no quantification (DTU, SUPPA2) without any error"
need "(BAM input uses no index: skip this part; \`{STAR_INDEX}\` is empty.)"
need "FASTA, GTF and the STAR index are reused"
need "not used: rnasplice always builds its own Salmon index"
need "SALMON_INDEX: 6 CPUs, 36 GB and 8 h"
forbid "likely takes more than an hour (not measured by this skill's verification)"
need "in the real-data test (human, Ensembl 116) it took 33 min with a peak of 19.8 GB (a mouse genome was not measured)"
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
forbid "did not exercise the reuse path (unverified)"
need "The real-data test exercised the reuse path with a human index from /nfcore-rnaseq-setup (\`versionGenome\` 2.7.4a, \`sjdbOverhang\` 93): the pipeline copied the 28.5 GB index into \`work/\` in about 1 minute (process STAR_GENOMEPARAMS_UPGRADE), and 88.5-92.7% of the reads of each sample mapped uniquely."
need "STAR_GENOMEGENERATE needs about 32 GB of memory and an hour or more (not measured by this skill: the real-data test reused an index)"
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

# --- Task 5 additions (fixtures of the pinned revision)
# every option key of every module group of the recorded schema is written exactly once (also for modules that are
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
# BAM input switches DTU and SUPPA2 off; no Salmon runs then (gate G5: zero Salmon tasks).
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

# --- Task 6 (Steps 10-11), gate branch B:
# no `nextflow pull` (a pinned revision downloads on first use), no salmon_index, Salmon selector '.*:SALMON_QUANT.*',
# no max_* params (HAS_MAX_PARAMS=no) but process.resourceLimits (HAS_RESOURCE_LIMITS=yes), no Nextflow environment helper.
need "## Step 10 — MultiQC title, output directory, nextflow.config"
need "## Step 11 — Params file, submission script and helper script"
for a in "**nextflow.config template.**" "**Params file template.**" "**Submission script.**" "**Genome download helper.**"; do anchor_once "$a"; done
need "1. \`results/{TODAY_ISO}_{WD_NAME}\` (default)"
need "If it exists, do not overwrite it"
need "\`multiqc_title\` is always double-quoted"
need "\`{PARAMS_YAML}\` = \`{SHEET_PREFIX}_params.yaml\`"
need "delete the \`star_index\` line when \`{STAR_INDEX}\` is empty"
need "sbatch --dependency=afterok:"
need "#SBATCH -t 48:00:00"
need "The head job only coordinates the pipeline (2 CPUs, 8 GB) but must outlive every task, hence 48 h."
need "the job reads the version from \`nextflow -version\` on the compute node"
need "The script fetches nothing itself: Nextflow downloads the pinned revision from GitHub the first time a job runs it"
need "Never run these scripts on the login node"
SUB=$(bash "$CUT" "$SKILL" "**Submission script.**" 2>/dev/null)
[ -n "$SUB" ] || { echo "FAIL: submission script block not found"; fail=1; }
launch=$(printf '%s\n' "$SUB" | awk '/^nextflow run nf-core\/rnasplice/ {p = 1} p {l = l $0; if ($0 !~ /\\[ \t]*$/) {print l; exit}}' | sed 's/\\[ \t]*/ /g')
[ "$launch" = "nextflow run nf-core/rnasplice -r {VERSION} -c nextflow.config -profile slurm,singularity -params-file {PARAMS_YAML}" ] \
  || { echo "FAIL: launch line must be exactly 'nextflow run nf-core/rnasplice -r {VERSION} -c nextflow.config -profile slurm,singularity -params-file {PARAMS_YAML}', found: '$launch'"; fail=1; }
[ "$(printf '%s\n' "$SUB" | grep -v '^[ \t]*#' | grep -v '^[ \t]*$' | tail -n 1)" = "$launch" ] \
  || { echo "FAIL: the launch line must be the last command of the submission script (its exit status is the job's)"; fail=1; }
all_nf=$(grep -cE '(^|[;&|(] *)(NXF_[A-Z_]+=[^ ]+ +)?nextflow +(run|pull|-version)' "$SKILL")
sub_nf=$(printf '%s\n' "$SUB" | grep -cE '(^|[;&|(] *)(NXF_[A-Z_]+=[^ ]+ +)?nextflow +(run|pull|-version)')
{ [ "$all_nf" -eq "$sub_nf" ] && [ "$sub_nf" -ge 2 ]; } \
  || { echo "FAIL: nextflow is invoked outside the submission script ($all_nf lines in the skill, $sub_nf in the script); the wizard must never run nextflow on the login node"; fail=1; }
while IFS= read -r l; do
  case "$l" in *"||"*) ;; *) echo "FAIL: unguarded module load (a failed module add silently breaks later ones): $l"; fail=1 ;; esac
done < <(grep -E '^[ \t]*module (add|load) ' "$SKILL")
# The module and conda lines that worked in every gate job, in that order, each guarded.
for l in 'module add miniconda3/v4 || { echo "ERROR: cannot load module miniconda3/v4" >&2; exit 1; }' \
         'source /home/software/conda/miniconda3/bin/condainit || { echo "ERROR: cannot source condainit" >&2; exit 1; }' \
         'conda activate {CONDA_ENV} || { echo "ERROR: cannot activate conda environment {CONDA_ENV}" >&2; exit 1; }' \
         'module add singularity/3.10.4 || { echo "ERROR: cannot load module singularity/3.10.4" >&2; exit 1; }' \
         'export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-$HOME/.singularity/cache}"'; do
  printf '%s\n' "$SUB" | grep -qxF -- "$l" || { echo "FAIL: the submission script must contain the line: $l"; fail=1; }
done
# no pull, and NXF_OFFLINE is neither set nor unset by the script (it may be set in ~/.bashrc; the first download works anyway).
! printf '%s\n' "$SUB" | grep -q 'NXF_OFFLINE' || { echo "FAIL: the submission script must not set or unset NXF_OFFLINE"; fail=1; }
printf '%s\n' "$SUB" | grep -qxF "NF_MIN=\"$(gv NEXTFLOW_MIN)\"; NF_MAX_EXCL=\"$(gv NEXTFLOW_MAX_EXCL)\"" || { echo "FAIL: the submission script must carry the verified Nextflow range"; fail=1; }
printf '%s\n' "$SUB" | grep -qxF '#SBATCH -n 2' && printf '%s\n' "$SUB" | grep -qxF '#SBATCH --mem=8G' && printf '%s\n' "$SUB" | grep -qxF '#SBATCH -p bcc' \
  || { echo "FAIL: the head job must request -n 2, --mem=8G and -p bcc (the head job only coordinates the pipeline)"; fail=1; }
# nextflow.config: the selector table checked against every process name of the verification runs and the
# resources the runs requested and used. A row gives '.*:ROW' when every real process name ending in a component that starts
# with ROW ends in exactly ROW; '.*:ROW[XY]' when every such component is ROW plus one letter (SUPPA2: DIFFSPLICE_IOE,
# DIFFSPLICE_IOI); and '.*:ROW.*' otherwise (Salmon: SALMON_QUANT_SALMON, SALMON_QUANT_STAR); rows that match no process are
# left out. Nextflow matches withName regexes against the whole process name (grep -Ex here).
# The fifth column of SEL_ROWS is the evidence of the row's values:
#   gate      unchanged since the verification gate: every process it matches requested exactly these values in the reference
#             run G3 (trace_resources_g3.tsv) and in the real-data test (realdata_resource_table.txt), so it was applied as written;
#   realdata  set from the real-data test (2026-10-02, one human dataset, 6 samples; not yet applied in a run): every process it
#             matches ran in that test, and its observed peak memory and longest run time are below the requested values.
# Both kinds must cover the real-data observations (peak_rss < memory, realtime < time); the one exception is the DEXSEQ_EXON
# peak_rss of 187 GB, which counts forked BiocParallel workers separately (the task was not killed at 32 GB; sstat 8.4 GB).
CFG=$(bash "$CUT" "$SKILL" "**nextflow.config template.**" 2>/dev/null)
[ -n "$CFG" ] || { echo "FAIL: nextflow.config template not found"; fail=1; }
# The resource selectors live in the slurm profile (CFGS); the only other withName block is the STAR_ALIGN ext.args override
# of the top-level process block (checked in the real-data fix round section below).
CFGS=$(printf '%s\n' "$CFG" | awk '/^    slurm \{$/ {s = 1} s {print} s && /^    \}$/ {exit}')
SEL_ROWS="STAR_GENOMEGENERATE	8	64 GB	8h	gate
STAR_ALIGN	8	48 GB	4h	realdata
RMATS_PREP	2	8 GB	4h	realdata
RMATS_POST	8	16 GB	8h	realdata
DEXSEQ_COUNT	1	4 GB	8h	realdata
DEXSEQ_EXON	8	32 GB	8h	gate
DEXSEQ_DTU	8	16 GB	4h	realdata
SALMON_QUANT	8	32 GB	4h	realdata
DIFFSPLICE_IO	1	8 GB	16h	realdata
MAKE_TRANSCRIPTS_FASTA	2	8 GB	2h	realdata"
RD_OBS="$FIX/realdata_resource_table.txt"
# realdata_resource_table.txt (copy of the real-data test's evidence/resource_table.txt): process n cpus memGB timeh peakGB mem%
# maxmin time% meancpu maxcpu att, one row per short process name, plus its header row.
[ "$(awk 'NF == 12 && $1 != "process" && $2 ~ /^[0-9]+$/ && $6 ~ /^[0-9.]+$/ && $8 ~ /^[0-9.]+$/' "$RD_OBS" | wc -l)" -ge 35 ] \
  || { echo "fixture realdata_resource_table.txt broken (process n cpus memGB timeh peakGB mem% maxmin time% meancpu maxcpu att)"; exit 2; }
NAMES="$FIX/trace_process_names.txt"; OBS="$FIX/trace_resources_g3.tsv"
[ "$(wc -l < "$NAMES")" -eq 75 ] || { echo "fixture trace_process_names.txt must hold the 75 process names of the verification runs"; exit 2; }
# selector<TAB>cpus<TAB>memory<TAB>time for each withName block of the template
sels=$(printf '%s\n' "$CFGS" | awk -v q="'" '
  $0 ~ "^[ \t]*withName: " q "[^" q "]+" q " [{]$" { s = $0; sub("^[ \t]*withName: " q, "", s); sub(q " [{]$", "", s); c = m = t = ""; next }
  s != "" && /^[ \t]*cpus = [0-9]+$/ { c = $3 }
  s != "" && $0 ~ "^[ \t]*memory = " q { m = $0; sub("^[ \t]*memory = " q, "", m); sub(q "$", "", m) }
  s != "" && $0 ~ "^[ \t]*time = " q { t = $0; sub("^[ \t]*time = " q, "", t); sub(q "$", "", t) }
  s != "" && /^[ \t]*[}]$/ { print s "\t" c "\t" m "\t" t; s = "" }')
[ -n "$sels" ] || { echo "FAIL: no withName selectors in the nextflow.config template"; fail=1; }
while IFS=$'\t' read -r re c m t; do
  [ -n "$re" ] || continue
  row=$(printf '%s\n' "$re" | sed -nE 's/^\.\*:([A-Z0-9_]+)(\.\*|\[[A-Z]+\])?$/\1/p')
  rres=$(printf '%s\n' "$SEL_ROWS" | awk -F'\t' -v r="$row" '$1 == r {print $2 "\t" $3 "\t" $4}')
  [ -n "$row" ] && [ -n "$rres" ] || { echo "FAIL: selector $re is not a row of the selector table ('.*:ROW', '.*:ROW[XY]' or '.*:ROW.*')"; fail=1; continue; }
  [ "$c	$m	$t" = "$rres" ] || { echo "FAIL: selector $re requests $c CPUs, $m, $t; the selector table says $(echo "$rres" | tr '\t' ' ')"; fail=1; }
  hits=$(grep -Ex -- "$re" "$NAMES")
  [ -n "$hits" ] || { echo "FAIL: selector $re matches no process of the verification runs ($(wc -l < "$NAMES") names)"; fail=1; continue; }
  while IFS= read -r h; do
    case "${h##*:}" in "$row"*) ;; *) echo "FAIL: selector $re also matches $h, which is not a $row process"; fail=1 ;; esac
  done <<< "$hits"
done <<< "$sels"
while IFS=$'\t' read -r row c m t ev; do
  comps=$(sed 's/.*://' "$NAMES" | grep -E "^$row" | sort -u)
  [ -n "$comps" ] || continue
  if [ "$(printf '%s\n' "$comps" | grep -vx -- "$row")" = "" ]; then want=".*:$row"
  elif [ "$(printf '%s\n' "$comps" | grep -vxE -- "$row[A-Z]")" = "" ]; then want=".*:$row[$(printf '%s\n' "$comps" | sed "s/^$row//" | sort | tr -d '\n')]"
  else want=".*:$row.*"; fi
  printf '%s\n' "$sels" | cut -f1 | grep -qxF -- "$want" \
    || { echo "FAIL: process $(echo $comps | tr ' ' ',') ran in the verification runs but has no withName selector '$want'"; fail=1; }
  while IFS= read -r n; do
    k=0; while IFS= read -r re; do [ -n "$re" ] && printf '%s\n' "$n" | grep -Eqx -- "$re" && k=$((k + 1)); done < <(printf '%s\n' "$sels" | cut -f1)
    [ "$k" -eq 1 ] || { echo "FAIL: process $n is matched by $k withName selectors (expected exactly 1)"; fail=1; }
  done < <(awk -F: -v r="$row" 'index($NF, r) == 1' "$NAMES")
done <<< "$SEL_ROWS"
# Requested resources of the reference run G3 and of the real-data test: every 'gate' row must have been applied as written in
# both (its processes requested exactly the template's values), and no process of G3 ran on the bare process default of the
# template (that would be a process with no label and no selector). 'realdata' rows are compared with the observations below.
sel_ev() { printf '%s\n' "$SEL_ROWS" | awk -F'\t' -v r="$(printf '%s\n' "$1" | sed -nE 's/^\.\*:([A-Z0-9_]+)(\.\*|\[[A-Z]+\])?$/\1/p')" '$1 == r {print $5}'; }
dflt=$(printf '%s\n' "$CFGS" | awk -v q="'" '/withName/ {exit} /^[ \t]*cpus = [0-9]+$/ {c = $3} $0 ~ "^[ \t]*memory = " q {m = $3 " " $4} $0 ~ "^[ \t]*time = " q {t = $3} END {gsub(q, "", m); gsub(q, "", t); print c "\t" m "\t" t}')
[ "$dflt" = "2	8 GB	4h" ] || { echo "FAIL: the process default of the slurm profile must be 2 CPUs, 8 GB, 4h (found: $(echo "$dflt" | tr '\t' ' '))"; fail=1; }
[ "$(awk -F'\t' 'NF == 4' "$OBS" | wc -l)" -ge 40 ] || { echo "fixture trace_resources_g3.tsv broken (process<TAB>cpus<TAB>memory<TAB>time)"; exit 2; }
# selector<TAB>cpus<TAB>memory<TAB>time<TAB>evidence (gate or realdata) for each withName block of the slurm profile
sels_ev=$(while IFS=$'\t' read -r re c m t; do [ -n "$re" ] && printf '%s\t%s\t%s\t%s\t%s\n' "$re" "$c" "$m" "$t" "$(sel_ev "$re")"; done <<< "$sels")
while IFS=$'\t' read -r n oc om ot; do
  while IFS=$'\t' read -r re c m t ev; do
    [ "$ev" = gate ] || continue
    printf '%s\n' "$n" | grep -Eqx -- "$re" || continue
    [ "$oc	$om	$ot" = "$c	$m	$t" ] || { echo "FAIL: selector $re: the verification run requested $oc CPUs, $om, $ot for $n, the template says $c, $m, $t"; fail=1; }
  done <<< "$sels_ev"
  [ "$oc	$om	$ot" != "$dflt" ] || { echo "FAIL: process $n ran on the bare process default ($oc CPUs, $om, $ot): no label and no selector"; fail=1; }
done < "$OBS"
# The real-data table names processes by their last component; a selector matches it as 'X:<name>' (one awk pass).
rd_fail=$(printf '%s\n' "$sels_ev" | awk -F'\t' '
  FNR == NR { if ($1 == "") next; ns++; re[ns] = $1; c[ns] = $2; m[ns] = $3; t[ns] = $4; ev[ns] = $5; next }
  { split($0, a, " "); if (a[1] == "process" || a[1] == "") next
    for (s = 1; s <= ns; s++) {
      if (("X:" a[1]) !~ ("^(" re[s] ")$")) continue
      nrd++; seen[s] = 1
      if (ev[s] == "gate" && (a[3] "\t" a[4] " GB\t" a[5] "h") != (c[s] "\t" m[s] "\t" t[s]))
        print "FAIL: selector " re[s] ": the real-data test requested " a[3] " CPUs, " a[4] " GB, " a[5] "h for " a[1] ", the template says " c[s] ", " m[s] ", " t[s]
      mg = m[s]; sub(/ GB$/, "", mg); th = t[s]; sub(/h$/, "", th)
      if (a[1] == "DEXSEQ_EXON") { if (a[6] != "187.44") print "FAIL: the DEXSEQ_EXON peak_rss exception covers the recorded 187.44 GB only (found " a[6] ")" }
      else if (!(a[6] + 0 < mg + 0)) print "FAIL: selector " re[s] " requests " m[s] ", but " a[1] " used up to " a[6] " GB in the real-data test"
      if (!(a[8] + 0 < th * 60)) print "FAIL: selector " re[s] " requests " t[s] ", but " a[1] " ran up to " a[8] " min in the real-data test"
    } }
  END {
    for (s = 1; s <= ns; s++) if (ev[s] == "realdata" && !(s in seen)) print "FAIL: selector " re[s] " is a realdata row, but none of its processes ran in the real-data test"
    if (nrd != 10) print "FAIL: " nrd + 0 " real-data processes are matched by a selector (expected 10)"
  }' - "$RD_OBS")
[ -z "$rd_fail" ] || { printf '%s\n' "$rd_fail"; fail=1; }
[ "$(printf '%s\n' "$CFG" | grep -cE '^(timeline|report|trace|dag) +\{ enabled = true; overwrite = true;')" -eq 4 ] || { echo "FAIL: overwrite = true must be set for timeline, report, trace and dag"; fail=1; }
printf '%s\n' "$CFG" | grep -qF "fields = 'task_id,hash,native_id,name,status,exit,attempt,cpus,memory,time,realtime,%cpu,peak_rss'" || { echo "FAIL: trace fields must include the requested cpus, memory and time, the attempt and %cpu (D4 of the real-data test)"; fail=1; }
if [ "$(gv HAS_RESOURCE_LIMITS)" = yes ]; then
  printf '%s\n' "$CFG" | grep -qxF "    resourceLimits = [ cpus: 16, memory: '64 GB', time: '24h' ]" || { echo "FAIL: the resourceLimits line is missing"; fail=1; }
else
  ! printf '%s\n' "$CFG" | grep -q resourceLimits || { echo "FAIL: resourceLimits is not known to Nextflow $(gv NEXTFLOW_TESTED)"; fail=1; }
fi
[ "$(gv HAS_RESOURCE_LIMITS)" = yes ] && need "\`resourceLimits\` caps every task at 16 CPUs, 64 GB and 24 h"
need "\`SALMON_QUANT_SALMON\` and \`SALMON_QUANT_STAR\`"
need "no SALMON_QUANT_STAR task runs; the pattern covers both names."
forbid "selector was only matched against the process names recorded there"
forbid "is first applied in this skill's cluster acceptance run"
forbid "checked on the nf-core test data of this skill's verification only"
# Params file template
BASEY=$(bash "$CUT" "$SKILL" "**Params file template.**" 2>/dev/null)
[ "$(printf '%s\n' "$BASEY" | grep -cx '{MODULE_PARAMS}')" -eq 1 ] || { echo "FAIL: the params template must contain one {MODULE_PARAMS} line"; fail=1; }
for kv in 'input: "{SAMPLESHEET_CSV}"' 'contrasts: "{CONTRASTS_CSV}"' 'source: "{SOURCE}"' 'outdir: "{OUTDIR}"' 'multiqc_title: "{MULTIQC_TITLE}"' \
          'fasta: "{FASTA_PATH}"' 'gtf: "{GTF_PATH}"' 'star_index: "{STAR_INDEX}"' 'gencode: false' 'save_reference: false'; do
  printf '%s\n' "$BASEY" | grep -qxF -- "$kv" || { echo "FAIL: the params template must contain the line: $kv"; fail=1; }
done
! printf '%s\n' "$BASEY" | grep -q '^salmon_index:' || { echo "FAIL: the params template must not write salmon_index (the pipeline always builds its own Salmon index)"; fail=1; }
if [ "$(gv HAS_MAX_PARAMS)" = no ]; then
  ! printf '%s\n' "$BASEY" | grep -qE '^max_(cpus|memory|time):' || { echo "FAIL: max_cpus/max_memory/max_time are not parameters of this revision (Nextflow stops at launch)"; fail=1; }
fi
dups=$(printf '%s\n%s\n' "$BASEY" "$MOD" | grep -oE '^[A-Za-z_0-9]+:' | sort | uniq -d | tr '\n' ' ')
[ -z "$dups" ] || { echo "FAIL: key written in both the params template and the module block: $dups"; fail=1; }
defaults_equal "$BASEY"
types_ok "$BASEY"
# Genome download helper: compute node, resumable wget, gunzip -c, set -u, absolute directory, no rm.
HELP=$(bash "$CUT" "$SKILL" "**Genome download helper.**" 2>/dev/null)
[ -n "$HELP" ] || { echo "FAIL: genome download helper not found"; fail=1; }
printf '%s\n' "$HELP" | grep -qx 'set -u' || { echo "FAIL: the download helper must run with set -u"; fail=1; }
printf '%s\n' "$HELP" | grep -q '^#SBATCH -p bcc$' || { echo "FAIL: the download helper must be an sbatch script for partition bcc"; fail=1; }
! printf '%s\n' "$HELP" | grep -qE '(^|[^A-Za-z_])rm( |$)' || { echo "FAIL: the download helper must not remove files (no rm with variable paths)"; fail=1; }
! printf '%s\n' "$HELP" | grep -q 'gunzip -k' || { echo "FAIL: gunzip -k is not available on CentOS 7; use gunzip -c file.gz > file"; fail=1; }
printf '%s\n' "$HELP" | grep -q 'wget -c ' && printf '%s\n' "$HELP" | grep -q 'gunzip -c ' || { echo "FAIL: the download helper must use wget -c and gunzip -c"; fail=1; }
need "write \`download_genome_{REF_TAG}.sh\`"
# Step 7 promises the helper (Step 11 now also names the file, so the Task 4 need alone no longer pins the Step 7 sentence).
need "and generate \`download_genome_{REF_TAG}.sh\` (Step 11)"
forbid "create_nextflow_env_"
forbid "build_star_index"
# --- end Task 6

# --- Task 6 review fixes
# I1: the job changes to {CWD} first (every path it uses is relative to {CWD}), then checks its files before anything else.
for l in 'cd "{CWD}" || { echo "ERROR: cannot change to {CWD}" >&2; exit 1; }' \
         'for f in nextflow.config {PARAMS_YAML} {SAMPLESHEET_CSV} {CONTRASTS_CSV}; do [ -s "$f" ] || { echo "ERROR: $f is missing or empty in {CWD}" >&2; exit 1; }; done'; do
  printf '%s\n' "$SUB" | grep -qxF -- "$l" || { echo "FAIL: the submission script must contain the line: $l"; fail=1; }
done
[ "$(printf '%s\n' "$SUB" | grep -v '^#' | grep -v '^[ \t]*$' | head -n 1)" = 'cd "{CWD}" || { echo "ERROR: cannot change to {CWD}" >&2; exit 1; }' ] \
  || { echo "FAIL: cd \"{CWD}\" must be the first command of the submission script"; fail=1; }
need "it changes to \`{CWD}\` first, so it works wherever it is submitted from"
# I2: executor, queue, Singularity, mail and node lines of the config and both scripts; --parsable in the dependency line.
for l in "            executor = 'slurm'" "            queue = 'bcc'" "            enabled = true" "            autoMounts = true"; do
  printf '%s\n' "$CFG" | grep -qxF -- "$l" || { echo "FAIL: the nextflow.config template must contain the line: $l"; fail=1; }
done
printf '%s\n' "$CFG" | awk '/^    singularity \{$/ {s = 1} s && /^            enabled = true$/ {e = 1} s && /^            autoMounts = true$/ {a = 1} s && /^    \}$/ {exit} END {exit !(e && a)}' \
  || { echo "FAIL: the singularity profile must set enabled = true and autoMounts = true"; fail=1; }
printf '%s\n' "$CFG" | awk '/^    slurm \{$/ {s = 1} s && /^            executor = .slurm.$/ {e = 1} s && /^            queue = .bcc.$/ {q = 1} s && /^    \}$/ {exit} END {exit !(e && q)}' \
  || { echo "FAIL: the slurm profile must set executor = 'slurm' and queue = 'bcc'"; fail=1; }
for l in '#SBATCH -N 1' '#SBATCH --mail-type=END,FAIL' '#SBATCH --mail-user={USER_EMAIL}'; do
  printf '%s\n' "$SUB" | grep -qxF -- "$l" || { echo "FAIL: the submission script must contain the line: $l"; fail=1; }
  printf '%s\n' "$HELP" | grep -qxF -- "$l" || { echo "FAIL: the download helper must contain the line: $l"; fail=1; }
done
for l in '#SBATCH -n 2' '#SBATCH --mem=8G' '#SBATCH -t 4:00:00'; do
  printf '%s\n' "$HELP" | grep -qxF -- "$l" || { echo "FAIL: the download helper must contain the line: $l"; fail=1; }
done
need "jid=\$(sbatch --parsable download_genome_{REF_TAG}.sh) && sbatch --dependency=afterok:\$jid nf-core_rnasplice_{VERSION_TAG}.sh"
# M2: wizard text of Steps 10-11
need "**MultiQC title** (numbered): 1. \`{SEQ_DATE}_{WD_NAME}\` · 2. \`{TODAY_YYMMDD}_{WD_NAME}\` · 3. Custom. Store \`{MULTIQC_TITLE}\`."
need "· 2. \`results/{TODAY_ISO}_{SHEET_PREFIX}\` · 3. Custom. Store \`{OUTDIR}\`."
need "if it exists and is not empty, ask (numbered): 1. use it anyway (the pipeline adds to it and may overwrite files) · 2. choose another."
need "Check for an existing config: \`ls \"{CWD}/nextflow.config\"\`."
need "To submit:  sbatch nf-core_rnasplice_{VERSION_TAG}.sh"
need "Show every written file in full."
# M3: custom title and directory characters (a double quote or backslash would break the double-quoted YAML value)
need "A custom title or output directory may contain only letters, digits, \`.\`, \`_\`, \`-\` and (in the directory) \`/\`"
need "contains a character the params file cannot hold"
# M5, M6: fixed resources do not grow on retry; -resume only in the documented procedure, never in the generated script
need "Fixed \`withName\` resources do not grow when the pipeline retries a task"
need "raise that selector's \`memory\` or \`time\` in nextflow.config (up to the \`resourceLimits\` caps) and resubmit as described under **Resuming a run**"
anchor_once "**Resuming a run.**"
need "\`-resume\` is a Nextflow option, not a pipeline parameter; the wizard never writes it into the generated script."
resume_out=$(grep -n -- '-resume' "$SKILL" | grep -v '^[0-9]*:\*\*Resuming a run\.\*\*' | cut -d: -f1 | tr '\n' ' ')
[ -z "$resume_out" ] || { echo "FAIL: -resume appears outside the **Resuming a run.** paragraph (line $resume_out)"; fail=1; }
! printf '%s\n' "$SUB" | grep -q -- '-resume' || { echo "FAIL: the generated submission script must not carry -resume"; fail=1; }
# M8: the helper reports the move of an invalid download only when the move succeeded
need "could not be moved aside: rename it by hand, then run the helper again"
# --- end Task 6 review fixes

# --- Task 7 (Step 12, Notes)
need "## Step 12 — Summary and hand-off"
anchor_once "**Hand-off note.**"
HO=$(bash "$CUT" "$SKILL" "**Hand-off note.**" 2>/dev/null)
n=0
while IFS= read -r p; do
  n=$((n + 1)); q=${p//\{CONTRAST\}/$(gv TEST_CONTRAST)}; q=${q%/}
  grep -qxF -- "./$q" "$FIX/output_tree.txt" || { echo "FAIL: hand-off path $p not in the recorded output tree"; fail=1; }
done < <(printf '%s\n' "$HO" | sed -n 's/^  \([^ ]\{1,\}\)  .*$/\1/p')
[ "$n" -ge 6 ] || { echo "FAIL: the hand-off note lists $n output paths (at least 6 expected)"; fail=1; }
! printf '%s\n' "$HO" | grep -q '<' || { echo "FAIL: the hand-off note still contains a <...> placeholder"; fail=1; }
need ".MATS.JC.txt"
need ".MATS.JCEC.txt"
case "$(gv RMATS_B1_GROUP)" in
  treatment) need "the treatment samples are rMATS's \`b1\`" ;;
  control) need "this pipeline passes the control samples as rMATS's \`b1\`" ;;
esac
need "**rnasplice or \`/bulk-rnaseq-pipeline\`?**"
need "do not run both on the same samples"
need "This skill does not analyse the results"
need "## Notes for the assistant"
need "Raw FASTQ/BAM files are read-only"
need "Never pass pipeline parameters on the \`nextflow run\` line"
need "Never run Nextflow, conda, Java, Python or R on the login node"
[ "$(gv SALMON_ROUTE)" = star_salmon_both ] && need "use the \`star_salmon\` ones"
# cut-offs are the user's choice (Step 8 gives Akerberg et al. 2022 as an example only);
# DTU/SUPPA2 naming and the reversed DTU fold-change sign (gate_report.md, output tree); BAM input lists only the three STAR modules;
# the wizard submits nothing itself and says how to follow the run; the pinned unreleased revision; test data only.
forbid "Usual filters"
forbid "validated on real data"
need "This skill sets no cut-offs: choose your own"
need "\`IncLevelDifference\` = mean(\`IncLevel1\`) - mean(\`IncLevel2\`)"
need "DTU and SUPPA2 files are named \`{TREATMENT}-{CONTROL}\` (for example \`KO-WT\`), not after the contrast name"
need "the fold-change column is \`log2fold_{CONTROL}_{TREATMENT}\` (control over treatment)"
# (Task 7 review C1: the need on "The sign is reversed relative to rMATS and SUPPA2" was wrong and is now forbidden below.)
need "With BAM input the note has only the rMATS, DEXSeq exon usage and edgeR exon usage lines"
need "With a paired design (\`rmats_paired_stats: true\`) the rMATS directory is \`star/rmats/{CONTRAST}_paired/\`"
need "The wizard does not start the pipeline itself"
need "squeue -u \$USER"
need "tail -f nf-core_rnasplice_{VERSION_TAG}.{JOBID}.log"
forbid "verified on nf-core's test data only"
need "it was verified on nf-core's test data and run once on one real human dataset (6 samples, see the README); do not present its settings or outputs as validated beyond these two datasets."
if [ "$(gv GATE_OUTCOME)" = B ]; then
  need "is an unreleased development commit of nf-core/rnasplice (\`$(gv VERSION_TAG)\`), pinned by this skill"
  need "needs this skill's verification gate to be run again"
fi
# The hand-off note lists one line per analysis, in this order, plus MultiQC (directory names from output_tree.txt).
ho_paths=$(printf '%s\n' "$HO" | sed -n 's/^  \([^ ]\{1,\}\)  .*$/\1/p' | tr '\n' ' ')
[ "$ho_paths" = "star/rmats/{CONTRAST}/ star/rmats/{CONTRAST}/SE.MATS.JC.txt star/dexseq_exon/results/ star/edger/ salmon/dexseq_dtu/results/ salmon/suppa/ multiqc/ " ] \
  || { echo "FAIL: the hand-off note must list rMATS, SE.MATS.JC.txt, DEXSeq exon, edgeR exon, DTU, SUPPA2 and MultiQC in this order (found: $ho_paths)"; fail=1; }
# --- end Task 7

# --- Task 7 review fixes
# C1: SUPPA2 dPSI = mean(control) - mean(treatment) although the header reads treatment-control (gate_report.md, corrected after
# the Task 7 review: 1098/1098 local and 1854/1854 isoform events of G3). Same direction as DTU, opposite to rMATS.
forbid "The sign is reversed relative to rMATS and SUPPA2"
forbid "SUPPA2 dPSI is treatment minus control"
need "SUPPA2 dPSI = mean PSI of the control minus mean PSI of the treatment, although the column header reads \`local_{TREATMENT}-local_{CONTROL}_dPSI\`"
need "a positive dPSI means more inclusion in the control"
need "⚠️ The sign is reversed relative to rMATS: DTU fold changes, SUPPA2 dPSI and DEXSeq exon fold changes are positive when the control has more; rMATS \`IncLevelDifference\` and edgeR \`logFC\` are positive when the treatment has more."
# I1: every plain-language sign clause is pinned (each flip fails).
case "$(gv RMATS_B1_GROUP)" in
  treatment)
    need "= inclusion level of the treatment group minus that of the control group (in this pipeline"
    need "positive values mean more inclusion in the treatment."
    need "hold one comma-separated value per treatment sample, in samplesheet order (the \`_2\` columns: the control samples)" ;;
  control)
    need "= inclusion level of the control group minus that of the treatment group (this pipeline"
    need "positive values mean more inclusion in the control, so flip the sign to read it as treatment minus control" ;;
esac
need "so a positive value means a larger share of the transcript in the control"
# I2: exon-usage directions (gate_report.md: DEXSeq exon from header + run_dexseq_exon.R; edgeR from run_edger_exon.R, code only).
forbid "read it from the column names"
need "the fold-change column is \`log2fold_{CONTROL}_{TREATMENT}\` (control relative to treatment, like DTU"
need "so a positive value means more usage of the exon bin in the control"
need "\`logFC\` = treatment minus control (exon usage relative to its gene; from the pinned \`run_edger_exon.R\`"
need "so a positive value means more usage of the exon in the treatment"
# M1: placeholders of the hand-off note and several contrasts.
need "\`{CONTRAST}\` = its \`contrast\` column, \`{TREATMENT}\` and \`{CONTROL}\` = its \`treatment\` and \`control\` columns"
need "print every line that contains one of them once per contrast"
# M2: standing rules in the Notes.
need "- Shell variables and functions do not persist between Bash tool calls: a procedure defines and calls its functions"
need "Never pre-fetch the pipeline (no Nextflow \`pull\` command"
# M3, M4, M8, M9: weak sentences pinned.
need "The output paths and sign conventions in this step were verified on nf-core's test data (this skill's verification run of the pinned revision) and on one real dataset, the real-data test (human, 6 samples): there every path of the hand-off note existed and every sign statement held on real rows (the edgeR direction on one exon); no other real dataset was run."
need "first confirm the directory names against the pipeline's \`docs/output.md\` for \`{VERSION}\` with \`WebFetch\`"
need "(and MultiQC): DTU and SUPPA2 are off."
need "a \`.MATS.JC.txt\` table (junction-spanning reads only) and a \`.MATS.JCEC.txt\` table (junction and exon-body reads)"
# M5: the tail line follows the #SBATCH -o line of the submission script; the log lands where sbatch ran.
slog=$(printf '%s\n' "$SUB" | sed -n 's/^#SBATCH -o \(.*\)$/\1/p')
jobid_ph='{JOBID}'
if [ -n "$slog" ]; then
  need "tail -f ${slog//%j/$jobid_ph}"
else
  echo "FAIL: the submission script has no #SBATCH -o line"; fail=1
fi
need "SLURM writes this log in the directory where \`sbatch\` was run"
# M6: no Nextflow command for the login node in the wizard text.
forbid "nextflow log"
# M7: the bulk DTU module runs the same steps on other inputs.
need "import and Salmon index differ"
# --- end Task 7 review fixes

# --- Task 8 (README, run_all_tests.sh), gate branch B
# (no explicit pull, STAR-only index reuse, both paired-test defaults, corrected directions, resources, BAM verified from code,
# unreleased pin, test runner with --full, validation on nf-core test data only, honest limitations).
forbidr() { ! grep -qF -- "$1" "$README" 2>/dev/null || { echo "FAIL: README has forbidden text: $1"; fail=1; }; }
needr "# \`/nfcore-rnasplice-setup\` — nf-core/rnasplice Differential Splicing Setup Skill"
needr "cp nfcore-rnasplice-setup.md ~/.claude/commands/"
needr "## Validation status"
needr "nf-core/rnasplice $(gv PIPELINE_REVISION)"
needr "Nextflow $(gv NEXTFLOW_TESTED)"
needr "4 paired-end human chrX samples"
forbidr "No real biological data"
needr "This is one dataset, not a validation on real data in general."
needr "## Known limitations"
needr "## rnasplice or \`/bulk-rnaseq-pipeline\`"
needr "sashimi_plot"
needr "tests/run_all_tests.sh"
needr "## Where the results are"
# no plan variant marker may survive in the skill or the README
VM='\[(A/C|A/B|B|C|BAM|NO-BAM|fixed|from_sheet|none|pseudo_only|star_salmon_only|star_salmon_both|treatment|control|sheet|sorted_by_name|unordered|all_samples|per_contrast)\]|\[C: '
! grep -nE "$VM" "$SKILL" || { echo "FAIL: a plan variant marker is left in the skill (lines above)"; fail=1; }
! grep -nE "$VM" "$README" 2>/dev/null || { echo "FAIL: a plan variant marker is left in the README (lines above)"; fail=1; }
# Task 9: acceptance DONE, as a whole line, and the PENDING line gone (as in ase-pipeline's checker).
grep -qE '^Cluster acceptance \(nf-core test data\): DONE \(2026-[0-9-]+\)$' "$README" 2>/dev/null && ! grep -qF 'Cluster acceptance (nf-core test data): PENDING' "$README" 2>/dev/null \
  || { echo "FAIL: README must have the whole line 'Cluster acceptance (nf-core test data): DONE (2026-MM-DD)' and no PENDING line"; fail=1; }
# Correction 9: every path line of the Step 12 hand-off note is a row of the README table, with the same description.
while IFS=$'\t' read -r p d; do
  needr "| \`$p\` | $d |"
done < <(printf '%s\n' "$HO" | sed -n 's/^  \([^ ]\{1,\}\)  \(.*\)$/\1\t\2/p')
# Correction 1: no explicit pull and no NXF_OFFLINE override; the compute nodes download the revision and the containers.
forbidr "nextflow pull"
forbidr "NXF_OFFLINE=false"
needr "| Internet on compute nodes | The job downloads the pinned revision of the pipeline from GitHub and the containers"
needr "| Nextflow $(gv NEXTFLOW_MIN) or newer | In a conda environment; verified with Nextflow $(gv NEXTFLOW_TESTED) (env \`$(gv CONDA_ENV_TESTED)\`)"
# Correction 2: only the STAR index is reused; the Salmon index is always built by the pipeline.
needr "An existing STAR index is reused only when compatible (\`versionGenome $(gv STAR_VERSION_GENOME)\`); the Salmon index is always built by the pipeline, never reused |"
forbidr "a path this skill's verification did not exercise: unverified"
needr "at this revision the pipeline copies a given index into its \`work/\` directory: in the real-data test, 28.5 GB in about 1 minute, with 88.5-92.7% of the reads of each sample mapped uniquely)"
needr "lacks the non-coding transcripts"
forbidr "Salmon indexes reused"
forbidr "compatible STAR/Salmon indexes reused"
# Correction 3: both paired-test defaults (recorded pipeline configuration) and the skill's rule.
needr "At this revision \`rmats_paired_stats\` defaults to \`$(cfg_default rmats_paired_stats)\` and SUPPA2's \`diffsplice_paired\` to \`$(cfg_default diffsplice_paired)\`. The skill writes both from the confirmed pairing"
needr "which is asked only with exactly two conditions of equal size"
needr "so that the i-th sample of one condition is paired with the i-th sample of the other"
# Correction 4: directions, as Step 12 states them (gate_report.md, SUPPA2 corrected after the Task 7 review).
forbidr "SUPPA2 matches rMATS"
case "$(gv RMATS_B1_GROUP)" in
  treatment) needr "| rMATS | \`IncLevelDifference\` = mean(\`IncLevel1\`) − mean(\`IncLevel2\`); \`b1\` = the treatment samples | more inclusion in the treatment |" ;;
  control) needr "| rMATS | \`IncLevelDifference\` = mean(\`IncLevel1\`) − mean(\`IncLevel2\`); \`b1\` = the control samples | more inclusion in the control |" ;;
esac
needr "| SUPPA2 | dPSI = mean PSI of the control − mean PSI of the treatment, although the header reads \`local_{TREATMENT}-local_{CONTROL}_dPSI\` | more inclusion in the control |"
needr "| DEXSeq DTU | \`log2fold_{CONTROL}_{TREATMENT}\` | a larger share of the transcript in the control |"
needr "| DEXSeq exon usage | \`log2fold_{CONTROL}_{TREATMENT}\` | more usage of the exon bin in the control |"
needr "| edgeR exon usage | \`logFC\` = treatment − control | more usage of the exon in the treatment | pinned code; checked on one exon of the real-data test (logFC −1.59 on the RBPMS2 exon skipped in the KO) |"
# Correction 5: resources; every selector row of the README table equals the Step 10 template; the head job equals Step 11.
while IFS= read -r row; do
  needr "$row"
done < <(printf '%s\n' "$CFGS" | awk -v q="'" '
  /withName:/ { n = $0; sub(/^[^:]*:[ \t]*/, "", n); gsub(q, "", n); sub(/^\.\*:/, "", n); sub(/[ \t]*\{.*$/, "", n) }
  n != "" && /^[ \t]*cpus[ \t]*=/ { c = $NF }
  n != "" && /^[ \t]*memory[ \t]*=/ { m = $0; sub(/^[^=]*=[ \t]*/, "", m); gsub(q, "", m) }
  n != "" && /^[ \t]*time[ \t]*=/ { t = $NF; gsub(q, "", t); print "| `" n "` | " c " | " m " | " t " |"; n = "" }')
rl=$(printf '%s\n' "$CFG" | sed -n "s/^ *resourceLimits = \[ cpus: \([0-9]*\), memory: '\([0-9]* GB\)', time: '\([0-9]*\)h' \]$/\1 CPUs, \2 and \3 h/p")
needr "caps every task at ${rl:-UNPARSED resourceLimits}"
hj_n=$(printf '%s\n' "$SUB" | sed -n 's/^#SBATCH -n //p'); hj_m=$(printf '%s\n' "$SUB" | sed -n 's/^#SBATCH --mem=//p'); hj_t=$(printf '%s\n' "$SUB" | sed -n 's/^#SBATCH -t //p')
needr "(\`-n $hj_n --mem=$hj_m -t $hj_t\`)"
needr "so it asks for more than the 4 h that this repository's skills use by default (a deliberate exception)"
forbidr "approved by the user"
forbidr "are NOT measured on a real genome"
# Correction 6: genome-BAM input; single-end and forward verified from code only.
needr "strandedness and read type are written to the BAM samplesheet (\`$(gv BAM_SHEET_HEADER)\`)"
needr "single-end BAM input and \`forward\` strandedness were verified from the pipeline code only, not by a run"
# Correction 7: unreleased pin, release 1.0.4, upstream rewrite, pin bump procedure.
needr "is a development commit (\`$(gv VERSION_TAG)\`), not a release"
needr "Release 1.0.4 does not launch under Nextflow $(gv NEXTFLOW_TESTED)"
needr "upstream PR #291"
needr "**Bumping the pin**"
for fx in rnasplice_schema.json config_params.txt trace_process_names.txt output_tree.txt command_lines.txt verified_urls.txt gate_values.tsv; do
  needr "\`$fx\`"
done
# Correction 8: the test runner (default and --full), the git requirement.
needr "bash tests/run_all_tests.sh --full"
needr "RNASPLICE_TEST_GIT_DIR"
needr "The default run takes about 17 minutes (bash, awk and sed only) and skips the proof-tool self-test (\`test_proof_tools.sh\`, 30 to 35 minutes)"
needr "\`bash tests/run_all_tests.sh --full\` includes it (about 50 minutes in all)."
forbidr "about 20 minutes"
# Correction 9: what was verified: the gate date, every run directory of RUN_DIRS, nf-core test data only.
needr "Verification gate ($(gv GATE_DATE))"
needr "nf-core's tiny test dataset"
needr "$(gv TEST_READ_LENGTH) bp reads"
needr "\`$(gv RUN_DIRS | sed 's/{.*$//')\`"
for rd in $(gv RUN_DIRS | sed -n 's/^.*{\(.*\)}$/\1/p' | tr ',' ' '); do
  needr "\`$rd\`"
done
# Correction 10: honest limitations.
needr "MISO (dropped: see Key design points)"
needr "MISO itself is unmaintained Python 2 software"
needr "there is no downstream report skill for the splicing tables yet (planned later)"
needr "No Salmon-results or transcriptome-BAM input, and no iGenomes \`genome\` key"
needr "Strandedness is asked, not detected"
needr "they were first used in the cluster acceptance run (run (a) below)."
forbidr "they are first used in the cluster acceptance run"
for k in ignore_tx_version miso_genes miso_read_len fig_height fig_width isoformswitchanalyzer_alpha isoformswitchanalyzer_dIF; do
  needr "\`$k\`"
done
needr "edgeR exon usage uses $(gv EDGER_DEU_FUNCTION) in this revision"
# Correction 8: run_all_tests.sh names every script of tests/ (TESTS, FULL_TESTS or NOT_TESTS), runs each under env -i, and
# says when the proof-tool self-test was skipped. CHECK_RUN_ALL overrides the file (used by prove_mutations.sh, target runall).
RUN_ALL=${CHECK_RUN_ALL:-$HERE/run_all_tests.sh}
if [ -s "$RUN_ALL" ]; then
  ra_tests=$(sed -n 's/^TESTS="\(.*\)"$/\1/p' "$RUN_ALL"); ra_full=$(sed -n 's/^FULL_TESTS="\(.*\)"$/\1/p' "$RUN_ALL"); ra_not=$(sed -n 's/^NOT_TESTS="\(.*\)"$/\1/p' "$RUN_ALL")
  for f in "$HERE"/*.sh; do
    b=${f##*/}
    case " $ra_tests $ra_full $ra_not " in *" $b "*) ;; *) echo "FAIL: tests/$b is not named in run_all_tests.sh (TESTS, FULL_TESTS or NOT_TESTS)"; fail=1 ;; esac
  done
  [ "$ra_not" = "cut_block.sh test_env.sh render_params.sh prove_red.sh run_all_tests.sh" ] \
    || { echo "FAIL: run_all_tests.sh NOT_TESTS must be exactly the helpers (cut_block.sh test_env.sh render_params.sh prove_red.sh run_all_tests.sh), found: $ra_not"; fail=1; }
  [ "$ra_full" = "test_proof_tools.sh" ] || { echo "FAIL: run_all_tests.sh FULL_TESTS must be test_proof_tools.sh, found: $ra_full"; fail=1; }
  grep -qF 'env -i HOME="${HOME:-/nonexistent}" PATH="$P" /bin/bash --noprofile --norc "$HERE/$1"' "$RUN_ALL" \
    || { echo "FAIL: run_all_tests.sh must run every test under env -i with an explicit PATH"; fail=1; }
  grep -qF 'test_proof_tools.sh was SKIPPED' "$RUN_ALL" || { echo "FAIL: run_all_tests.sh must say when the proof-tool self-test was skipped"; fail=1; }
  # Task 8 review I3: the version comparison itself (its behaviour is tested by test_run_all.sh on a fake tree).
  grep -qF "1.8.5 | sort -V | head -n1)\" != 1.8.5 ]; then" "$RUN_ALL" || { echo "FAIL: run_all_tests.sh must check for git 1.8.5 or newer before the proof-tool self-test"; fail=1; }
  case " $ra_tests " in *" test_run_all.sh "*) ;; *) echo "FAIL: run_all_tests.sh TESTS must run test_run_all.sh (the test of the runner itself)"; fail=1 ;; esac
else
  echo "FAIL: tests/run_all_tests.sh is missing or empty"; fail=1
fi
# --- end Task 8

# --- Task 8 review fixes (I1, I2, I4, minors)
# I1: the pipeline's labels take precedence over the process default; the maximum they requested comes from the gate trace.
forbidr "gives every process 2 CPUs, 8 GB and 4 h"
tr_max=$(awk -F'\t' '{ c = $2 + 0; m = $3; sub(/ GB$/, "", m); t = $4; sub(/h$/, "", t)
  if (c > mc) mc = c; if (m + 0 > mm) mm = m + 0; if (t + 0 > mt) mt = t + 0 } END { print mc " CPUs, " mm " GB and " mt " h" }' "$FIX/trace_resources_g3.tsv")
needr "Every other process keeps the resources of the pipeline's own process labels, which take precedence over that default: in the verification run they requested up to $tr_max"
needr "sets a process default of 2 CPUs, 8 GB and 4 h, raises the processes below with \`withName\` selectors"
# I2: the SALMON_QUANT.* selector was never applied in a run.
forbidr "they were checked on the nf-core test data only"
# (real-data fix round: the selector values changed after the real-data test; their evidence is pinned in that section below)
needr "With this skill's route no SALMON_QUANT_STAR task runs; the pattern covers both names."
forbidr "selector was only matched against the recorded process names"
forbidr "is first applied in the cluster acceptance run"
# I4: honesty statements (each deletion or inversion fails).
forbidr "and nothing in this skill is validated on real data."
forbidr "has been validated on real"
forbidr "validated on real RNA-seq data"
forbidr "tested on real data"
forbidr "likely takes more than an hour to build (unverified)."
needr "The pipeline's decoy-aware Salmon index of the human genome took 33 min (peak 19.8 GB) in the real-data test; a mouse genome was not measured, and STAR_GENOMEGENERATE was not measured at all (the real-data test reused a STAR index)."
needr "an unreleased development commit of the pipeline, not a release"
forbidr "the current release"
forbidr "mean PSI of the treatment − mean PSI of the control"
forbidr "mean PSI of the treatment - mean PSI of the control"
needr "| rMATS | \`IncLevelDifference\` = mean(\`IncLevel1\`) − mean(\`IncLevel2\`); \`b1\` = the treatment samples | more inclusion in the treatment | verification run (bamlists and values) |"
needr "| more inclusion in the control | verification run (every event, values) and pinned code |"
# M3: the DTU and DEXSeq exon directions were checked on every non-NA row (counts recorded in gate_report.md).
n_dtu=$(sed -n 's/^.*all \([0-9][0-9]*\) non-NA rows have sign(log2fold_YRI_GBR) = sign(YRI - GBR) of the fitted columns.*$/\1/p' "$FIX/gate_report.md")
n_dex=$(sed -n 's/^.*DEXSeq exon usage values .*all \([0-9][0-9]*\) non-NA rows have sign(log2fold_YRI_GBR) = sign(YRI - GBR)\.$/\1/p' "$FIX/gate_report.md")
[ -n "$n_dtu" ] && [ -n "$n_dex" ] || { echo "FAIL: gate_report.md must record the DTU and DEXSeq exon sign checks (all N non-NA rows)"; fail=1; }
needr "| DEXSeq DTU | \`log2fold_{CONTROL}_{TREATMENT}\` | a larger share of the transcript in the control | verification run (header plus checked rows: all ${n_dtu:-N} non-NA rows) and pinned code |"
needr "| DEXSeq exon usage | \`log2fold_{CONTROL}_{TREATMENT}\` | more usage of the exon bin in the control | verification run (header plus checked rows: all ${n_dex:-N} non-NA rows) and pinned code |"
grep -qF "Merge pull request #291" "$FIX/gate_report.md" || { echo "FAIL: gate_report.md must record the evidence of upstream PR #291"; fail=1; }
needr "until it was killed by hand (the retry then finished in 2 s); cause not verified."
needr "- One contrast was run with the skill's settings; several contrasts only in the pipeline's own test profile."
needr "- With FASTQ input, Salmon (index and quantification) always runs, also when neither DTU nor SUPPA2 is chosen."
needr "With FASTQ input, Salmon (index build and quantification) runs even in an rMATS-only run"
needr "(the gate ran paired-end samples only: \`unstranded\` and \`reverse\`, from FASTQ and from genome BAM); single-end FASTQ input was not run either."
forbidr "\`star_index\` (reuse of an existing STAR index) was not run either."
needr "\`star_index\` (reuse of an existing STAR index) was first run in the real-data test."
needr "): run it in one place only."
forbidr "run it in both places"
forbidr "no more than 4 h"
# Minors: Step 6 stop, prefix, run directories, option keys, PR #291 wording, the runner test.
needr "warns when read lengths differ, and stops by default when they differ between conditions |"
forbidr "{prefix}"
needr "\`{SHEET_PREFIX}\` is the file prefix chosen in Step 5"
needr "(on the author's cluster account; these run directories are not distributed, and the recorded evidence is in \`tests/fixtures/\`)"
needr "every option key that has a non-null pipeline default"
forbidr "every option of every module"
needr "which ported the rMATS subworkflow to the nf-core structure"
needr "a test of the test runner itself on a fake tree (\`test_run_all.sh\`"
# M8: no row may be added to the README tables: each has exactly the rows the skill defines.
tbl_rows() { awk -v h="$1" 'index($0, h) == 1 { f = 1; getline; next } f && /^\|/ { n++; next } f { exit } END { print n + 0 }' "$README" 2>/dev/null; }
want_sel=$(printf '%s\n' "$CFGS" | grep -c "withName:")
want_ho=$(printf '%s\n' "$HO" | grep -c '^  [^ ]\{1,\}  ')
want_steps=$(grep -cE '^## Step [0-9]+ ' "$SKILL")
for spec in "| Selector | CPUs | Memory | Time |:$want_sel" "| Path | Content |:$want_ho" "| Output | Column | Positive value means | Source |:5" \
            "| File | Description |:8" "| Step | Topic | Asked or detected |:$want_steps"; do
  h=${spec%:*}; w=${spec##*:}; got=$(tbl_rows "$h")
  [ "$got" = "$w" ] || { echo "FAIL: README table '$h' has $got rows, expected $w"; fail=1; }
done
# --- end Task 8 review fixes
# --- Task 9 fix round (D1-D4, U1, acceptance record)
# D1: an R2 mate is never a sample; single-end only without an R1 partner (fastq_layout, tested by test_sample_names.sh).
need "fastq_layout() {"
need "that R2 file is its mate, never a sample of its own. Only a file that is neither such an R1 file nor the R2 mate of one is single-end"
need "- **Rows:** one samplesheet row per \`paired\` or \`single\` line of \`fastq_layout\`"
forbid "and whose R2 file exists → paired-end; otherwise single-end."
# D2: BAM sample names lose .umi_dedup.sorted.bam, .markdup.sorted.bam, .sorted.bam, _sorted.bam or .bam (bam_sample_name).
need "bam_sample_name() {"
need "derive the sample names with \`bam_sample_name\` (block in Step 4a: it removes the longest of"
# D3: identical file-name / title options collapse to one.
need "When \`{SEQ_DATE}\` equals \`{TODAY_YYMMDD}\` (the input files have no date prefix), options 1 and 2 are the same: offer only 1. \`{TODAY_YYMMDD}_{WD_NAME}\` · 2. Custom prefix."
need "When \`{SEQ_DATE}\` equals \`{TODAY_YYMMDD}\`, options 1 and 2 are the same: offer only 1. \`{TODAY_YYMMDD}_{WD_NAME}\` · 2. Custom."
# D4: BAM path style as for FASTQ.
need "**Path style:** if \`{BAM_DIR}\` is inside \`{CWD}\`, use paths relative to \`{CWD}\`; otherwise absolute paths (as for FASTQ input)."
# Acceptance record (task-9-report.md numbers): three of them pinned, plus the not-exercised sentence.
needr "job 11380326, wall time 11 min 49 s, 80/80 tasks COMPLETED"
needr "job 11380459, wall time 6 min 16 s, 43/43 tasks COMPLETED"
needr "the rMATS sign on 202/202 rows of SE.MATS.JC.txt and the SUPPA2 sign on 1098/1098 local and 1854/1854 isoform events"
needr "- Not exercised: real data; the Ensembl download helper; a paired design; strandedness other than unstranded, and single-end input; a human- or mouse-size genome; the time and memory of any process on real data."
# --- end Task 9 fix round
# --- Final review fixes (I1, I3, M1, M2, M4, M6)
# I1: custom prefix and custom file names: Step 10's character rule; install_rnasplice_sheets refuses a split or odd name.
need "A custom prefix (and every custom file name of this step) may contain only letters, digits, \`.\`, \`_\` and \`-\`"
need "2. choose another filename (letters, digits, \`.\`, \`_\` and \`-\` only, as for the custom prefix;"
need "[ \$# -le 6 ] || { echo \"ERROR: \$# arguments, at most 6: a file name was split at a space"
need "''|*/*|*[!A-Za-z0-9._-]*) echo \"ERROR: target name '\$p' must be a plain file name in the current directory (letters, digits, ., _ and - only)\""
# I3: no tool that must run on a compute node is invoked in the wizard text outside the submission script (as for nextflow).
# Salmon may also run inside the strandedness helper script (Step 4, a compute-node job; its block is cut here).
STRH=$(bash "$CUT" "$SKILL" "**Strandedness helper script.**" 2>/dev/null)
for spec in 'conda (activate|create|install|run|env)( |$)' '(Rscript|R +(-e|--vanilla|-f|CMD))( |$)' 'STAR +--?[A-Za-z]' 'salmon +(index|quant|--?[A-Za-z])'; do
  pat="(^|[^A-Za-z_])$spec"; also=""
  n_all=$(grep -cE "$pat" "$SKILL"); n_sub=$(printf '%s\n' "$SUB" | grep -cE "$pat")
  case "$spec" in salmon*) n_sub=$((n_sub + $(printf '%s\n' "$STRH" | grep -cE "$pat"))); also=" and the strandedness helper" ;; esac
  [ "$n_all" -eq "$n_sub" ] || { echo "FAIL: '$spec' is invoked outside the submission script ($n_all lines in the skill, $n_sub in the script$also); the wizard must never run it on the login node"; fail=1; }
done
# M1: rMATS prep runs per sample at this revision (gate trace: 4 prep tasks, 1 post task).
need "At this revision rMATS prepares each sample once and runs one post step per contrast"
forbid "rMATS runs one prep/post pair per contrast"
# M2: the explicit fastq_layout call, unmatched patterns, subdirectories, bash 4.
need "\`shopt -s nullglob; fastq_layout \"{FASTQ_DIR}\"/*.fastq.gz \"{FASTQ_DIR}\"/*.fq.gz\`"
need "Only the files directly in \`{FASTQ_DIR}\` are used, not its subdirectories"
need "it needs bash 4 or newer"
need "is not a file (an unmatched pattern? call with shopt -s nullglob)"
# M4: contig names.
need "The FASTA and the GTF must use the same contig names"
need "tell the user to compare the \`@SQ\` names of \`samtools view -H\` of one BAM with the FASTA headers and the first column of the GTF, on a compute node"
needr "- FASTA, GTF and BAM contig names must match"
# M6: executor throttles of the config template.
printf '%s\n' "$CFG" | grep -qxF '            queueSize = 20' && printf '%s\n' "$CFG" | grep -qxF "            submitRateLimit = '10/1min'" \
  || { echo "FAIL: the nextflow.config template must keep queueSize = 20 and submitRateLimit = '10/1min'"; fail=1; }
# I2, M7: what the acceptance runs cover at the final skill.
needr "The runs and the static gate below were made on the skill text of commit bcb63e9."
needr "\`nextflow.config\` is byte-identical to the Step 10 template of that commit, and \`nf-core_rnasplice_dev-1b44723.sh\` and \`261002_{a,b,c}_params.yaml\` match the Step 11 templates (with the Step 8 block) line for line"
needr "- Static gate (commit bcb63e9):"
needr "its automatic sample names still ended in \`_sorted\` (made before Step 4b strips \`_sorted.bam\`; that fix is unit-tested, run (c) was not repeated)"
needr "| bash 4 or newer |"
# --- end Final review fixes

# --- Real-data fix round (real-data test 2026-10-02, human hiPSC-CM, 6 samples)
# D1: the STAR_ALIGN ext.args override of the top-level process block = the STAR_ALIGN ext.args block of the pinned
# conf/modules.config (fixture modules.config, byte-identical to the pinned commit) minus '--quantMode TranscriptomeSAM',
# '--quantTranscriptomeSAMoutput BanSingleEnd' and comment lines. Computed from the fixture, so a pin bump cannot drift silently.
QM="                '--quantMode TranscriptomeSAM',"; QB="                '--quantTranscriptomeSAMoutput BanSingleEnd',"
mc_block=$(awk '/^    withName: .STAR_ALIGN. \{$/ {f = 1} f {print} f && /^        \}$/ {exit}' "$FIX/modules.config")
{ [ "$(printf '%s\n' "$mc_block" | grep -cxF -- "$QM")" -eq 1 ] && [ "$(printf '%s\n' "$mc_block" | grep -cxF -- "$QB")" -eq 1 ] \
  && [ "$(printf '%s\n' "$mc_block" | sed -n 2p)" = "        ext.args   = {" ]; } \
  || { echo "fixture modules.config: the STAR_ALIGN ext.args block was not found, or lacks one of the two transcriptome-BAM flags"; exit 2; }
want_star=$(printf '%s\n' "$mc_block" | grep -vxF -- "$QM" | grep -vxF -- "$QB" | grep -vE '^[ \t]*//'; echo "    }")
got_star=$(printf '%s\n' "$CFG" | awk '/^process \{$/ {p = 1; next} p && /^\}$/ {p = 0} p && /^    withName: .STAR_ALIGN. \{$/ {f = 1} f {print} f && /^    \}$/ {exit}')
if [ -z "$got_star" ]; then
  echo "FAIL: the nextflow.config template must hold the STAR_ALIGN ext.args override (withName: 'STAR_ALIGN' in the top-level process block)"; fail=1
elif [ "$got_star" != "$want_star" ]; then
  echo "FAIL: the STAR_ALIGN ext.args override differs from the pinned conf/modules.config block minus --quantMode TranscriptomeSAM and --quantTranscriptomeSAMoutput BanSingleEnd:"
  diff <(printf '%s\n' "$want_star") <(printf '%s\n' "$got_star") | sed 's/^/    /'; fail=1
fi
# The selector string is the pinned one, used once; every other withName block is a resource selector of the slurm profile.
[ "$(printf '%s\n' "$CFG" | grep -c "withName: 'STAR_ALIGN'")" -eq 1 ] || { echo "FAIL: withName: 'STAR_ALIGN' must appear exactly once in the nextflow.config template"; fail=1; }
[ "$(printf '%s\n' "$CFG" | grep -c 'withName:')" -eq $(( $(printf '%s\n' "$CFGS" | grep -c 'withName:') + 1 )) ] \
  || { echo "FAIL: every withName block of the template except the STAR_ALIGN ext.args override belongs to the slurm profile"; fail=1; }
! printf '%s\n' "$CFG" | grep -qE 'quantMode|quantTranscriptomeSAMoutput' || { echo "FAIL: the nextflow.config template must not set --quantMode or --quantTranscriptomeSAMoutput"; fail=1; }
# Under this skill's params (seq_center and save_unaligned are never written: pipeline defaults null and false) the list
# evaluates to its quoted literals: they must equal the STAR arguments of all six STAR_ALIGN tasks of the real-data run
# (fixture realdata_star_align_args.txt, copied from their .command.sh).
[ "$(cfg_default seq_center)" = null ] && [ "$(cfg_default save_unaligned)" = false ] || { echo "config extraction broken: seq_center or save_unaligned default changed"; exit 2; }
for k in seq_center save_unaligned; do
  ! printf '%s\n%s\n' "$BASEY" "$MOD" | grep -q "^$k:" || { echo "FAIL: $k must not be written (the STAR_ALIGN override is checked for its default)"; fail=1; }
done
star_eval=$(printf '%s\n' "$got_star" | sed -n "s/^ *'\(--[^']*\)',\{0,1\}$/\1/p" | tr '\n' ' ' | sed 's/ $//')
[ "$star_eval" = "$(head -n 1 "$FIX/realdata_star_align_args.txt")" ] \
  || { echo "FAIL: the STAR_ALIGN override evaluates to '$star_eval', not to the arguments of the real-data run ($(head -n 1 "$FIX/realdata_star_align_args.txt"))"; fail=1; }
need "**STAR without the transcriptome BAM.**"
need "The \`withName: 'STAR_ALIGN'\` block of the template therefore repeats the pinned \`ext.args\` without \`--quantMode TranscriptomeSAM\` and \`--quantTranscriptomeSAMoutput BanSingleEnd\`"
need "it is part of \`nextflow.config\`, so the launch line does not change"
need "the same selector and arguments, given in a separate config file, were used for the resumed run of the real-data test (22 tasks cached, 0 failed)"
need "Not verified: an unmodified STAR_ALIGN run to completion at full size"
need "With BAM input no STAR_ALIGN task runs, so the block changes nothing there."
need "the \`STAR_ALIGN\` \`ext.args\` block below) so they can compare"
needr "- **No STAR transcriptome BAM.**"
needr "repeats the pinned STAR_ALIGN \`ext.args\` without \`--quantMode TranscriptomeSAM\` and \`--quantTranscriptomeSAMoutput BanSingleEnd\`; the launch line is unchanged."
needr "Verified on the real-data test with the same selector and arguments in a separate config file (resumed run: 22 tasks cached, 0 failed)"
needr "the form written into \`nextflow.config\` is to be verified by an acceptance re-run on the nf-core test data"
needr "Not verified: an unmodified STAR_ALIGN run to completion at full size. BAM input is unaffected (no STAR_ALIGN task)."
needr "The \`nextflow.config\` template was changed after the real-data test (see **Real-data test** below); these three runs did not use its new form."
forbidr "\`nextflow.config\` is byte-identical to the Step 10 template, and"
# Tiers (SEL_ROWS above, evidence = realdata_resource_table.txt), queueSize 20, trace fields (D4), DEXSEQ_EXON peak_rss (D7).
forbid "The resources of every selector are judgement: they are not measured on real data."
forbid "It requests 8 CPUs, 36 GB and 8 h"
need "The SUPPA2 selector \`'.*:DIFFSPLICE_IO[EI]'\` matches its two differential-splicing processes, DIFFSPLICE_IOE and DIFFSPLICE_IOI."
need "**Where the selector values come from.** \`STAR_GENOMEGENERATE\` and \`DEXSEQ_EXON\` keep the values that this skill's verification run on the nf-core test data applied."
need "one human dataset (6 samples of about 38 M read pairs of 78 bp, one contrast, all five analyses), so they rest on that one dataset, and these new values have not yet been applied in a run."
need "- STAR_ALIGN: peak 37.7 GB (79% of 48 GB; the human index takes most of it), at most 15 min: 8 CPUs, 48 GB, 4 h."
need "- RMATS_PREP: single-threaded, 1.3 GB, 8 min: 2 CPUs, 8 GB, 4 h."
need "- RMATS_POST: 2.0 GB, 3 min for one contrast: 8 CPUs, 16 GB, 8 h"
need "- DEXSEQ_COUNT: single-threaded, 1.15 GB, 55 min: 1 CPU, 4 GB, 8 h"
need "- DEXSEQ_DTU: 11.7 GB, 2 min: 8 CPUs, 16 GB, 4 h."
need "- SALMON_QUANT_SALMON: 20.2 GB with the pipeline's decoy-aware human index, 13.5 min: 8 CPUs, 32 GB, 4 h."
need "- DIFFSPLICE_IOE/IOI (SUPPA2): single-threaded, 3.2-4.2 GB, up to 3 h 03 min (the longest task of the run): 1 CPU, 8 GB, 16 h"
need "- MAKE_TRANSCRIPTS_FASTA: 4.2 GB, 2 min: 2 CPUs, 8 GB, 2 h"
need "- DEXSEQ_EXON: 83 min; its memory is not known (see below): 8 CPUs, 32 GB, 8 h kept."
need "- STAR_GENOMEGENERATE: not run in the real-data test (the STAR index was reused): 8 CPUs, 64 GB, 8 h kept."
need "\`queueSize\` is 20 because with 10, in the real-data test, RMATS_POST waited 11 min and two Salmon tasks about 13 min for a free slot"
need "\`peak_rss\`, \`realtime\` and \`%cpu\` what it used, and \`attempt\` whether it was retried"
need "**DEXSEQ_EXON memory:** its \`peak_rss\` in the trace is misleading. In the real-data test the trace reported 187 GB for this task, which requested 32 GB and was not killed; \`sstat\` showed 8.4 GB during its serial phase."
need "Its true peak was not measured: do not raise its selector because of that column alone; compare with SLURM's MaxRSS"
needr "| Selector | CPUs | Memory | Time | Evidence |"
needr "The values marked \"real-data test\" were set from what each process used in the real-data test (below): ONE human dataset, 6 samples of about 38 M read pairs of 78 bp, one contrast, all five analyses. They rest on that one dataset, and these values have not yet been applied in a run"
needr "\`STAR_GENOMEGENERATE\` and \`DEXSEQ_EXON\` keep the values of the verification gate."
needr "The executor's \`queueSize\` is 20"
needr "for \`DEXSEQ_EXON\`, the trace's \`peak_rss\` is misleading (187 GB reported in the real-data test for a 32 GB task that was not killed; \`sstat\` showed 8.4 GB; probably forked workers counted separately): compare with SLURM's MaxRSS instead."
# each README selector row names its evidence: real-data test for 'realdata' rows, the verification gate for 'gate' rows
while IFS=$'\t' read -r re c m t; do
  [ -n "$re" ] || continue
  n=${re#.\*:}
  case "$(sel_ev "$re")" in
    realdata) needr "| \`$n\` | $c | $m | $t | real-data test: " ;;
    gate) needr "| \`$n\` | $c | $m | $t | verification gate" ;;
  esac
done <<< "$sels"
# D3: strandedness from a Salmon subsample (FASTQ input): an optional compute-node helper; the parser uses the RSeQC rule's
# thresholds (tested by test_strand_salmon.sh); the index is never passed to the pipeline; a kit name never decides.
for a in "**Strandedness from a Salmon subsample** (FASTQ input; optional)." "**Salmon index for the strandedness helper.**" "**Strandedness helper script.**" "**Reading the Salmon results.**"; do anchor_once "$a"; done
need "A kit name alone does not determine the direction"
need "in this skill's real-data test the GEO text said \"NEB Ultra II protocol\", without \"Directional\", and the reads were clearly reverse-stranded: Salmon found ISR), so never derive the answer from a kit name"
forbid_re 'Ultra II[^.]*(implies|means|so it is) (stranded|unstranded|reverse|forward)' "a strandedness derived from a kit name"
need "without such a directory, with FASTQ input, offer the Salmon helper"
need "only when a Salmon index of the same organism exists in the shared genome folder"
need "The index serves only this detection: the pipeline never receives it (\`{SALMON_INDEX}\` stays empty, Step 7), and the wizard never guesses from the result: it proposes, and the user confirms."
need "No \`SALMON_INDEX\` line: the helper is not offered"
need "at least two, one from each condition the user will compare"
need "if Step 4 already asked the organism and this directory for the strandedness helper, reuse both answers and do not ask again."
need "applies the rule of the RSeQC block above to Salmon's counts: forward fraction = ISF / (ISF + ISR + IU), reverse fraction = ISR / (ISF + ISR + IU) for paired-end reads (SF, SR and U for single-end reads); forward if the forward fraction is at least 0.8, reverse if the reverse fraction is at least 0.8, unstranded if the two fractions differ by at most 0.1, otherwise unclear."
need "The result is unclear as well when it disagrees with the library type Salmon detected itself"
need "when Salmon's compatible-fragment ratio is below 0.8"
need "any \`unclear\` result, or samples that disagree: show the lines, explain that the library type cannot be read reliably, and ask the user (no proposal); do not continue until they answer 1, 2 or 3"
need "salmon_strandedness() {"
need "verify the container URL of the script in this session with a HEAD request"
STRH_FULL=$(bash "$CUT" "$SKILL" "**Strandedness helper script.**" 2>/dev/null)
[ -n "$STRH_FULL" ] || { echo "FAIL: strandedness helper script block not found"; fail=1; }
for l in '#!/bin/bash' '#SBATCH -N 1' '#SBATCH -n 8' '#SBATCH --mem=32G' '#SBATCH -t 2:00:00' '#SBATCH -p bcc' '#SBATCH --mail-type=END,FAIL' '#SBATCH --mail-user={USER_EMAIL}' \
         '#SBATCH -o infer_strandedness_salmon_{WD_NAME}.%j.log' 'set -u' 'IDX="{STRAND_INDEX}"' \
         'module add singularity/3.10.4 || { echo "ERROR: cannot load module singularity/3.10.4" >&2; exit 1; }' \
         'command -v singularity >/dev/null || { echo "ERROR: singularity is not on PATH" >&2; exit 1; }' \
         'run_sample "{SAMPLE}" "{FASTQ_1}" "{FASTQ_2}"'; do
  printf '%s\n' "$STRH_FULL" | grep -qxF -- "$l" || { echo "FAIL: the strandedness helper must contain the line: $l"; fail=1; }
done
[ "$(printf '%s\n' "$STRH_FULL" | grep -v '^#' | grep -v '^[ \t]*$' | sed -n 2p)" = 'cd "{CWD}" || { echo "ERROR: cannot change to {CWD}" >&2; exit 1; }' ] \
  || { echo "FAIL: cd \"{CWD}\" must be the first command of the strandedness helper (after set -u)"; fail=1; }
! printf '%s\n' "$STRH_FULL" | grep -qE '(^|[^A-Za-z_])rm( |$)' || { echo "FAIL: the strandedness helper must not remove files"; fail=1; }
! printf '%s\n' "$STRH_FULL" | grep -qE 'nextflow|conda' || { echo "FAIL: the strandedness helper must not run nextflow or conda"; fail=1; }
[ "$(printf '%s\n' "$STRH_FULL" | grep -c 'salmon quant -i "\$IDX" -l A ')" -eq 2 ] || { echo "FAIL: the strandedness helper must run salmon quant -i \"\$IDX\" -l A for paired-end and single-end samples"; fail=1; }
printf '%s\n' "$STRH_FULL" | grep -qF 'head -n 4000000' || { echo "FAIL: the strandedness helper must subsample the first 1,000,000 reads (head -n 4000000)"; fail=1; }
# The container: the pipeline's own Salmon (gate SALMON_VERSION), a URL verified at the gate, and its cache file name.
s_url=$(printf '%s\n' "$STRH_FULL" | grep -oE 'https://depot\.galaxyproject\.org/singularity/salmon:[^"]+' | sort -u)
[ "$(printf '%s\n' "$s_url" | grep -c .)" -eq 1 ] && case "$s_url" in *"salmon:$(gv SALMON_VERSION)--"*) true ;; *) false ;; esac \
  || { echo "FAIL: the strandedness helper must download one Salmon $(gv SALMON_VERSION) container (found: $s_url)"; fail=1; }
grep -qE "^$(printf '%s' "$s_url" | sed 's/[.]/\\./g') 200 " "$FIX/verified_urls.txt" || { echo "FAIL: the Salmon container URL $s_url is not a verified URL of fixtures/verified_urls.txt"; fail=1; }
s_sif=$(printf '%s' "${s_url#https://}" | tr '/:' '--').img
printf '%s\n' "$STRH_FULL" | grep -qF "SIF=\"\${NXF_SINGULARITY_CACHEDIR:-\$HOME/.singularity/cache}/$s_sif\"" \
  || { echo "FAIL: the strandedness helper must look for the cached image $s_sif in the Singularity cache"; fail=1; }
needr "a Salmon helper that detects the library type from 1,000,000 reads of a few samples"
needr "A kit name alone does not determine the direction"
needr "The Salmon helper was run by hand in the real-data test (reverse-stranded paired-end data); its generated form was tested with stubs only."
# D2: disk estimate (report D2: about 40 GB + 18 GB per sample with the STAR_ALIGN override, about 38 GB per sample without it),
# where the output directory is asked and in the README prerequisites; work/ is never deleted by the skill.
need "**Disk space.** When asking for the output directory, tell the user how much space the run needs in \`{CWD}\`"
need "a run needed about 40 GB plus 18 GB per sample for \`work/\` and the output directory together, with the \`STAR_ALIGN\` block of the \`nextflow.config\` template (about 38 GB per sample without it), not counting the FASTQ files"
need "This is an estimate from one dataset; it grows with the number of reads per sample."
need "It is needed to resume a run (**Resuming a run**, Step 11); after a successful run the user can delete it to free the space. The wizard and the generated scripts never delete it."
forbid_re 'rm +-[a-z]*r[a-z]* +[^ ]*work' "a command that deletes work/ (the user deletes it, never the wizard)"
needr "| Disk space | About 40 GB plus 18 GB per sample for \`work/\` and the results together (FASTQ files not counted), measured in the real-data test on one human dataset of about 38 M read pairs of 78 bp per sample"
needr "(about 38 GB per sample without it); it grows with the reads."
needr "resuming a run needs it, and after a successful run it can be deleted; the skill never deletes it |"
# D5: stopping a run without orphan tasks, and without cancelling other runs.
need "- To stop a run: cancel the head job with \`scancel {JOBID}\`."
need "Check each one with \`scontrol show job <id>\` (its \`WorkDir\` must be under \`{CWD}/work\`, so that it belongs to this run), and cancel them one at a time with \`scancel <id>\`."
need "Never cancel by a name pattern or all of the user's jobs at once: that can stop other runs."
forbid_re 'scancel +(-[A-Za-z]|--[a-z])' "a scancel command that selects jobs by user or name (it can stop other runs)"
# D8: the rMATS chr column.
need "With an Ensembl genome, rMATS writes the \`chr\` column with a \`chr\` prefix that the FASTA and the GTF do not have"
need "strip it before joining the rMATS tables with the GTF or with other tables."
# D6: edgeR direction now checked numerically on one exon of the real-data test (report §6(4)).
need "checked numerically on one exon of the real-data test: the RBPMS2 exon skipped in the knockout has logFC −1.59 for KO vs WT"
forbid "which builds the contrast as treatment-control; not checked numerically"
# Item 5: the README record of the real-data test; numbers from the real-data report; honesty (one dataset, the 2.6x difference
# from the paper not explained), deviations and what was not verified.
grep -qE '^Real-data test \(one human dataset\): DONE \(2026-[0-9-]+\), with the deviations listed below$' "$README" 2>/dev/null \
  || { echo "FAIL: README must have the whole line 'Real-data test (one human dataset): DONE (2026-MM-DD), with the deviations listed below'"; fail=1; }
needr "This is one dataset, not a validation on real data in general."
needr "Akerberg et al. 2022 (GEO GSE207681)"
needr "WT_1, WT_2, WT_3 = SRR20021261, SRR20021260, SRR20021259 (GSM6307608-GSM6307610) and KO_1, KO_2, KO_3 = SRR20021257, SRR20021256, SRR20021255 (GSM6307612-GSM6307614)"
needr "- What ran: all five analyses (rMATS, SUPPA2, DEXSeq and edgeR exon usage, DEXSeq DTU), one contrast (\`KO_vs_WT\`), no paired design."
needr "81 tasks completed and 22 cached, 0 failed, 0 retried; the resumed run took 4 h 11 min"
needr "RBPMS2 expression (Salmon gene TPM) was 67% lower in KO"
needr "was found by all four splicing analyses"
needr "PC1 (65.9% of the variance) separated WT from KO"
needr "rMATS found 7,031 significant events (junction-count tables), about 2.6 times the 2,679 events the paper reports. The reason was not determined"
needr "every sign statement of Step 12 held on real rows"
needr "edgeR (logFC = KO − WT) on one exon"
needr "- Deviations of the test: the first head job was cancelled after 22 min"
needr "seven of its task jobs, still running, were cancelled by hand"
needr "The run was resumed with a separate config file (\`-c\`) holding the STAR_ALIGN \`ext.args\` override"
needr "dead intermediate files (the killed task directories, the unsorted BAM files, the STAR index copy) were deleted by hand during the run."
needr "- Not verified by this test: STAR_GENOMEGENERATE (the index was reused); an unmodified STAR_ALIGN run to completion; the true memory peak of DEXSEQ_EXON; the new selector values and the new \`nextflow.config\` template (the run used the earlier values and a separate file); a paired design; single-end input; \`forward\` strandedness; genome-BAM input; the Ensembl download helper; the generated Salmon strandedness helper"
forbidr "validated on several"
forbidr "validated on multiple"
forbidr "fully validated"
forbidr "reproduces the paper's"
forbidr "matches the paper's"
forbidr "explained by the sample number"
needr "the stub dry runs of the submission script, the download helper and the Salmon strandedness helper"
# --- end Real-data fix round
[ $fail -eq 0 ] && echo "PASS" || exit 1
