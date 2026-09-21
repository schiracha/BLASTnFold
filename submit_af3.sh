#!/usr/bin/env bash
# =============================================================================
# submit_af3.sh
# =============================================================================
# PURPOSE
#   Orchestrates the full AF3 → USalign pipeline by submitting each stage as
#   a SLURM job array with automatic sizing and dependency chaining.
#   Run directly on the login node: bash submit_af3.sh ...
#   Do NOT sbatch this script.
#
# CONCURRENCY THROTTLES  (the %N in --array=1-X%N)
#   data      : default 16  — CPU/memory only; run as many as the cluster allows
#   inference : default 4   — GPU-constrained; limit concurrent GPU jobs
#   align     : default 16  — CPU only, lightweight; high concurrency is fine
#
#   Override with --throttle-data, --throttle-infer, --throttle-align.
#   Set to 0 to remove the throttle entirely (SLURM will run all tasks at once).
#
# ARRAY SIZING STRATEGY
#   Priority 1: count actual files on disk  (used when resuming with --from)
#   Priority 2: estimate from input file line count (requires --chain-estimate)
#   If neither is possible, the script stops and prints the re-run command.
#
#   Recommended two-step workflow for a fresh run:
#     Step 1:  bash submit_af3.sh INPUT --mode accession [--source-tag TAG]
#              → submits make_json; prints the step-2 command
#     Step 2:  bash submit_af3.sh --from data --chain-min N --chain-max N
#              → counts finished JSONs, submits data+inference+align+summarize
#
#   One-shot (uses estimates — array sizes may be slightly off if sequences
#   are missing from the DB, which SLURM handles gracefully):
#     bash submit_af3.sh INPUT --mode accession --chain-estimate \
#         --chain-min 11 --chain-max 11
#
# PIPELINE STAGES
#   1. make_json    make_json_unified.sbatch   → AF3_Inputs/*_1mer.json
#   2. data         AF3_data.sbatch            → AF3_Data/*_Nmer_data.json
#   3. inference    AF3_inference.sbatch       → AF3_Inference/*/*_N.cif
#   4. align        US_align.sbatch            → Align_Out/*_stats.out
#   5. summarize    summarize.sbatch           → *_SummaryHits.csv
#
# USAGE
#   bash submit_af3.sh INPUT_FILE [OPTIONS]
#   bash submit_af3.sh --from STAGE [OPTIONS]
#
# OPTIONS
#   --from  STAGE           resume from this stage (default: make_json)
#                           valid: make_json | data | inference | align | summarize
#   --mode  MODE            input flavour for make_json (fasta|region|accession|hits)
#                           "regions" is accepted as an alias for "region"
#   --source-tag  TAG       provenance label for manifest
#   --chain-min  N          minimum oligomer size (default: 5)
#   --chain-max  N          maximum oligomer size (default: 15)
#   --min-tm-avg  F         hits mode TM threshold (default: 0.2)
#   --skip-existing         pass --skip-existing to make_json
#   --chain-estimate        size arrays from input file line count (one-shot mode)
#   --throttle-data   N     max concurrent data tasks    (default: 16, 0=unlimited)
#   --throttle-infer  N     max concurrent inference tasks (default: 4,  0=unlimited)
#   --throttle-align  N     max concurrent align tasks   (default: 16, 0=unlimited)
#   --dry-run               print sbatch commands without submitting
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
INPUT_FILE=""
FROM_STAGE="make_json"
MODE=""
SOURCE_TAG=""
CHAIN_MIN=5
CHAIN_MAX=15
MIN_TM_AVG=0.2
SKIP_EXISTING_FLAG=""
CHAIN_ESTIMATE=false
DRY_RUN=false
LIPID_TOKENS=258   # default ~250 tokens: 2×(POPC+POPE+CLR) = 258 tokens; set 0 to disable

# Concurrency throttles (0 = no throttle)
THROTTLE_DATA=16
THROTTLE_INFER=4
THROTTLE_ALIGN=16

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --from)            FROM_STAGE="$2";        shift 2 ;;
    --mode)            MODE="$2";              shift 2 ;;
    --source-tag)      SOURCE_TAG="$2";        shift 2 ;;
    --chain-min)       CHAIN_MIN="$2";         shift 2 ;;
    --chain-max)       CHAIN_MAX="$2";         shift 2 ;;
    --min-tm-avg)      MIN_TM_AVG="$2";        shift 2 ;;
    --skip-existing)   SKIP_EXISTING_FLAG="--skip-existing"; shift ;;
    --chain-estimate)  CHAIN_ESTIMATE=true;    shift ;;
    --throttle-data)   THROTTLE_DATA="$2";     shift 2 ;;
    --throttle-infer)  THROTTLE_INFER="$2";    shift 2 ;;
    --throttle-align)  THROTTLE_ALIGN="$2";    shift 2 ;;
    --lipid-tokens)    LIPID_TOKENS="$2";      shift 2 ;;
    --dry-run)         DRY_RUN=true;           shift ;;
    -*)
      echo "ERROR: Unknown option: $1" >&2
      echo "Usage: bash submit_af3.sh INPUT_FILE [OPTIONS]" >&2
      exit 1 ;;
    *)
      [[ -z "$INPUT_FILE" ]] && INPUT_FILE="$1" || \
        { echo "ERROR: Unexpected argument: $1" >&2; exit 1; }
      shift ;;
  esac
done

# Validate stage name
VALID_STAGES=(make_json data inference align summarize)
stage_valid=false
for s in "${VALID_STAGES[@]}"; do
  [[ "$s" == "$FROM_STAGE" ]] && stage_valid=true
done
if [[ "$stage_valid" == "false" ]]; then
  echo "ERROR: --from must be one of: ${VALID_STAGES[*]}" >&2; exit 1
fi

if [[ -z "$INPUT_FILE" && "$FROM_STAGE" == "make_json" ]]; then
  echo "ERROR: INPUT_FILE required when starting from make_json stage." >&2
  echo "Usage: bash submit_af3.sh INPUT_FILE [OPTIONS]" >&2
  exit 1
fi

# Normalize and validate --mode if provided
if [[ -n "$MODE" ]]; then
  case "$MODE" in
    regions) MODE="region" ;;  # accept plural alias
  esac
  case "$MODE" in
    fasta|region|accession|hits) ;;
    *) echo "ERROR: --mode '$MODE' is not valid. Use: fasta region accession hits" >&2; exit 1 ;;
  esac
fi

if [[ "$CHAIN_ESTIMATE" == "true" && -z "$INPUT_FILE" ]]; then
  echo "ERROR: --chain-estimate requires INPUT_FILE to estimate array sizes." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# Build the %N throttle suffix for --array.
# If throttle is 0 or empty, omit the suffix entirely (no limit).
# ---------------------------------------------------------------------------
throttle_suffix() {
  local t="$1"
  if [[ "$t" -gt 0 ]]; then
    echo "%${t}"
  else
    echo ""
  fi
}

# ---------------------------------------------------------------------------
# Helper: count files — always returns a clean integer
# ---------------------------------------------------------------------------
count_files() {
  local n
  n=$(find "$@" 2>/dev/null | wc -l | tr -d '[:space:]')
  echo "${n:-0}"
}

# ---------------------------------------------------------------------------
# Helper: count sequences from INPUT_FILE
#   FASTA  (.fa/.fasta or first char '>'): count '>' header lines
#   TSV/accession list / SummaryHits CSV : count non-blank, non-comment,
#                                          non-header lines
# ---------------------------------------------------------------------------
estimate_input_count() {
  local f="$1"
  local n
  if grep -qm1 '^>' "$f" 2>/dev/null; then
    # FASTA: count header lines
    n=$(grep -c '^>' "$f" 2>/dev/null | tr -d '[:space:]' || echo 0)
  else
    # Accession list / SummaryHits CSV
    n=$(grep -cEv '^\s*(#|$)|TMscore|accession' "$f" 2>/dev/null \
        | tr -d '[:space:]' || echo 0)
  fi
  echo "${n:-0}"
}

# ---------------------------------------------------------------------------
# Helper: submit a sbatch job; sets LAST_JOB_ID
# ---------------------------------------------------------------------------
LAST_JOB_ID=""

submit() {
  local label="$1"; shift
  if [[ "$DRY_RUN" == "true" ]]; then
    echo "[DRY-RUN] ${label}: sbatch $*"
    LAST_JOB_ID="DRY_RUN_${label}"
  else
    local out
    out=$(sbatch "$@")
    LAST_JOB_ID=$(echo "$out" | grep -oE '[0-9]+$' || true)
    echo "Submitted ${label}: job ${LAST_JOB_ID}  (${out})"
  fi
}

# ---------------------------------------------------------------------------
# Stage index
# ---------------------------------------------------------------------------
stage_index() {
  case "$1" in
    make_json)  echo 1 ;;
    data)       echo 2 ;;
    inference)  echo 3 ;;
    align)      echo 4 ;;
    summarize)  echo 5 ;;
  esac
}
FROM_IDX=$(stage_index "$FROM_STAGE")

mkdir -p logs

echo "===== submit_af3.sh ====="
echo "Input file      : ${INPUT_FILE:-<not used>}"
echo "From stage      : ${FROM_STAGE} (index ${FROM_IDX})"
echo "Chain range     : ${CHAIN_MIN}–${CHAIN_MAX}"
echo "Lipid tokens    : ${LIPID_TOKENS} (0=disabled)"
echo "Throttle data   : ${THROTTLE_DATA} (0=unlimited)"
echo "Throttle infer  : ${THROTTLE_INFER} (0=unlimited)"
echo "Throttle align  : ${THROTTLE_ALIGN} (0=unlimited)"
echo "Chain estimate  : ${CHAIN_ESTIMATE}"
echo "Dry run         : ${DRY_RUN}"
echo "========================="

# ---------------------------------------------------------------------------
# STAGE 1: make_json
# ---------------------------------------------------------------------------
MKJSON_JID=""
if [[ "$FROM_IDX" -le 1 ]]; then
  MJSON_ARGS=("$INPUT_FILE")
  [[ -n "$MODE"               ]] && MJSON_ARGS+=(--mode "$MODE")
  [[ -n "$SOURCE_TAG"         ]] && MJSON_ARGS+=(--source-tag "$SOURCE_TAG")
  [[ -n "$SKIP_EXISTING_FLAG" ]] && MJSON_ARGS+=("$SKIP_EXISTING_FLAG")
  MJSON_ARGS+=(--min-tm-avg "$MIN_TM_AVG")

  submit "make_json" \
    "${SCRIPT_DIR}/make_json_unified.sbatch" \
    "${MJSON_ARGS[@]}"
  MKJSON_JID="$LAST_JOB_ID"
fi

# ---------------------------------------------------------------------------
# Determine array sizes for stages 2-4
#
# KEY RULE: each stage uses its own output directory count ONLY when
# resuming FROM that stage (files already exist).  When the stage is
# being submitted fresh (predecessor not yet run), derive from the
# predecessor count.  This prevents stale files from old runs inflating
# array sizes.
#
#   FROM_IDX=2 (data)      → count AF3_Inputs/  (the inputs to data)
#   FROM_IDX=3 (inference) → count AF3_Data/    (the inputs to inference)
#   FROM_IDX=4 (align)     → count AF3_Inference/ (the inputs to align)
# ---------------------------------------------------------------------------
N_OLIGOMERS=$(( CHAIN_MAX - CHAIN_MIN + 1 ))
SLURM_ARRAY_MAX=2000   # Rivanna hard limit on array size

# ---------------------------------------------------------------------------
# Stage 2 sizing: always count AF3_Inputs/*_1mer.json
# (these are the actual inputs to the data stage regardless of FROM_IDX)
# ---------------------------------------------------------------------------
N_1MER=0; N_1MER_SOURCE="none"

# Priority 1: count input file directly when provided and starting from
# make_json (FROM_IDX == 1). AF3_Inputs/ may have stale JSONs from a prior
# run — do not trust disk counts until make_json has actually run.
if [[ "$FROM_IDX" -eq 1 && -n "$INPUT_FILE" && -f "$INPUT_FILE" ]]; then
  N_1MER=$(estimate_input_count "$INPUT_FILE")
  [[ "$N_1MER" -gt 0 ]] && N_1MER_SOURCE="input_file"
fi

# Priority 2: disk count — only reliable when resuming from data or later
if [[ "$N_1MER" -eq 0 && "$FROM_IDX" -ge 2 ]]; then
  N_1MER=$(count_files AF3_Inputs -maxdepth 1 -name '*_1mer.json')
  [[ "$N_1MER" -gt 0 ]] && N_1MER_SOURCE="disk"
fi

# Priority 3: downstream job dependency (make_json just submitted, no file yet)
if [[ "$N_1MER" -eq 0 && -n "$MKJSON_JID" ]]; then
  N_1MER_SOURCE="job_dep"
fi

# Priority 4: estimate from input file (fallback when --chain-estimate set)
if [[ "$N_1MER" -eq 0 && "$CHAIN_ESTIMATE" == "true" && -n "$INPUT_FILE" ]]; then
  N_1MER=$(estimate_input_count "$INPUT_FILE")
  [[ "$N_1MER" -gt 0 ]] && N_1MER_SOURCE="estimate"
fi

# ---------------------------------------------------------------------------
# Stage 3 sizing: use AF3_Data/ count ONLY when resuming FROM inference
# (i.e. data stage already ran).  Otherwise derive from N_1MER.
# ---------------------------------------------------------------------------
N_DATA=0; N_DATA_SOURCE="none"
if [[ "$FROM_IDX" -ge 3 ]]; then
  # Resuming from inference or later — data files should already exist
  N_DATA=$(count_files AF3_Data -maxdepth 1 -name '*_data.json')
  [[ "$N_DATA" -gt 0 ]] && N_DATA_SOURCE="disk"
fi
if [[ "$N_DATA" -eq 0 && "$N_1MER" -gt 0 ]]; then
  N_DATA=$(( N_1MER * N_OLIGOMERS ))
  N_DATA_SOURCE="${N_1MER_SOURCE}_derived"
fi

# ---------------------------------------------------------------------------
# Stage 4 sizing: one array task per protein subdirectory in AF3_Inference/.
# Each task aligns all model CIFs for that protein against all references.
# When resuming from align, count subdirs on disk; otherwise derive from N_DATA.
# ---------------------------------------------------------------------------
N_CIF=0; N_CIF_SOURCE="none"
if [[ "$FROM_IDX" -ge 4 ]]; then
  # Resuming from align — subdirectories should already exist
  N_CIF=$(find AF3_Inference -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l)
  [[ "$N_CIF" -gt 0 ]] && N_CIF_SOURCE="disk"
fi
if [[ "$N_CIF" -eq 0 && "$N_DATA" -gt 0 ]]; then
  # One subdir per input sequence (same count as data jobs)
  N_CIF="$N_DATA"
  N_CIF_SOURCE="${N_DATA_SOURCE}_derived"
fi

# ---------------------------------------------------------------------------
# Guard: warn if any array would exceed Rivanna's limit
# ---------------------------------------------------------------------------
for stage_name in "data:${N_1MER}" "inference:${N_DATA}" "align:${N_CIF}"; do
  sname="${stage_name%%:*}"; sval="${stage_name##*:}"
  if [[ "$sval" -gt "$SLURM_ARRAY_MAX" ]]; then
    echo "WARNING: ${sname} array size ${sval} exceeds SLURM_ARRAY_MAX=${SLURM_ARRAY_MAX}." >&2
    echo "         Rivanna may reject this submission." >&2
    echo "         Consider processing in batches or contact RC for a limit increase." >&2
  fi
done

echo "Array size estimates:"
echo "  Stage 2 (data)      : ${N_1MER}  [${N_1MER_SOURCE}]"
echo "  Stage 3 (inference) : ${N_DATA}  [${N_DATA_SOURCE}]"
echo "  Stage 4 (align)     : ${N_CIF}   [${N_CIF_SOURCE}] (jobs = proteins, not CIFs)"

# ---------------------------------------------------------------------------
# STAGE 2: data pipeline
# ---------------------------------------------------------------------------
DATA_JID=""
if [[ "$FROM_IDX" -le 2 ]]; then
  if [[ "$N_1MER" -eq 0 && "$N_1MER_SOURCE" == "job_dep" ]]; then
    echo ""
    echo "NOTE: make_json was just submitted (job ${MKJSON_JID})."
    echo "      Array size for data stage is unknown until make_json finishes."
    echo "      After make_json completes, run:"
    echo "        bash ${BASH_SOURCE[0]} --from data --chain-min ${CHAIN_MIN} --chain-max ${CHAIN_MAX}"
    echo "      Skipping stages 2-5 for now."
    FROM_IDX=99
  elif [[ "$N_1MER" -eq 0 ]]; then
    echo "ERROR: No *_1mer.json files found in AF3_Inputs/ and no predecessor job." >&2
    echo "       Run make_json_unified.sbatch first, then re-run with --from data." >&2
    exit 1
  else
    TSUF=$(throttle_suffix "$THROTTLE_DATA")
    DATA_SCRIPT_TMP="$(mktemp --suffix=.sbatch)"
    sed \
      -e "s/^CHAIN_MIN=.*/CHAIN_MIN=${CHAIN_MIN}/" \
      -e "s/^CHAIN_MAX=.*/CHAIN_MAX=${CHAIN_MAX}/" \
      -e "s/^LIPID_TOKENS=.*/LIPID_TOKENS=${LIPID_TOKENS}/" \
      -e "s/#SBATCH --array=.*/#SBATCH --array=1-${N_1MER}${TSUF}/" \
      "${SCRIPT_DIR}/AF3_data.sbatch" > "$DATA_SCRIPT_TMP"

    DEP_DATA=""
    if [[ -n "$MKJSON_JID" && "$MKJSON_JID" != DRY_RUN* ]]; then
      DEP_DATA="--dependency=afterok:${MKJSON_JID}"
    fi

    submit "data" ${DEP_DATA:+"$DEP_DATA"} "$DATA_SCRIPT_TMP"
    DATA_JID="$LAST_JOB_ID"
    [[ "$DRY_RUN" == "false" ]] && rm -f "$DATA_SCRIPT_TMP"
  fi
fi

# ---------------------------------------------------------------------------
# STAGE 3: inference
# ---------------------------------------------------------------------------
INFER_JID=""
if [[ "$FROM_IDX" -le 3 ]]; then
  if [[ "$N_DATA" -eq 0 ]]; then
    echo "ERROR: Cannot determine inference array size." >&2
    echo "       After data stage completes, re-run with --from inference." >&2
    exit 1
  fi

  TSUF=$(throttle_suffix "$THROTTLE_INFER")
  INFER_SCRIPT_TMP="$(mktemp --suffix=.sbatch)"
  sed \
    -e "s/#SBATCH --array=.*/#SBATCH --array=1-${N_DATA}${TSUF}/" \
    "${SCRIPT_DIR}/AF3_inference.sbatch" > "$INFER_SCRIPT_TMP"

  DEP_INFER=""
  if [[ -n "$DATA_JID" && "$DATA_JID" != DRY_RUN* ]]; then
    DEP_INFER="--dependency=afterok:${DATA_JID}"
  fi

  submit "inference" ${DEP_INFER:+"$DEP_INFER"} "$INFER_SCRIPT_TMP"
  INFER_JID="$LAST_JOB_ID"
  [[ "$DRY_RUN" == "false" ]] && rm -f "$INFER_SCRIPT_TMP"
fi

# ---------------------------------------------------------------------------
# STAGE 4: US-align
# ---------------------------------------------------------------------------
ALIGN_JID=""
if [[ "$FROM_IDX" -le 4 ]]; then
  if [[ "$N_CIF" -eq 0 ]]; then
    echo "ERROR: Cannot determine align array size." >&2
    echo "       After inference stage completes, re-run with --from align." >&2
    exit 1
  fi

  TSUF=$(throttle_suffix "$THROTTLE_ALIGN")
  ALIGN_SCRIPT_TMP="$(mktemp --suffix=.sbatch)"
  sed \
    -e "s/#SBATCH --array=.*/#SBATCH --array=1-${N_CIF}${TSUF}/" \
    "${SCRIPT_DIR}/US_align.sbatch" > "$ALIGN_SCRIPT_TMP"

  DEP_ALIGN=""
  if [[ -n "$INFER_JID" && "$INFER_JID" != DRY_RUN* ]]; then
    # afterany (not afterok) so align runs even if some inference tasks failed.
    # Missing CIFs are logged to Align_Out/skipped.log by US_align.sbatch.
    DEP_ALIGN="--dependency=afterany:${INFER_JID}"
  fi

  submit "align" ${DEP_ALIGN:+"$DEP_ALIGN"} "$ALIGN_SCRIPT_TMP"
  ALIGN_JID="$LAST_JOB_ID"
  [[ "$DRY_RUN" == "false" ]] && rm -f "$ALIGN_SCRIPT_TMP"
fi

# ---------------------------------------------------------------------------
# STAGE 5: summarize
# ---------------------------------------------------------------------------
if [[ "$FROM_IDX" -le 5 ]]; then
  DEP_SUM=""
  if [[ -n "$ALIGN_JID" && "$ALIGN_JID" != DRY_RUN* ]]; then
    DEP_SUM="--dependency=afterok:${ALIGN_JID}"
  fi

  submit "summarize" ${DEP_SUM:+"$DEP_SUM"} "${SCRIPT_DIR}/summarize.sbatch"
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "===== submit_af3.sh complete ====="
echo "Input file  : ${INPUT_FILE:-<not used>}"
echo "From stage  : ${FROM_STAGE}"
echo "Chain range : ${CHAIN_MIN}–${CHAIN_MAX}"
echo "Dry run     : ${DRY_RUN}"
echo "=================================="
