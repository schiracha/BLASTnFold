#!/usr/bin/env bash
#SBATCH -A cavorigins
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH -t 1:00:00
#SBATCH --job-name=US_align
#SBATCH --mail-user=rts3dd@virginia.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=logs/usalign_%A_%a.out
#SBATCH --error=logs/usalign_%A_%a.err
#SBATCH --array=1-1%50   # set upper bound to number of protein subdirs; see USAGE below

# =============================================================================
# US_align.sbatch
# =============================================================================
# PURPOSE
#   Runs USalign for one protein (AF3_Inference subdirectory) per array task.
#   All model CIFs for that protein (typically 10: 2 seeds × 5 samples) are
#   aligned against all six reference structures within the single job.
#   USalign is lightweight — one CPU and ~1 hour is ample per protein.
#
#   Output per task:
#     Align_Out/<REF>_<basename>_stats.out   — alignment metrics (parsed by summarize.sh)
#     Align_Out/<REF>_<basename>_stats.err   — stderr
#     Align_Out/<REF>_<basename>.txt         — rotation matrix
#     Align_Out/<REF>_<basename>.cif         — superimposed structure
#
# INPUT
#   AF3_Inference/*/   subdirectories, each containing *_[0-9]*.cif model files
#   Reference CIF files in working directory (auto-downloaded if absent)
#
# USAGE
#   N=$(find AF3_Inference -mindepth 1 -maxdepth 1 -type d | wc -l)
#   sbatch US_align.sbatch   # after editing --array to 1-$N%50
#   # Or use submit_af3.sh for automatic sizing and dependency chaining.
#
# REFERENCES
#   Positive controls : 7SC0, 9DN1, 8BJM
#   Negative controls : 7D60, 7RSL
#   Neutral (monomer) : AF-Q03135  (AlphaFold2 human Cav-1; no -mm 1)
#   Neutral (complex) : AF-Complex (optional; place AF-Complex.cif in workdir)
#
# NOTE ON AF-Q03135 ALIGNMENT
#   AF-Q03135 is a monomer. When aligning a multi-chain query against it,
#   USalign without -mm 1 compares the best-matching single chain of the query
#   to the full monomer. This is the intended "neutral monomer" comparison.
#   -mm 1 is intentionally OMITTED for this reference only.
#
# NOTE ON REFERENCE DOWNLOADS
#   Task 1 downloads any missing PDB/AlphaFold CIFs and writes a sentinel
#   (refs/.downloads_complete). All other tasks wait via flock poll.
# =============================================================================

set -euo pipefail
trap 'echo "ERROR at line $LINENO (exit $?)" >&2' ERR

export OMP_NUM_THREADS=1

WORKDIR=$(pwd -P)
INPUT_DIR="${WORKDIR}/AF3_Inference"
OUTPUT_DIR="${WORKDIR}/Align_Out"
REFS_DIR="${WORKDIR}/refs"
mkdir -p "$OUTPUT_DIR" "$REFS_DIR" logs

# Shared log for CIFs that could not be aligned because the input file is
# missing (i.e. the inference task failed).  Written with flock so concurrent
# array tasks do not interleave lines.
SKIP_LOG="${OUTPUT_DIR}/skipped.log"

# ---------------------------------------------------------------------------
# Reference structure definitions
#   REF_FILE  : local CIF filename (relative to WORKDIR)
#   REF_MONO  : "yes" → omit -mm 1 (monomer ref); "no" → include -mm 1 (complex ref)
# ---------------------------------------------------------------------------
declare -A REF_FILE REF_MONO
REF_FILE=(
  [7SC0]="7SC0.cif"
  [9DN1]="9DN1.cif"
  [8BJM]="8BJM.cif"
  [7D60]="7D60.cif"
  [7RSL]="7RSL.cif"
  [AF-Q03135]="AF-Q03135.cif"
  [AFComp]="AF-Complex.cif"
)
REF_MONO=(
  [7SC0]="no"
  [9DN1]="no"
  [8BJM]="no"
  [7D60]="no"
  [7RSL]="no"
  [AF-Q03135]="yes"
  [AFComp]="no"
)
# Ordered: positive controls first, then negative, then neutral
REF_ORDER=(7SC0 9DN1 8BJM 7D60 7RSL AF-Q03135 AFComp)

# ---------------------------------------------------------------------------
# Download missing reference structures (task 1 only).
# flock sentinel prevents race conditions across parallel array tasks.
# ---------------------------------------------------------------------------
SENTINEL="${REFS_DIR}/.downloads_complete"
LOCKFILE="${REFS_DIR}/.download_lock"

if [[ "$SLURM_ARRAY_TASK_ID" -eq 1 ]]; then
  (
    flock -x 200
    if [[ ! -f "$SENTINEL" ]]; then
      declare -A REF_URL=(
        [7SC0.cif]="https://files.rcsb.org/download/7SC0.cif"
        [9DN1.cif]="https://files.rcsb.org/download/9DN1.cif"
        [8BJM.cif]="https://files.rcsb.org/download/8BJM.cif"
        [7D60.cif]="https://files.rcsb.org/download/7D60.cif"
        [7RSL.cif]="https://files.rcsb.org/download/7RSL.cif"
        [AF-Q03135.cif]="https://alphafold.ebi.ac.uk/files/AF-Q03135-2-F1-model_v6.cif"
      )
      for fname in "${!REF_URL[@]}"; do
        if [[ ! -f "$fname" ]]; then
          echo "Downloading ${fname} ..."
          wget -q -O "$fname" "${REF_URL[$fname]}"
          echo "  Retrieved: ${fname}"
        else
          echo "  Present  : ${fname}"
        fi
      done
      if [[ ! -f "AF-Complex.cif" ]]; then
        echo "WARNING: AF-Complex.cif not found. Place it in ${WORKDIR}/ before running."
        echo "  (AFComp alignments will be skipped for all tasks.)"
      else
        echo "  Present  : AF-Complex.cif"
      fi
      touch "$SENTINEL"
    fi
  ) 200>"$LOCKFILE"
else
  waited=0; max_wait=120
  while [[ ! -f "$SENTINEL" && $waited -lt $max_wait ]]; do
    sleep 2; (( waited += 2 )) || true
  done
  if [[ ! -f "$SENTINEL" ]]; then
    echo "WARNING: Download sentinel not found after ${max_wait}s; proceeding anyway." >&2
  fi
fi

# ---------------------------------------------------------------------------
# Build deterministic sorted list of protein subdirectories.
# Each subdirectory contains all model CIFs for one input sequence.
# ---------------------------------------------------------------------------
readarray -t INPUT_DIRS < <(
  find "$INPUT_DIR" -mindepth 1 -maxdepth 1 -type d | sort
)
TOTAL=${#INPUT_DIRS[@]}

if (( TOTAL == 0 )); then
  echo "ERROR: No protein subdirectories found under ${INPUT_DIR}" >&2
  exit 1
fi

INDEX=$(( SLURM_ARRAY_TASK_ID - 1 ))
if (( INDEX < 0 || INDEX >= TOTAL )); then
  echo "ERROR: Array task ID ${SLURM_ARRAY_TASK_ID} out of range (1..${TOTAL})" >&2
  exit 1
fi

PROTEIN_DIR="${INPUT_DIRS[$INDEX]}"
PROTEIN_NAME="$(basename "$PROTEIN_DIR")"

# Collect all model CIFs for this protein
readarray -t MODEL_FILES < <(
  find "$PROTEIN_DIR" -maxdepth 1 -type f -name '*_[0-9]*.cif' | sort
)

echo "===== USalign ====="
echo "Protein  : ${PROTEIN_NAME}"
echo "Task     : ${SLURM_ARRAY_TASK_ID} / ${TOTAL}"
echo "Models   : ${#MODEL_FILES[@]}"
echo "Node     : $(hostname)"
echo "OMP      : ${OMP_NUM_THREADS}"
echo "Started  : $(date -u '+%Y-%m-%dT%H:%M:%SZ')"

if (( ${#MODEL_FILES[@]} == 0 )); then
  echo "WARNING: No model CIFs found in ${PROTEIN_DIR} — skipping."
  exit 0
fi

# ---------------------------------------------------------------------------
# Outer loop: each model CIF for this protein
# Inner loop: each reference structure
# Missing CIFs (failed inference) are logged to skipped.log and skipped.
# ---------------------------------------------------------------------------
for INPUT_FILE in "${MODEL_FILES[@]}"; do
  BASENAME="$(basename "${INPUT_FILE}" .cif)"

  if [[ ! -f "$INPUT_FILE" ]]; then
    MSG="$(date -u '+%Y-%m-%dT%H:%M:%SZ')\tMISSING_CIF\t${SLURM_ARRAY_TASK_ID}\t${INPUT_FILE}"
    (
      flock -x 201
      printf '%s\n' "$MSG" >> "$SKIP_LOG"
    ) 201>"${SKIP_LOG}.lock"
    echo "  SKIPPED (missing): ${INPUT_FILE}"
    continue
  fi

  echo "  --- Model: ${BASENAME} ---"

  # Emit Basename as a parseable field for summarize.sh
  echo "Basename: ${BASENAME}"

  for REF_KEY in "${REF_ORDER[@]}"; do
    REF_CIF="${REF_FILE[$REF_KEY]}"
    IS_MONO="${REF_MONO[$REF_KEY]}"

    if [[ ! -f "$REF_CIF" ]]; then
      echo "    SKIP ${REF_KEY}: '${REF_CIF}' not present."
      continue
    fi

    STATS_OUT="${OUTPUT_DIR}/${REF_KEY}_${BASENAME}_stats.out"
    if [[ -f "$STATS_OUT" ]]; then
      echo "    SKIP ${REF_KEY}: already done."
      continue
    fi

    USALIGN_EXTRA_FLAGS=()
    if [[ "$IS_MONO" == "no" ]]; then
      USALIGN_EXTRA_FLAGS+=( -mm 1 )
    fi

    USalign \
      "$INPUT_FILE" \
      "$REF_CIF" \
      -mol  prot \
      "${USALIGN_EXTRA_FLAGS[@]+"${USALIGN_EXTRA_FLAGS[@]}"}" \
      -ter  0 \
      -TMscore 0 \
      -m    "${OUTPUT_DIR}/${REF_KEY}_${BASENAME}.txt" \
      -do \
      -split 2 \
      -infmt1 3 \
      -infmt2 3 \
      -a    T \
      -outfmt -1 \
      -o    "${OUTPUT_DIR}/${REF_KEY}_${BASENAME}" \
      > "${STATS_OUT}" \
      2>"${OUTPUT_DIR}/${REF_KEY}_${BASENAME}_stats.err"

    echo "    Done: ${REF_KEY}"
  done
done

echo "===== USalign DONE ====="
echo "Finished : $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
