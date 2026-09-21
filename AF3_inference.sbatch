#!/bin/bash
#SBATCH -A cavorigins
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH -C "a6000|a100"
#SBATCH -t 15:00:00
#SBATCH --job-name=AF3_infer
#SBATCH --mail-user=rts3dd@virginia.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=logs/af3infer_%A_%a.out
#SBATCH --error=logs/af3infer_%A_%a.err
#SBATCH --array=1-1%4   # set upper bound to $(ls AF3_Data/*_data.json | wc -l) before submitting

# =============================================================================
# AF3_inference.sbatch
# =============================================================================
# PURPOSE
#   Runs AF3 structure prediction (inference only) for one *_data.json per
#   array task.  Default: 2 seeds × 5 diffusion samples = 10 models per run.
#
# OUTPUT NAMING CONVENTION
#   Single-seed runs  (modelSeeds has 1 entry):
#     AF3_Inference/<basename>/<basename>_0.cif  ..  <basename>_4.cif
#
#   Multi-seed runs   (modelSeeds has >1 entry):
#     AF3_Inference/<basename>/<basename>_seed1_sample0.cif
#     AF3_Inference/<basename>/<basename>_seed1_sample1.cif  ...
#     AF3_Inference/<basename>/<basename>_seed2_sample0.cif  ...
#
#   The script auto-detects which convention applies from the AF3 output
#   directory structure and uses the appropriate naming to avoid collisions.
#
# INPUT
#   AF3_Data/*_data.json       (produced by AF3_data.sbatch, flat directory)
#
# OUTPUT
#   AF3_Inference/<basename>/<basename>_*.cif
#   AF3_Inference/inference.log    append-only provenance log
#
# USAGE
#   N=$(ls AF3_Data/*_data.json 2>/dev/null | wc -l)
#   # Edit --array to 1-$N%4, then submit. Or use submit_af3.sh.
# =============================================================================

set -euo pipefail
trap 'echo "ERROR at line $LINENO (exit $?)" >&2' ERR

module purge
module load gcc/11.4.0 alphafold/3.0.0

# Prevent Python from picking up user-installed packages in ~/.local.
# Without this, a previously installed alphafold3==3.0.1 wheel in
# ~/.local/lib/python3.11/site-packages/ takes precedence over the module
# and causes FileNotFoundError on pickle files whose paths changed between
# versions.  Setting PYTHONNOUSERSITE=1 forces Python to ignore ~/.local
# entirely for this job, using only the module's site-packages.
export PYTHONNOUSERSITE=1

: "${EBROOTALPHAFOLD:?EBROOTALPHAFOLD not set (AF3 install root)}"
: "${ALPHAFOLD_DATA_PATH:?ALPHAFOLD_DATA_PATH not set (AF3 databases path)}"

WORKDIR=$(pwd -P)
INPUT_DIR="${WORKDIR}/AF3_Data"
OUTPUT_DIR="${WORKDIR}/AF3_Inference"
mkdir -p "$OUTPUT_DIR" logs

INFER_LOG="${OUTPUT_DIR}/inference.log"

# ---------------------------------------------------------------------------
# Build deterministic sorted input list
# ---------------------------------------------------------------------------
readarray -t INPUT_FILES < <(find "$INPUT_DIR" -maxdepth 1 -type f -name '*_data.json' | sort)
TOTAL=${#INPUT_FILES[@]}

if (( TOTAL == 0 )); then
  echo "ERROR: No *_data.json files found in ${INPUT_DIR}" >&2
  exit 1
fi

INDEX=$(( SLURM_ARRAY_TASK_ID - 1 ))
if (( INDEX < 0 || INDEX >= TOTAL )); then
  echo "ERROR: Array task ID ${SLURM_ARRAY_TASK_ID} out of range (1..${TOTAL})" >&2
  exit 1
fi

INPUT_FILE="${INPUT_FILES[$INDEX]}"
BASENAME="$(basename "$INPUT_FILE" _data.json)"

# ---------------------------------------------------------------------------
# Detect number of seeds from the JSON modelSeeds array
# Used for naming convention selection and memory tuning
# ---------------------------------------------------------------------------
N_SEEDS=$(python3 -c "
import json, sys
d = json.load(open('${INPUT_FILE}'))
print(len(d.get('modelSeeds', [1])))
" 2>/dev/null || echo 1)

# ---------------------------------------------------------------------------
# Extract oligomer size for memory tuning
# ---------------------------------------------------------------------------
if [[ "$BASENAME" =~ _([0-9]+)mer$ ]]; then
  N_CHAINS="${BASH_REMATCH[1]}"
else
  echo "WARNING: Could not parse oligomer size from '${BASENAME}'; defaulting to 11" >&2
  N_CHAINS=11
fi

# Tune JAX memory fraction by oligomer size
if   (( N_CHAINS <= 8  )); then export XLA_CLIENT_MEM_FRACTION=2.5
elif (( N_CHAINS >= 12 )); then export XLA_CLIENT_MEM_FRACTION=3.5
else                             export XLA_CLIENT_MEM_FRACTION=3.2
fi

export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TF_FORCE_UNIFIED_MEMORY=true
export OMP_NUM_THREADS=1

echo "===== AF3 inference ====="
echo "Task        : ${SLURM_ARRAY_TASK_ID} / ${TOTAL}"
echo "Input       : ${INPUT_FILE}"
echo "Basename    : ${BASENAME}   (${N_CHAINS}-mer, ${N_SEEDS} seed(s))"
echo "Node        : $(hostname)"
echo "GPU         : ${CUDA_VISIBLE_DEVICES:-<unset>}"
echo "Mem fraction: ${XLA_CLIENT_MEM_FRACTION}"
echo "Started     : $(date -u '+%Y-%m-%dT%H:%M:%SZ')"

# ---------------------------------------------------------------------------
# Idempotency: skip if output already has at least one CIF
# ---------------------------------------------------------------------------
DEST="${OUTPUT_DIR}/${BASENAME}"
if ls "${DEST}/${BASENAME}_"*.cif 2>/dev/null | grep -q .; then
  echo "Output CIFs already present in ${DEST}/ — skipping inference."
  exit 0
fi

# ---------------------------------------------------------------------------
# Run inference on scratch, then rsync to permanent storage
# ---------------------------------------------------------------------------
SCRATCH="${TMPDIR:-/tmp}/af3infer_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$SCRATCH"

export JAX_COMPILATION_CACHE_DIR="/scratch/${USER}/jaxcache/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$JAX_COMPILATION_CACHE_DIR"

python "${EBROOTALPHAFOLD}/app/alphafold/run_alphafold.py" \
    --norun_data_pipeline \
    --db_dir="${ALPHAFOLD_DATA_PATH}" \
    --jax_compilation_cache_dir="${JAX_COMPILATION_CACHE_DIR}" \
    --model_dir="/scratch/rts3dd/AF3" \
    --json_path="${INPUT_FILE}" \
    --output_dir="${SCRATCH}"

# Locate AF3 output subdirectory (AF3 lowercases the name field)
AF3_OUT_SUBDIR=""
while IFS= read -r -d '' d; do
  AF3_OUT_SUBDIR="$d"; break
done < <(find "$SCRATCH" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z)

if [[ -z "$AF3_OUT_SUBDIR" ]]; then
  echo "ERROR: AF3 produced no output subdirectory under ${SCRATCH}" >&2
  exit 1
fi

echo "AF3 output subdir: ${AF3_OUT_SUBDIR}"
mkdir -p "$DEST"
rsync -a --no-owner --no-group "${AF3_OUT_SUBDIR}/" "${DEST}/"

# ---------------------------------------------------------------------------
# Flatten and rename model files.
#
# AF3 output layout (seed-S_sample-M subdirectories, AF3 3.0.x standard):
#   seed-1_sample-0/model.cif
#   seed-1_sample-1/model.cif  ...
#   seed-2_sample-0/model.cif  (if 2 seeds)  ...
#
# SINGLE-SEED naming  → <basename>_<sample>.cif
#   e.g. basename_0.cif, basename_1.cif ... basename_4.cif
#
# MULTI-SEED naming   → <basename>_seed<S>_sample<M>.cif
#   e.g. basename_seed1_sample0.cif, basename_seed2_sample0.cif ...
#   This avoids collisions when sample indices repeat across seeds.
#
# Convention is chosen from N_SEEDS detected above.
# ---------------------------------------------------------------------------
n_renamed=0

for seed_dir in "${DEST}"/seed-*_sample-*/; do
  [[ -d "$seed_dir" ]] || continue
  [[ -f "${seed_dir}model.cif" ]] || continue
  dir_name="$(basename "$seed_dir")"

  if [[ "$dir_name" =~ seed-([0-9]+)_sample-([0-9]+) ]]; then
    seed_num="${BASH_REMATCH[1]}"
    sample_num="${BASH_REMATCH[2]}"

    if [[ "$N_SEEDS" -gt 1 ]]; then
      # Multi-seed: encode both seed and sample to avoid collisions
      target="${DEST}/${BASENAME}_seed${seed_num}_sample${sample_num}.cif"
    else
      # Single-seed: simple sample index only
      target="${DEST}/${BASENAME}_${sample_num}.cif"
    fi

    if [[ ! -f "$target" ]]; then
      mv "${seed_dir}model.cif" "$target"
      (( n_renamed++ )) || true
    fi
    rm -rf "$seed_dir"
  else
    echo "WARNING: Could not parse seed/sample from '${dir_name}'" >&2
  fi
done

# Fallback layouts (some AF3 minor versions)
if [[ -f "${DEST}/${BASENAME}_model.cif" && ! -f "${DEST}/${BASENAME}_0.cif" ]]; then
  mv "${DEST}/${BASENAME}_model.cif" "${DEST}/${BASENAME}_0.cif"
  (( n_renamed++ )) || true
fi
if [[ -f "${DEST}/model.cif" && ! -f "${DEST}/${BASENAME}_0.cif" ]]; then
  mv "${DEST}/model.cif" "${DEST}/${BASENAME}_0.cif"
  (( n_renamed++ )) || true
fi

if (( n_renamed == 0 )); then
  echo "WARNING: No model CIF files found to rename in ${DEST}/" >&2
fi

# ---------------------------------------------------------------------------
# Provenance log
# ---------------------------------------------------------------------------
{
  printf '%s\t%s\t%s\t%s\t%s\t%d\n' \
    "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" \
    "$BASENAME" \
    "$N_CHAINS" \
    "$N_SEEDS" \
    "$(hostname)" \
    "$n_renamed"
} >> "$INFER_LOG"

rm -rf "$SCRATCH" "$JAX_COMPILATION_CACHE_DIR"

echo "Renamed ${n_renamed} model CIF(s) → ${DEST}/"
echo "Naming: $( [[ "$N_SEEDS" -gt 1 ]] && echo "multi-seed (seed+sample)" || echo "single-seed (sample only)" )"
echo "===== AF3 inference DONE ====="
echo "Finished : $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
