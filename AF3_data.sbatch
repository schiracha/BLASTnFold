#!/bin/bash
#SBATCH -A cavorigins
#SBATCH -p standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH --job-name=AF3_data
#SBATCH --mail-user=rts3dd@virginia.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=logs/af3data_%A_%a.out
#SBATCH --error=logs/af3data_%A_%a.err
#SBATCH --array=1-1%8   # set upper bound to $(ls AF3_Inputs/*_1mer.json | wc -l) before submitting

# =============================================================================
# AF3_data.sbatch
# =============================================================================
# PURPOSE
#   Runs the AF3 data pipeline (MSA + template search) for one 1-mer JSON per
#   array task, then fans out the resulting *_data.json into one copy per
#   oligomer size in [CHAIN_MIN, CHAIN_MAX].  All fan-out copies land flat in
#   AF3_Data/ so AF3_inference.sbatch can find them with a single glob.
#
# INPUT
#   AF3_Inputs/*_1mer.json   (produced by make_json_unified.sbatch)
#
# OUTPUT
#   AF3_Data/<basename>_<N>mer_data.json   for N in [CHAIN_MIN, CHAIN_MAX]
#   AF3_Data/data_pipeline.log             append-only provenance log
#
# USAGE
#   # 1. Count inputs
#   N=$(ls AF3_Inputs/*_1mer.json 2>/dev/null | wc -l)
#   # 2. Edit --array above to 1-$N%8, then:
#   sbatch AF3_data.sbatch
#
# NOTE: The --array upper limit must be set to the number of 1-mer JSONs
#       BEFORE submission (SLURM does not allow dynamic array sizes).
#       Use submit_af3.sbatch for automatic sizing.
# =============================================================================

# ===== CONFIGURATION =====
CHAIN_MIN=5
CHAIN_MAX=15
INPUT_DIR_NAME="AF3_Inputs"    # relative to working directory
OUTPUT_DIR_NAME="AF3_Data"     # flat output directory

# Lipid environment for membrane protein modelling.
# AF3 tokenises ligands at 1 token per heavy atom.  The three lipids used are:
#   POPC (CCD: PLW) — 52 tokens/molecule  — zwitterionic PC head, common in eukaryotes
#   POPE (CCD: EPE) — 49 tokens/molecule  — zwitterionic PE head, common in all kingdoms
#   Cholesterol (CCD: CLR) — 28 tokens/molecule — sterol, modulates membrane fluidity
# One "set" = 1 POPC + 1 POPE + 1 CLR = 129 tokens.
# LIPID_TOKENS controls the total token budget for lipids.
# The number of sets is floor(LIPID_TOKENS / 129); fractional sets are discarded.
#   LIPID_TOKENS=0   → no lipids added (disable entirely)
#   LIPID_TOKENS=129 → 1 set  (129 tokens)
#   LIPID_TOKENS=258 → 2 sets (258 tokens)  ← default, ~250 tokens as requested
#   LIPID_TOKENS=387 → 3 sets (387 tokens)
# Set to 0 to model proteins without a lipid environment.
# Override at submission time with: --lipid-tokens N  (via submit_af3.sh)
LIPID_TOKENS=258
# =========================

set -euo pipefail
trap 'echo "ERROR at line $LINENO (exit $?)" >&2' ERR

module purge
module load gcc/11.4.0 alphafold/3.0.0

# Prevent user-installed alphafold3 packages in ~/.local from shadowing the
# module.  See AF3_inference.sbatch for full explanation.
export PYTHONNOUSERSITE=1

: "${EBROOTALPHAFOLD:?EBROOTALPHAFOLD not set}"
: "${ALPHAFOLD_DATA_PATH:?ALPHAFOLD_DATA_PATH not set}"

WORKDIR=$(pwd -P)
INPUT_DIR="${WORKDIR}/${INPUT_DIR_NAME}"
OUTPUT_DIR="${WORKDIR}/${OUTPUT_DIR_NAME}"
mkdir -p "$OUTPUT_DIR" logs

PIPELINE_LOG="${OUTPUT_DIR}/data_pipeline.log"

# ---------------------------------------------------------------------------
# Build deterministic, sorted input list
# ---------------------------------------------------------------------------
readarray -t INPUT_FILES < <(find "$INPUT_DIR" -maxdepth 1 -type f -name '*_1mer.json' | sort)
TOTAL=${#INPUT_FILES[@]}

if (( TOTAL == 0 )); then
  echo "ERROR: No *_1mer.json files found in ${INPUT_DIR}" >&2
  exit 1
fi

INDEX=$(( SLURM_ARRAY_TASK_ID - 1 ))
if (( INDEX < 0 || INDEX >= TOTAL )); then
  echo "ERROR: Array task ID ${SLURM_ARRAY_TASK_ID} out of range (1..${TOTAL})" >&2
  exit 1
fi

INPUT_FILE="${INPUT_FILES[$INDEX]}"
# Preserve original case for filename; AF3 will lowercase internally
BASENAME="$(basename "$INPUT_FILE" .json)"         # e.g. WP_123_100-200_1mer
BASE_NAME_NO_MER="${BASENAME%_1mer}"               # e.g. WP_123_100-200

echo "===== AF3 data pipeline ====="
echo "Task        : ${SLURM_ARRAY_TASK_ID} / ${TOTAL}"
echo "Input file  : ${INPUT_FILE}"
echo "Base name   : ${BASE_NAME_NO_MER}"
echo "Chain range : ${CHAIN_MIN}..${CHAIN_MAX}"
echo "Node        : $(hostname)"
echo "Started     : $(date -u '+%Y-%m-%dT%H:%M:%SZ')"

# ---------------------------------------------------------------------------
# Skip if all fan-out files already exist (idempotent re-runs)
# ---------------------------------------------------------------------------
all_exist=true
for N in $(seq "$CHAIN_MIN" "$CHAIN_MAX"); do
  [[ $N -eq 1 ]] && continue
  if [[ ! -f "${OUTPUT_DIR}/${BASE_NAME_NO_MER}_${N}mer_data.json" ]]; then
    all_exist=false; break
  fi
done

if [[ "$all_exist" == "true" ]]; then
  echo "All fan-out files already present — skipping data pipeline (use --skip-existing=false to force)."
  exit 0
fi

# ---------------------------------------------------------------------------
# Run AF3 data pipeline (MSA + template search)
# ---------------------------------------------------------------------------
SCRATCH_AF3="${TMPDIR:-/tmp}/af3data_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$SCRATCH_AF3"

python "${EBROOTALPHAFOLD}/app/alphafold/run_alphafold.py" \
    --norun_inference \
    --db_dir="${ALPHAFOLD_DATA_PATH}" \
    --jackhmmer_n_cpu=16 \
    --nhmmer_n_cpu=16 \
    --jax_compilation_cache_dir="/scratch/${USER}/jaxcache" \
    --model_dir="/scratch/rts3dd/AF3" \
    --json_path="${INPUT_FILE}" \
    --output_dir="${SCRATCH_AF3}"

# ---------------------------------------------------------------------------
# Locate the generated *_data.json
# AF3 lowercases the "name" field and uses it as the subdirectory name.
# We search for the file rather than hardcoding the path.
# ---------------------------------------------------------------------------
DATA_FILE=""
while IFS= read -r -d '' f; do
  DATA_FILE="$f"
  break
done < <(find "$SCRATCH_AF3" -type f -name '*_data.json' -print0 | sort -z)

if [[ -z "$DATA_FILE" || ! -f "$DATA_FILE" ]]; then
  echo "ERROR: AF3 data pipeline produced no *_data.json under ${SCRATCH_AF3}" >&2
  exit 1
fi

echo "Data pipeline complete. Found: ${DATA_FILE}"

# ---------------------------------------------------------------------------
# Fan out: one *_data.json per oligomer size, written flat to OUTPUT_DIR
# ---------------------------------------------------------------------------
python3 - "$DATA_FILE" "$BASE_NAME_NO_MER" "$CHAIN_MIN" "$CHAIN_MAX" "$OUTPUT_DIR" "$LIPID_TOKENS" <<'PYEOF'
import json, sys, os

data_file    = sys.argv[1]
base_no_mer  = sys.argv[2]
chain_min    = int(sys.argv[3])
chain_max    = int(sys.argv[4])
output_dir   = sys.argv[5]
lipid_tokens = int(sys.argv[6]) if len(sys.argv) > 6 else 0

TOKEN_LIMIT  = 5000   # AF3 hard limit for total residues + ligand atoms

# ---------------------------------------------------------------------------
# Lipid definitions
# CCD codes present in AF3's bundled CCD database (confirmed in PDB CCD).
# Token cost = number of heavy atoms per molecule.
#   PLW = POPC  (1-palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine)   52 heavy atoms
#   EPE = POPE  (1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine) 49 heavy atoms
#   CLR = Cholesterol                                                  28 heavy atoms
# One "set" = 1 POPC + 1 POPE + 1 CLR = 129 tokens.
# ---------------------------------------------------------------------------
LIPID_SET = [
    {"ccd": "PLW", "name": "POPC", "tokens": 52},
    {"ccd": "EPE", "name": "POPE", "tokens": 49},
    {"ccd": "CLR", "name": "Cholesterol", "tokens": 28},
]
TOKENS_PER_SET = sum(l["tokens"] for l in LIPID_SET)  # 129

# Compute how many complete sets fit in the lipid token budget
n_lipid_sets = lipid_tokens // TOKENS_PER_SET if lipid_tokens > 0 else 0
actual_lipid_tokens = n_lipid_sets * TOKENS_PER_SET

if lipid_tokens > 0 and n_lipid_sets == 0:
    print(f"WARNING: LIPID_TOKENS={lipid_tokens} is less than one set "
          f"({TOKENS_PER_SET} tokens). No lipids will be added.",
          file=sys.stderr)

if n_lipid_sets > 0:
    print(f"Lipid environment: {n_lipid_sets} set(s) × "
          f"(1 POPC + 1 POPE + 1 CLR) = {actual_lipid_tokens} tokens")
else:
    print("Lipid environment: disabled (LIPID_TOKENS=0 or too small)")

with open(data_file) as fh:
    template = json.load(fh)

chains_alphabet = [chr(65 + i) for i in range(26)]

# Determine sequence length from the first protein entry in the template
seq_len = 0
for seq_entry in template.get("sequences", []):
    if "protein" in seq_entry:
        seq_len = len(seq_entry["protein"].get("sequence", ""))
        break

if seq_len == 0:
    print("ERROR: Could not determine sequence length from data file.", file=sys.stderr)
    sys.exit(1)

# ---------------------------------------------------------------------------
# Token budget check: protein chains + lipids must not exceed TOKEN_LIMIT.
# Cap chain_max so that  n * seq_len + actual_lipid_tokens <= TOKEN_LIMIT.
# If even a monomer + lipids exceeds the limit, reduce lipids first, then
# fall back to no lipids, then error.
# ---------------------------------------------------------------------------
original_chain_max = chain_max

# Adjust lipid tokens down if they alone would eat too much budget
available_for_protein = TOKEN_LIMIT - actual_lipid_tokens
if available_for_protein < seq_len:
    # Can't even fit a monomer — try with no lipids
    print(f"WARNING: {actual_lipid_tokens} lipid tokens leave only "
          f"{available_for_protein} tokens for protein (seq_len={seq_len}). "
          f"Disabling lipids for {base_no_mer}.", file=sys.stderr)
    n_lipid_sets = 0
    actual_lipid_tokens = 0
    available_for_protein = TOKEN_LIMIT

max_feasible = available_for_protein // seq_len

if max_feasible < 1:
    print(f"ERROR: Sequence length {seq_len} aa exceeds token limit even as monomer. "
          f"Skipping {base_no_mer}.", file=sys.stderr)
    sys.exit(1)

if chain_max > max_feasible:
    print(f"WARNING: chain_max={chain_max} would produce "
          f"{chain_max * seq_len + actual_lipid_tokens} total tokens "
          f"(limit {TOKEN_LIMIT}). Capping chain_max at {max_feasible} "
          f"({max_feasible * seq_len + actual_lipid_tokens} tokens) for {base_no_mer}.")
    chain_max = max_feasible

if chain_min > chain_max:
    print(f"WARNING: chain_min={chain_min} > capped chain_max={chain_max} for "
          f"{base_no_mer}. No fan-out files written.")
    sys.exit(0)

print(f"Sequence length : {seq_len} aa")
print(f"Token budget    : chain_min={chain_min}, chain_max={chain_max} "
      f"(requested max was {original_chain_max})")

# ---------------------------------------------------------------------------
# Build lipid sequence entries for injection.
# Each lipid molecule gets its own unique chain ID, assigned after the
# protein chains.  Chain IDs start at A for proteins; lipids continue
# from the next available letter.
# For n protein chains + n_lipid_sets * 3 lipid molecules, we need
# n + n_lipid_sets * 3 ≤ 26 chain IDs total.
# ---------------------------------------------------------------------------
def build_lipid_entries(n_protein_chains, n_sets):
    """Return a list of AF3 JSON sequence entries for the lipid molecules."""
    if n_sets == 0:
        return []
    entries = []
    chain_idx = n_protein_chains  # start after protein chains
    for _ in range(n_sets):
        for lipid in LIPID_SET:
            if chain_idx >= 26:
                print(f"WARNING: Ran out of chain IDs at chain_idx={chain_idx}. "
                      f"Truncating lipid entries.", file=sys.stderr)
                return entries
            chain_id = chains_alphabet[chain_idx]
            entries.append({
                "ligand": {
                    "id": [chain_id],
                    "ccdCodes": [lipid["ccd"]]
                }
            })
            chain_idx += 1
    return entries

# ---------------------------------------------------------------------------
# Fan-out loop
# ---------------------------------------------------------------------------
for n in range(chain_min, chain_max + 1):
    if n == 1:
        continue
    if n > 26:
        print(f"WARNING: n={n} exceeds 26-chain alphabet, skipping.", file=sys.stderr)
        continue

    # Check that n protein chains + lipid molecules fit in 26-chain limit
    n_lipid_molecules = n_lipid_sets * len(LIPID_SET)
    if n + n_lipid_molecules > 26:
        # Reduce lipid sets to fit
        max_sets = (26 - n) // len(LIPID_SET)
        if max_sets < n_lipid_sets:
            print(f"  NOTE: {n}-mer + {n_lipid_sets} lipid sets = "
                  f"{n + n_lipid_molecules} chains > 26 limit. "
                  f"Reducing to {max_sets} lipid sets for this oligomer.")
        n_sets_this = max_sets
    else:
        n_sets_this = n_lipid_sets

    data = json.loads(json.dumps(template))   # deep copy

    new_name = f"{base_no_mer}_{n}mer"
    data["name"] = new_name

    # Update protein chain IDs
    for seq_entry in data.get("sequences", []):
        if "protein" in seq_entry:
            seq_entry["protein"]["id"] = chains_alphabet[:n]

    # Inject lipid entries after protein sequences
    lipid_entries = build_lipid_entries(n, n_sets_this)
    if lipid_entries:
        # Remove any existing lipid entries from template (shouldn't exist
        # in 1-mer inputs, but defensive)
        data["sequences"] = [e for e in data["sequences"] if "ligand" not in e]
        data["sequences"].extend(lipid_entries)

    total_tokens = n * seq_len + n_sets_this * TOKENS_PER_SET
    out_path = os.path.join(output_dir, f"{new_name}_data.json")
    with open(out_path, "w") as fh:
        json.dump(data, fh, indent=2)
    lipid_str = (f" + {n_sets_this}x(POPC+POPE+CLR)={n_sets_this*TOKENS_PER_SET} lipid tokens"
                 if n_sets_this > 0 else " (no lipids)")
    print(f"  Created: {out_path}  "
          f"({n}×{seq_len}={n*seq_len} protein tokens{lipid_str} = {total_tokens} total)")

print(f"Fan-out complete: {chain_min}-mer to {chain_max}-mer")
PYEOF

# ---------------------------------------------------------------------------
# Append to provenance log
# ---------------------------------------------------------------------------
{
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" \
    "$(basename "$INPUT_FILE")" \
    "$BASE_NAME_NO_MER" \
    "$CHAIN_MIN" \
    "$CHAIN_MAX"
} >> "$PIPELINE_LOG"

# Cleanup scratch
rm -rf "$SCRATCH_AF3"

echo "===== AF3 data pipeline DONE ====="
echo "Outputs in : ${OUTPUT_DIR}/"
echo "Finished   : $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
