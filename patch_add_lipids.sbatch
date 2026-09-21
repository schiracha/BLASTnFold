#!/usr/bin/env bash
#SBATCH -A cavorigins
#SBATCH -p standard
#SBATCH --cpus-per-task=4
#SBATCH -t 1:00:00
#SBATCH --job-name=patch_lipids
#SBATCH --mail-user=rts3dd@virginia.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=logs/patch_lipids_%j.out
#SBATCH --error=logs/patch_lipids_%j.err

# =============================================================================
# patch_add_lipids.sbatch
# =============================================================================
# PURPOSE
#   Adds lipid environment entries to all existing *_data.json files in
#   AF3_Data/ (or a specified directory).  Idempotent: files that already
#   contain lipid entries are skipped.
#
# LIPIDS ADDED (equal mix, standard eukaryotic membrane composition)
#   PLW = POPC        52 tokens/molecule
#   EPE = POPE        49 tokens/molecule
#   CLR = Cholesterol 28 tokens/molecule
#   One "set" = 1×POPC + 1×POPE + 1×CLR = 129 tokens
#
# USAGE
#   sbatch patch_add_lipids.sbatch [OPTIONS]
#
# OPTIONS
#   --lipid-tokens N    total lipid token budget (default: 258 = 2 sets)
#                       set to 0 to remove existing lipids instead
#   --data-dir DIR      directory containing *_data.json files
#                       (default: ./AF3_Data)
#   --dry-run           print what would be done without modifying files
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
LIPID_TOKENS=258
DATA_DIR="$(pwd -P)/AF3_Data"
DRY_RUN=false

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --lipid-tokens) LIPID_TOKENS="$2"; shift 2 ;;
    --data-dir)     DATA_DIR="$2";     shift 2 ;;
    --dry-run)      DRY_RUN=true;      shift ;;
    *) echo "ERROR: Unknown option: $1" >&2; exit 1 ;;
  esac
done

if [[ ! -d "$DATA_DIR" ]]; then
  echo "ERROR: Directory not found: $DATA_DIR" >&2
  exit 1
fi

mkdir -p logs

echo "===== patch_add_lipids ====="
echo "Data dir     : $DATA_DIR"
echo "Lipid tokens : $LIPID_TOKENS"
echo "Dry run      : $DRY_RUN"
echo "Started      : $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
echo ""

# ---------------------------------------------------------------------------
# Count input files
# ---------------------------------------------------------------------------
mapfile -t JSON_FILES < <(find "$DATA_DIR" -maxdepth 1 -type f -name '*_data.json' | sort)
TOTAL=${#JSON_FILES[@]}

if (( TOTAL == 0 )); then
  echo "No *_data.json files found in $DATA_DIR"
  exit 0
fi

echo "Found $TOTAL *_data.json files"
echo ""

# ---------------------------------------------------------------------------
# Python patch — processes all files in one interpreter call for efficiency
# ---------------------------------------------------------------------------
python3 - "$DATA_DIR" "$LIPID_TOKENS" "$DRY_RUN" "${JSON_FILES[@]}" << 'PYEOF'
import json, sys, os, shutil
from pathlib import Path

data_dir     = sys.argv[1]
lipid_tokens = int(sys.argv[2])
dry_run      = sys.argv[3].lower() == "true"
json_files   = sys.argv[4:]

# ---------------------------------------------------------------------------
# Lipid definitions
# ---------------------------------------------------------------------------
LIPID_SET = [
    {"ccd": "PLW", "name": "POPC",         "tokens": 52},
    {"ccd": "EPE", "name": "POPE",         "tokens": 49},
    {"ccd": "CLR", "name": "Cholesterol",  "tokens": 28},
]
TOKENS_PER_SET = sum(l["tokens"] for l in LIPID_SET)  # 129
CHAINS = [chr(65 + i) for i in range(26)]

n_sets = lipid_tokens // TOKENS_PER_SET if lipid_tokens > 0 else 0

if lipid_tokens > 0 and n_sets == 0:
    print(f"WARNING: LIPID_TOKENS={lipid_tokens} < one set ({TOKENS_PER_SET}). "
          f"No lipids will be added.")
    sys.exit(0)

if lipid_tokens == 0:
    print(f"LIPID_TOKENS=0: will REMOVE any existing lipid entries.")
else:
    print(f"Adding {n_sets} set(s) × (POPC+POPE+CLR) = "
          f"{n_sets * TOKENS_PER_SET} tokens per file")

print()

n_patched  = 0
n_skipped  = 0
n_removed  = 0
n_error    = 0

for fpath in json_files:
    fname = os.path.basename(fpath)
    try:
        with open(fpath) as fh:
            data = json.load(fh)

        existing_ligands = [e for e in data.get("sequences", []) if "ligand" in e]
        protein_entries  = [e for e in data.get("sequences", []) if "protein" in e]

        # ---- REMOVE mode (lipid_tokens == 0) ----
        if lipid_tokens == 0:
            if not existing_ligands:
                print(f"  SKIP  {fname}  (no ligands present)")
                n_skipped += 1
                continue
            if not dry_run:
                data["sequences"] = protein_entries
                with open(fpath, "w") as fh:
                    json.dump(data, fh, indent=2)
            print(f"  {'WOULD REMOVE' if dry_run else 'REMOVED'} "
                  f"{len(existing_ligands)} ligand entry(s) from {fname}")
            n_removed += 1
            continue

        # ---- ADD mode ----
        # Idempotency: skip if correct lipid entries already present
        existing_ccds = [e["ligand"]["ccdCodes"][0]
                         for e in existing_ligands if "ligand" in e]
        expected_ccds = [l["ccd"] for l in LIPID_SET] * n_sets
        if existing_ccds == expected_ccds:
            print(f"  SKIP  {fname}  (already has correct lipids)")
            n_skipped += 1
            continue

        # Determine number of protein chains from the protein id list
        n_protein_chains = 0
        for e in protein_entries:
            ids = e["protein"].get("id", [])
            n_protein_chains = max(n_protein_chains, len(ids))

        # Adjust sets down if we'd exceed 26 chains
        n_lipid_molecules = n_sets * len(LIPID_SET)
        effective_sets = n_sets
        if n_protein_chains + n_lipid_molecules > 26:
            effective_sets = (26 - n_protein_chains) // len(LIPID_SET)
            if effective_sets == 0:
                print(f"  SKIP  {fname}  "
                      f"(no room for lipids: {n_protein_chains} protein chains "
                      f"already uses {n_protein_chains}/26 chain IDs)")
                n_skipped += 1
                continue
            print(f"  NOTE  {fname}  "
                  f"reducing to {effective_sets} lipid sets to stay within 26 chains")

        # Build lipid entries
        chain_idx = n_protein_chains
        lipid_entries = []
        for _ in range(effective_sets):
            for lipid in LIPID_SET:
                lipid_entries.append({
                    "ligand": {
                        "id": [CHAINS[chain_idx]],
                        "ccdCodes": [lipid["ccd"]]
                    }
                })
                chain_idx += 1

        if not dry_run:
            # Replace any existing ligand entries with the new ones
            data["sequences"] = protein_entries + lipid_entries
            with open(fpath, "w") as fh:
                json.dump(data, fh, indent=2)

        action = "WOULD PATCH" if dry_run else "PATCHED"
        replaced = f" (replaced {len(existing_ligands)} existing)" if existing_ligands else ""
        print(f"  {action} {fname}  "
              f"(+{len(lipid_entries)} lipid entries{replaced})")
        n_patched += 1

    except Exception as e:
        print(f"  ERROR {fname}: {e}", file=sys.stderr)
        n_error += 1

print()
print("===== Summary =====")
if lipid_tokens == 0:
    print(f"  Removed  : {n_removed}")
else:
    print(f"  Patched  : {n_patched}")
print(f"  Skipped  : {n_skipped}  (already correct or no room)")
print(f"  Errors   : {n_error}")
if dry_run:
    print("  (dry run — no files modified)")
PYEOF

echo ""
echo "Finished : $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
