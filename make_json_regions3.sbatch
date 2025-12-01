#!/usr/bin/env bash
#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=standard
#SBATCH -J "AF3Regions"
#SBATCH --mail-user=rts3dd@virginia.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --account=cavorigins

set -euo pipefail

# Usage: sbatch make_af3_jsons.sh input.txt
# Input format: ACCESSION:START-END (e.g., WP_251811466.1:209-369)

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 input.txt"
  exit 1
fi

INPUT_FILE="$1"

# ---- Output directory ----
OUTPUT_DIR="$(pwd -P)/AF3_Seqs"
mkdir -p "$OUTPUT_DIR"

missing_log="$OUTPUT_DIR/missing_accessions.log"
found_log="$OUTPUT_DIR/found_accessions.log"
: > "$missing_log"
: > "$found_log"

# ---- BLAST+ setup ----
module purge
module load blast
command -v blastdbcmd >/dev/null 2>&1 || { echo "ERROR: blastdbcmd not found."; exit 2; }

DB_DIR="${DB_DIR:-/scratch/rts3dd/BLAST/reference_db}"
export BLASTDB="${DB_DIR}${BLASTDB:+:$BLASTDB}"
DBS="${DBS:-nr env_nr tsa_nr}"

# ---- Extract triplets: accession (no version), start, end ----
triplets="$(mktemp)"

awk -F '[:\\-]' '{
  acc = $1;
  sub(/\.[0-9]+$/, "", acc);  # remove version suffix
  print acc, $2, $3;
}' "$INPUT_FILE" > "$triplets"

# ---- Sanitize accession for filenames ----
sanitize() {
  echo "$1" | tr -c 'A-Za-z0-9.-' '_' | sed 's/_$//'
}

# ---- Process one triplet ----
process_one() {
  local acc="$1" s="$2" e="$3"

  if (( s > e )); then
    local tmp="$s"; s="$e"; e="$tmp"
  fi

  local seq
  seq="$(blastdbcmd -db "$DBS" -dbtype prot -entry "$acc" -range "${s}-${e}" -outfmt "%f" -line_length 1000000 2>/dev/null || true)"

  if [[ -z "${seq:-}" ]]; then
    echo "Warning: $acc:$s-$e not found or out-of-range." >> "$missing_log"
    return 0
  fi

  local seq_only
  seq_only="$(printf "%s" "$seq" | sed -n '2,$p' | tr -d '\n\r')"

  if [[ -z "$seq_only" ]]; then
    echo "Warning: $acc:$s-$e yielded empty sequence." >> "$missing_log"
    return 0
  fi

  local chains=(A B C D E F G H I J K)
  local chain_json
  chain_json="$(printf '"%s", ' "${chains[@]}")"
  chain_json="[${chain_json%, }]"

  local safe_acc; safe_acc="$(sanitize "$acc")"
  local json_file="$OUTPUT_DIR/${safe_acc}.json"

  cat > "$json_file" <<EOF
{
  "name": "${acc}:${s}-${e}",
  "modelSeeds": [1, 2],
  "sequences": [
    {
      "protein": {
        "id": $chain_json,
        "sequence": "$seq_only",
        "unpairedMsa": "",
        "pairedMsa": "",
        "templates": []
      }
    }
  ],
  "dialect": "alphafold3",
  "version": 1
}
EOF

  echo "$acc:$s-$e (len=${#seq_only}, chains=11) -> $(basename "$json_file")" >> "$found_log"
}

export -f process_one
export OUTPUT_DIR BLASTDB DBS missing_log found_log
export -f sanitize

# ---- Parallel execution ----
xargs -n 3 -P "${SLURM_CPUS_PER_TASK:-8}" bash -c 'process_one "$0" "$1" "$2"' < "$triplets"
rm -f "$triplets"

# ---- Summary ----
echo "Done. JSON files saved in $OUTPUT_DIR/"
echo "Found list: $found_log ($(wc -l < "$found_log") entries)"
echo "Missing list: $missing_log ($(wc -l < "$missing_log") entries)"
