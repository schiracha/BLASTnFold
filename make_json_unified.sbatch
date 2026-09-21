#!/usr/bin/env bash
#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=standard
#SBATCH -J "AF3_mkjson"
#SBATCH --mail-user=rts3dd@virginia.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --account=cavorigins
#SBATCH --output=logs/mkjson_%j.out
#SBATCH --error=logs/mkjson_%j.err

# =============================================================================
# make_json_unified.sbatch
# =============================================================================
# PURPOSE
#   Unified front-end for the AF3 pipeline.  Accepts four input flavours and
#   emits one 1-mer JSON per entry into AF3_Inputs/, plus a provenance manifest
#   (AF3_Inputs/manifest.tsv).
#
# INPUT FLAVOURS  (auto-detected; override with --mode)
#   fasta       FASTA file (.fasta/.fa/.afa) — sequences already in hand,
#               no BLAST fetch required.  Accession parsed from header line:
#               the first word after '>'.  Designed for output from
#               ClusterHits.sbatch (cluster_out/clustered_sequences.fasta).
#
#   region      ACCESSION:START-END lines  (e.g. WP_251811466.1:209-369)
#               Also accepts tab-separated  ACCESSION  START  END
#               Fetches the specified region from the BLAST databases.
#
#   accession   Plain accession list, one per line — fetches full-length
#               sequences from the BLAST databases.
#
#   hits        SummaryHits CSV produced by summarize.sh — fetches full-length
#               sequences for hits passing TM-score thresholds.
#
# USAGE
#   sbatch make_json_unified.sbatch  INPUT_FILE  [OPTIONS]
#
# OPTIONS
#   --mode  fasta|region|accession|hits   force input flavour (default: auto-detect)
#   --outdir  DIR                   output directory (default: AF3_Inputs)
#   --min-tm-avg  FLOAT             hits mode: min TMscore_avgL  (default 0.2)
#   --min-tm-refl FLOAT             hits mode: min TMscore_refL  (default 0.0)
#   --source-tag  TAG               free-text provenance label written to manifest
#   --skip-existing                 skip JSONs that already exist (resume mode)
#
# TYPICAL WORKFLOWS
#   # From ClusterHits output (recommended for large inputs)
#   sbatch make_json_unified.sbatch cluster_out/clustered_sequences.fasta \
#       --source-tag ACC000_clustered
#
#   # From BLAST region list
#   sbatch make_json_unified.sbatch HITS_WP123_01_cd.out --mode region
#
#   # From plain accession list
#   sbatch make_json_unified.sbatch ACCESSIONS_master_01.txt --mode accession
#
#   # From SummaryHits (iterative re-modelling)
#   sbatch make_json_unified.sbatch Run01_SummaryHits.csv --min-tm-avg 0.3
#
# OUTPUT
#   AF3_Inputs/<safe_acc>_1mer.json          — one per FASTA entry or accession
#   AF3_Inputs/<safe_acc>_START-END_1mer.json — one per region entry
#   AF3_Inputs/manifest.tsv                  — provenance table (appended)
#   AF3_Inputs/found_accessions.log
#   AF3_Inputs/missing_accessions.log        — only populated for DB-fetch modes
#
# MANIFEST COLUMNS
#   json_basename  accession  region  seq_len  source_tag  input_file  timestamp
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# 0. Defaults
# ---------------------------------------------------------------------------
MODE=""
OUTPUT_DIR=""
MIN_TM_AVG=0.2
MIN_TM_REFL=0.0
SOURCE_TAG=""
SKIP_EXISTING=false
INPUT_FILE=""

# ---------------------------------------------------------------------------
# 1. Argument parsing
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)          MODE="$2";         shift 2 ;;
    --outdir)        OUTPUT_DIR="$2";   shift 2 ;;
    --min-tm-avg)    MIN_TM_AVG="$2";  shift 2 ;;
    --min-tm-refl)   MIN_TM_REFL="$2"; shift 2 ;;
    --source-tag)    SOURCE_TAG="$2";   shift 2 ;;
    --skip-existing) SKIP_EXISTING=true; shift ;;
    -*)
      echo "ERROR: Unknown option: $1"
      echo "Usage: sbatch $0 INPUT_FILE [--mode fasta|region|accession|hits]"
      exit 1 ;;
    *)
      if [[ -z "$INPUT_FILE" ]]; then INPUT_FILE="$1"
      else echo "ERROR: Unexpected positional argument: $1"; exit 1; fi
      shift ;;
  esac
done

if [[ -z "$INPUT_FILE" ]]; then
  echo "ERROR: No input file specified."
  echo "Usage: sbatch $0 INPUT_FILE [OPTIONS]"
  exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "ERROR: Input file not found: $INPUT_FILE"
  exit 1
fi

OUTPUT_DIR="${OUTPUT_DIR:-$(pwd -P)/AF3_Inputs}"
SOURCE_TAG="${SOURCE_TAG:-$(basename "$INPUT_FILE")}"
mkdir -p "$OUTPUT_DIR" logs

MANIFEST="$OUTPUT_DIR/manifest.tsv"
MISSING_LOG="$OUTPUT_DIR/missing_accessions.log"
FOUND_LOG="$OUTPUT_DIR/found_accessions.log"

# Manifest is append-only (provenance log across multiple runs).
# Missing/found logs are per-run diagnostics — always overwrite so the
# summary count at the end of this run is accurate.
[[ -f "$MANIFEST" ]] || printf 'json_basename\taccession\tregion\tseq_len\tsource_tag\tinput_file\ttimestamp\n' > "$MANIFEST"
: > "$MISSING_LOG"
: > "$FOUND_LOG"

# ---------------------------------------------------------------------------
# 2. Auto-detect input flavour
# ---------------------------------------------------------------------------
detect_mode() {
  local f="$1"
  # FASTA: first non-blank line starts with >
  local first_line
  first_line=$(grep -m1 -Ev '^\s*$' "$f" || true)
  if [[ "$first_line" == ">"* ]]; then
    echo "fasta"; return
  fi
  # SummaryHits CSV: header contains TMscore
  if head -1 "$f" | grep -qi 'TMscore\|tms'; then
    echo "hits"; return
  fi
  # Tab-separated 3-column with integer coords → region (_cd.tsv)
  if echo "$first_line" | awk -F'\t' 'NF==3 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/{exit 0} {exit 1}'; then
    echo "region"; return
  fi
  # 12-column consolidated TSV: header row = "kingdom\ttaxid\tsskingdom\taccession..."
  # OR data row: col 7 and 8 are integers and col 4 is an accession-like string
  if echo "$first_line" | awk -F'\t' '
      NF >= 8 && $1=="kingdom" && $4=="accession" && $7=="sstart" { exit 0 }
      NF >= 8 && $7~/^[0-9]+$/ && $8~/^[0-9]+$/ && $4~/^[A-Za-z0-9_]/ { exit 0 }
      { exit 1 }'; then
    echo "region"; return
  fi
  # Colon+dash range pattern → region (colon-style)
  if echo "$first_line" | grep -qE '[A-Za-z0-9_]+(\.[0-9]+)?:[0-9]+-[0-9]+'; then
    echo "region"; return
  fi
  echo "accession"
}

if [[ -z "$MODE" ]]; then
  MODE=$(detect_mode "$INPUT_FILE")
  echo "Auto-detected input mode: $MODE"
else
  # Normalize aliases (e.g. "regions" → "region")
  case "$MODE" in
    regions) MODE="region" ;;
  esac
  # Validate
  case "$MODE" in
    fasta|region|accession|hits) ;;
    *) echo "ERROR: Unknown --mode '$MODE'. Valid: fasta region accession hits"; exit 1 ;;
  esac
  echo "Using forced input mode: $MODE"
fi

echo "Input file   : $INPUT_FILE"
echo "Output dir   : $OUTPUT_DIR"
echo "Source tag   : $SOURCE_TAG"
echo "Mode         : $MODE"
echo "Skip existing: $SKIP_EXISTING"

# ---------------------------------------------------------------------------
# 3. BLAST+ setup — only needed for non-fasta modes
# ---------------------------------------------------------------------------
if [[ "$MODE" != "fasta" ]]; then
  module purge
  module load blast
  command -v blastdbcmd >/dev/null 2>&1 || \
    { echo "ERROR: blastdbcmd not found after loading blast module."; exit 2; }
  DB_DIR="${DB_DIR:-/scratch/rts3dd/BLAST/reference_db}"
  export BLASTDB="${DB_DIR}${BLASTDB:+:$BLASTDB}"
  DBS="${DBS:-nr env_nr tsa_nr}"
else
  # Export dummies so exported functions don't reference unbound variables
  export BLASTDB="" DBS=""
fi

# ---------------------------------------------------------------------------
# 4. Shared helpers
# ---------------------------------------------------------------------------
sanitize() {
  printf '%s' "$1" | tr -c 'A-Za-z0-9.-' '_' | sed 's/_$//'
}
export -f sanitize

TIMESTAMP=$(date -u '+%Y-%m-%dT%H:%M:%SZ')
export OUTPUT_DIR DBS BLASTDB MISSING_LOG FOUND_LOG MANIFEST \
       SOURCE_TAG INPUT_FILE TIMESTAMP SKIP_EXISTING

WORKER_TMP_DIR="$(mktemp -d)"
export WORKER_TMP_DIR

write_json() {
  local acc="$1" s="$2" e="$3" seq_only="$4"

  local seq_len="${#seq_only}"
  local safe_acc; safe_acc="$(sanitize "$acc")"

  local region_tag json_name
  if [[ -n "$s" && -n "$e" ]]; then
    region_tag="${s}-${e}"
    json_name="${safe_acc}_${region_tag}_1mer"
  else
    region_tag="full"
    json_name="${safe_acc}_1mer"
  fi

  local json_file="$OUTPUT_DIR/${json_name}.json"

  if [[ "$SKIP_EXISTING" == "true" && -f "$json_file" ]]; then
    echo "SKIP (exists): $json_file" >&2
    return 0
  fi

  local internal_name="${acc}"
  [[ -n "$s" && -n "$e" ]] && internal_name="${acc}_${region_tag}"
  internal_name="${internal_name}_1mer"

  cat > "$json_file" <<ENDJSON
{
  "name": "${internal_name}",
  "modelSeeds": [1, 2],
  "sequences": [
    {
      "protein": {
        "id": ["A"],
        "sequence": "${seq_only}",
        "unpairedMsa": "",
        "pairedMsa": "",
        "templates": []
      }
    }
  ],
  "dialect": "alphafold3",
  "version": 1
}
ENDJSON

  local wid="$$"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(basename "$json_file")" \
    "$acc" \
    "$region_tag" \
    "$seq_len" \
    "$SOURCE_TAG" \
    "$(basename "$INPUT_FILE")" \
    "$TIMESTAMP" \
    >> "${WORKER_TMP_DIR}/manifest_${wid}.tsv"

  echo "${acc}|region=${region_tag}|len=${seq_len}|json=$(basename "$json_file")" \
    >> "${WORKER_TMP_DIR}/found_${wid}.log"
}
export -f write_json

# ---------------------------------------------------------------------------
# 5. Per-mode processing functions
# ---------------------------------------------------------------------------

# --- 5a. FASTA mode: parse sequences directly, no DB fetch ---
process_fasta_entry() {
  # Called with: HEADER_LINE (the full >... line)
  # Sequence is read from a temp file written by the FASTA splitter below.
  # Args: acc seq_only
  local acc="$1"
  local seq_only="$2"

  if [[ -z "$seq_only" ]]; then
    echo "EMPTY_SEQ: ${acc}" >> "${WORKER_TMP_DIR}/missing_$$.log"
    return 0
  fi

  # Sanitize any gap characters introduced by alignment (clustered_sequences
  # may be aligned FASTA with '-' gaps from MAFFT; strip them)
  seq_only="${seq_only//-/}"
  seq_only="${seq_only//./.}"   # keep dots as-is (rare in protein FASTA)

  write_json "$acc" "" "" "$seq_only"
}
export -f process_fasta_entry

# --- 5b. Region mode ---
process_region() {
  local acc="$1" s="$2" e="$3"
  local acc_nover="${acc%%.*}"
  if (( s > e )); then local tmp=$s; s=$e; e=$tmp; fi

  local seq
  seq="$(blastdbcmd -db "$DBS" -dbtype prot -entry "$acc_nover" \
           -range "${s}-${e}" -outfmt "%f" -line_length 1000000 2>/dev/null || true)"

  if [[ -z "${seq:-}" ]]; then
    echo "MISSING: ${acc}:${s}-${e}" >> "${WORKER_TMP_DIR}/missing_$$.log"
    return 0
  fi

  local seq_only
  seq_only="$(printf '%s' "$seq" | sed -n '2,$p' | tr -d '\n\r')"
  [[ -z "$seq_only" ]] && { echo "EMPTY_SEQ: ${acc}:${s}-${e}" >> "${WORKER_TMP_DIR}/missing_$$.log"; return 0; }

  write_json "$acc_nover" "$s" "$e" "$seq_only"
}
export -f process_region

# --- 5c. Accession mode (also used by hits mode) ---
process_accession() {
  local acc="$1"
  local acc_nover="${acc%%.*}"

  local seq
  seq="$(blastdbcmd -db "$DBS" -dbtype prot -entry "$acc_nover" \
           -outfmt "%f" -line_length 1000000 2>/dev/null || true)"

  if [[ -z "${seq:-}" ]]; then
    echo "MISSING: ${acc}" >> "${WORKER_TMP_DIR}/missing_$$.log"
    return 0
  fi

  local seq_only
  seq_only="$(printf '%s' "$seq" | sed -n '2,$p' | tr -d '\n\r')"
  [[ -z "$seq_only" ]] && { echo "EMPTY_SEQ: ${acc}" >> "${WORKER_TMP_DIR}/missing_$$.log"; return 0; }

  write_json "$acc_nover" "" "" "$seq_only"
}
export -f process_accession

# ---------------------------------------------------------------------------
# 6. Build work list and dispatch
# ---------------------------------------------------------------------------
CPUS="${SLURM_CPUS_PER_TASK:-8}"
WORK="$(mktemp)"

case "$MODE" in

  fasta)
    # Parse FASTA into (accession, sequence) pairs using Python — handles
    # multi-line sequences cleanly and is much safer than awk for this.
    # Writes a tab-separated work file: ACC\tSEQUENCE
    # Accession = first whitespace-delimited token after '>' (version kept)
    echo "Parsing FASTA entries..."
    python3 - "$INPUT_FILE" "$WORK" << 'PYEOF'
import sys, re

fasta_file = sys.argv[1]
work_file  = sys.argv[2]

entries = []
acc = None
seq_parts = []

with open(fasta_file) as fh:
    for line in fh:
        line = line.rstrip('\n\r')
        if line.startswith('>'):
            if acc is not None:
                seq = ''.join(seq_parts).replace('-', '').replace(' ', '')
                if seq:
                    entries.append((acc, seq))
            # First token after '>' is the accession; strip version suffix
            raw = line[1:].split()[0] if line[1:].strip() else "unknown"
            acc = re.sub(r'\.\d+$', '', raw)
            seq_parts = []
        elif acc is not None:
            seq_parts.append(line.strip())

# flush last entry
if acc is not None:
    seq = ''.join(seq_parts).replace('-', '').replace(' ', '')
    if seq:
        entries.append((acc, seq))

with open(work_file, 'w') as fh:
    for a, s in entries:
        fh.write(f"{a}\t{s}\n")

print(f"Parsed {len(entries)} FASTA entries", flush=True)
PYEOF

    total=$(wc -l < "$WORK")
    echo "FASTA entries to process: $total"

    # Process in parallel: each line is ACC\tSEQ
    # -d '\n' ensures xargs splits only on newlines, not tabs — without it,
    # the tab between ACC and SEQ would be treated as a delimiter, passing
    # the accession and sequence as separate arguments (seq always empty).
    export -f process_fasta_entry
    xargs -d '\n' -n 1 -P "$CPUS" bash -c '
      IFS=$'"'"'\t'"'"' read -r acc seq <<< "$0"
      process_fasta_entry "$acc" "$seq"
    ' < "$WORK"
    ;;

  region)
    # Three supported sub-formats — auto-detected per line:
    #   12-col consolidated TSV  (kingdom taxid sskingdom ACC evalue bitscore SSTART SEND ...)
    #     → acc=$4  sstart=$7  send=$8
    #   3-col _cd.tsv            (accession sstart send)
    #     → acc=$1  sstart=$2  send=$3
    #   Colon-style              (ACC:START-END)
    #     → split on ':' then '-'
    # Header lines are skipped (column 1 is a word, column 2 is not a number).
    awk -F'\t' '
      /^\s*(#|$)/ { next }
      NR == 1 && $2 !~ /^[0-9]+$/ { next }   # skip header row
      {
        acc = ""; s = ""; e = ""
        if (NF == 1 && index($1, ":") > 0) {
          # Colon-style: ACC:START-END
          n = split($1, a, ":")
          if (n == 2) { m = split(a[2], b, "-"); if (m == 2) { acc=a[1]; s=b[1]; e=b[2] } }
        } else if (NF >= 8 && $7 ~ /^[0-9]+$/ && $8 ~ /^[0-9]+$/ && $4 !~ /^[0-9]/) {
          # 12-column consolidated format (col 4=acc, 7=sstart, 8=send)
          acc=$4; s=$7; e=$8
        } else if (NF >= 3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/) {
          # 3-column _cd format (col 1=acc, 2=sstart, 3=send)
          acc=$1; s=$2; e=$3
        }
        if (acc != "" && s ~ /^[0-9]+$/ && e ~ /^[0-9]+$/) {
          sub(/\.[0-9]+$/, "", acc)
          print acc, s, e
        }
      }
    ' "$INPUT_FILE" > "$WORK"

    total=$(wc -l < "$WORK")
    echo "Region entries to process: $total"
    xargs -n 3 -P "$CPUS" bash -c 'process_region "$0" "$1" "$2"' < "$WORK"
    ;;

  accession)
    grep -Ev '^\s*(#|$)' "$INPUT_FILE" | awk '{print $1}' > "$WORK"

    total=$(wc -l < "$WORK")
    echo "Accession entries to process: $total"
    xargs -n 1 -P "$CPUS" bash -c 'process_accession "$0"' < "$WORK"
    ;;

  hits)
    awk -F',' -v minavg="$MIN_TM_AVG" -v minrefl="$MIN_TM_REFL" '
      NR == 1 { next }
      NR == 2 { fmt = (NF >= 7) ? "new" : "old" }
      {
        acc = $1
        if (fmt == "new") { refl = $6+0; avg = $7+0 }
        else               { refl = $4+0; avg = $5+0 }
        if (minavg  > 0 && avg  < minavg)  next
        if (minrefl > 0 && refl < minrefl) next
        print acc
      }
    ' "$INPUT_FILE" | sort -u > "$WORK"

    total=$(wc -l < "$WORK")
    echo "Hits entries passing thresholds (min_avg=${MIN_TM_AVG}, min_refl=${MIN_TM_REFL}): $total"

    if (( total == 0 )); then
      echo "No entries passed thresholds. Exiting."
      rm -f "$WORK"; rm -rf "$WORKER_TMP_DIR"; exit 0
    fi
    xargs -n 1 -P "$CPUS" bash -c 'process_accession "$0"' < "$WORK"
    ;;

  *)
    echo "ERROR: Unknown mode '$MODE'"
    rm -f "$WORK"; rm -rf "$WORKER_TMP_DIR"; exit 1
    ;;
esac

rm -f "$WORK"

# ---------------------------------------------------------------------------
# 7. Merge per-worker temp files into shared logs/manifest
# ---------------------------------------------------------------------------
find "$WORKER_TMP_DIR" -name 'manifest_*.tsv' | sort | xargs -r cat >> "$MANIFEST"
find "$WORKER_TMP_DIR" -name 'found_*.log'    | sort | xargs -r cat >> "$FOUND_LOG"
find "$WORKER_TMP_DIR" -name 'missing_*.log'  | sort | xargs -r cat >> "$MISSING_LOG"
rm -rf "$WORKER_TMP_DIR"

# ---------------------------------------------------------------------------
# 8. Summary
# ---------------------------------------------------------------------------
n_found=$(wc -l < "$FOUND_LOG")
n_missing=$(wc -l < "$MISSING_LOG")
echo ""
echo "===== make_json_unified COMPLETE ====="
echo "Mode         : $MODE"
echo "Source       : $INPUT_FILE"
echo "Output dir   : $OUTPUT_DIR"
echo "JSONs written: $n_found"
echo "Missing      : $n_missing"
echo "Manifest     : $MANIFEST"
echo "======================================"
