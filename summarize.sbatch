#!/usr/bin/env bash
#SBATCH -A cavorigins
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH -t 0:30:00
#SBATCH --job-name=af3_summarize
#SBATCH --mail-user=rts3dd@virginia.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=logs/summarize_%j.out
#SBATCH --error=logs/summarize_%j.err

# =============================================================================
# summarize.sbatch
# =============================================================================
# PURPOSE
#   Parses all USalign stats files in Align_Out/ and emits:
#     - One per-reference CSV  (e.g. 7SC0_<PARDIR>.csv)
#     - A merged SummaryHits CSV filtered at TMscore_avgL > 0.2
#     - A merged SummaryHits_all CSV containing every alignment result
#
# USAGE
#   sbatch summarize.sbatch                  # Align_Out/ relative to CWD
#   sbatch summarize.sbatch /path/to/AlignOut
#   bash   summarize.sbatch                  # also runnable directly
#   # Or use submit_af3.sh for automatic dependency chaining.
#
# OUTPUT COLUMNS (per-reference CSVs)
#   basename, accession, region, oligomer_size, sample,
#   aligned_len, RMSD, Seq_ID,
#   TMscore_refL, TMscore_L1, TMscore_L2, TMscore_avgL
#
# OUTPUT COLUMNS (SummaryHits CSV)
#   accession, region, oligomer_size, sample,
#   reference, TMscore_refL, TMscore_avgL
#
# ACCESSION PARSING
#   Handles both naming conventions produced by make_json_unified.sbatch:
#     New (seed+sample) : <acc>_<N>mer_seed<S>_sample<M>          → region = full
#     New (region)      : <acc>_<start>-<end>_<N>mer_seed<S>_sample<M>
#     Old (single token): <acc>_<N>mer_<sample>                    → region = full
#     Old (region)      : <acc>_<start>-<end>_<N>mer_<sample>
# =============================================================================

set -euo pipefail

ALIGN_DIR="${1:-$(pwd -P)/Align_Out}"
PARDIR="$(basename "$(pwd -P)")"

if [[ ! -d "$ALIGN_DIR" ]]; then
  echo "ERROR: Align_Out directory not found: ${ALIGN_DIR}" >&2
  exit 1
fi

mkdir -p logs

# ---------------------------------------------------------------------------
# Output file names
# ---------------------------------------------------------------------------
out_7SC0="7SC0_${PARDIR}.csv"
out_9DN1="9DN1_${PARDIR}.csv"
out_8BJM="8BJM_${PARDIR}.csv"
out_7D60="7D60_${PARDIR}.csv"
out_7RSL="7RSL_${PARDIR}.csv"
out_AFQ="AF-Q03135_${PARDIR}.csv"
out_AFComp="AFComp_${PARDIR}.csv"
out_hits="${PARDIR}_SummaryHits.csv"
out_hits_all="${PARDIR}_SummaryHits_all.csv"

HEADER="basename,accession,region,oligomer_size,sample,aligned_len,RMSD,Seq_ID,TMscore_refL,TMscore_L1,TMscore_L2,TMscore_avgL"
for f in "$out_7SC0" "$out_9DN1" "$out_8BJM" "$out_7D60" "$out_7RSL" "$out_AFQ" "$out_AFComp"; do
  echo "$HEADER" > "$f"
done

# ---------------------------------------------------------------------------
# Collect stats files per reference
# ---------------------------------------------------------------------------
mapfile -t CANDS_7SC0    < <(find "$ALIGN_DIR" -maxdepth 1 -type f -name '7SC0_*_stats.out'      2>/dev/null | sort -u)
mapfile -t CANDS_9DN1    < <(find "$ALIGN_DIR" -maxdepth 1 -type f -name '9DN1_*_stats.out'      2>/dev/null | sort -u)
mapfile -t CANDS_8BJM    < <(find "$ALIGN_DIR" -maxdepth 1 -type f -name '8BJM_*_stats.out'      2>/dev/null | sort -u)
mapfile -t CANDS_7D60    < <(find "$ALIGN_DIR" -maxdepth 1 -type f -name '7D60_*_stats.out'      2>/dev/null | sort -u)
mapfile -t CANDS_7RSL    < <(find "$ALIGN_DIR" -maxdepth 1 -type f -name '7RSL_*_stats.out'      2>/dev/null | sort -u)
mapfile -t CANDS_AFQ     < <(find "$ALIGN_DIR" -maxdepth 1 -type f -name 'AF-Q03135_*_stats.out' 2>/dev/null | sort -u)
mapfile -t CANDS_AFComp  < <(find "$ALIGN_DIR" -maxdepth 1 -type f -name 'AFComp_*_stats.out'    2>/dev/null | sort -u)

total_files=$(( ${#CANDS_7SC0[@]} + ${#CANDS_9DN1[@]} + ${#CANDS_8BJM[@]} + \
               ${#CANDS_7D60[@]} + ${#CANDS_7RSL[@]}   + \
               ${#CANDS_AFQ[@]}  + ${#CANDS_AFComp[@]} ))

if (( total_files == 0 )); then
  echo "No stats files found in ${ALIGN_DIR}/" >&2
  exit 1
fi
echo "Found ${total_files} stats files across all references."

# ---------------------------------------------------------------------------
# extract_one  FILE
#   Parses one USalign stats file.
#   Emits a CSV row:
#     basename,accession,region,oligomer_size,sample,
#     aligned_len,RMSD,Seq_ID,
#     TMscore_refL,TMscore_L1,TMscore_L2,TMscore_avgL
# ---------------------------------------------------------------------------
extract_one() {
  local f="$1"

  # --- Derive raw basename ---
  local raw_basename=""

  # First try explicit Basename: header line (written by some wrappers)
  raw_basename=$(grep -m1 -E '^Basename:' "$f" \
    | sed -E 's/^Basename:[[:space:]]+//' | tr -d '\r') || true

  # Fallback: parse "Selected input" line written by US_align.sbatch
  if [[ -z "$raw_basename" ]]; then
    local inpath
    inpath=$(grep -m1 'Selected input' "$f" \
      | sed -E 's/.*Selected input[^:]*:[[:space:]]*//' | tr -d '\r ') || true
    [[ -n "$inpath" ]] && raw_basename="$(basename "$inpath" .cif)"
  fi

  # Final fallback: derive from filename (strip REF_ prefix and _stats.out)
  if [[ -z "$raw_basename" ]]; then
    local stem
    stem="$(basename "$f" _stats.out)"
    stem="${stem#7SC0_}"
    stem="${stem#9DN1_}"
    stem="${stem#8BJM_}"
    stem="${stem#7D60_}"
    stem="${stem#7RSL_}"
    stem="${stem#AF-Q03135_}"
    stem="${stem#AFComp_}"
    raw_basename="$stem"
  fi

  # --- Parse accession / region / oligomer_size / sample from basename ---
  # Handles two AF3 naming conventions:
  #
  #   New (seed+sample tokens, AF3 3.0+):
  #     <acc>_<N>mer_seed<S>_sample<M>
  #     <acc>_<start>-<end>_<N>mer_seed<S>_sample<M>
  #
  #   Old (single sample token):
  #     <acc>_<N>mer_<sample>
  #     <acc>_<start>-<end>_<N>mer_<sample>
  #
  # Strategy: detect trailing token style first, then work backwards to find
  # the <N>mer field, optional region field, and accession prefix.

  local accession="" region="" oligomer_size="" sample=""

  read -r accession region oligomer_size sample < <(
    echo "$raw_basename" | awk -F'_' '
    {
      n = NF

      # Detect new-style seed+sample suffix: last = "sample<M>", prev = "seed<S>"
      if (n >= 2 && $n ~ /^sample[0-9]+$/ && $(n-1) ~ /^seed[0-9]+$/) {
        s = $(n-1) "_" $n      # e.g. "seed1_sample4"
        nmer_idx = n - 2
      } else {
        # Old style: last field is the sample token
        s = $n
        nmer_idx = n - 1
      }

      nmer_field = $nmer_idx
      if (nmer_field ~ /^[0-9]+mer$/) {
        sub(/mer$/, "", nmer_field)
        oligomer = nmer_field

        # Check for region field immediately before the Nmer field
        if (nmer_idx >= 3 && $(nmer_idx-1) ~ /^[0-9]+-[0-9]+$/) {
          reg = $(nmer_idx-1)
          acc = ""
          for (i = 1; i <= nmer_idx-2; i++) acc = (acc == "") ? $i : acc "_" $i
        } else {
          reg = "full"
          acc = ""
          for (i = 1; i <= nmer_idx-1; i++) acc = (acc == "") ? $i : acc "_" $i
        }
      } else {
        # Unexpected pattern — surface the raw basename for manual inspection
        oligomer = "?"; reg = "?"; s = "?"; acc = $0
      }
      print acc, reg, oligomer, s
    }'
  )

  # --- Extract alignment metrics from USalign output ---
  local aligned="" rmsd="" seqid="" tms1="" tms2="" tmsavg="" tms_ref=""

  if grep -q 'Aligned length=' "$f"; then
    local line
    line=$(grep -m1 'Aligned length=' "$f" | tr -d '\r')

    aligned=$(echo "$line" | sed -E 's/.*Aligned length=[[:space:]]*([0-9]+).*/\1/') || true
    rmsd=$(echo "$line"    | sed -E 's/.*RMSD=[[:space:]]*([0-9.]+).*/\1/')          || true
    seqid=$(echo "$line"   | awk '{
        gsub(",", " ")
        for (i=1; i<=NF; i++) {
          if ($i ~ /Seq_ID/) {
            for (j=i+1; j<=NF; j++) {
              if ($j ~ /^[0-9]*\.[0-9]+$/ || $j ~ /^[0-9]+$/) { seq=$j }
            }
          }
        }
        if (seq != "") print seq
      }') || true
  fi

  tms1=$(grep -E 'TM-score=' "$f" | grep -m1 'Structure_1' \
    | sed -E 's/.*TM-score[= ]+([0-9.]+).*/\1/' | tr -d '\r') || true
  tms2=$(grep -E 'TM-score=' "$f" | grep -m1 'Structure_2' \
    | sed -E 's/.*TM-score[= ]+([0-9.]+).*/\1/' | tr -d '\r') || true
  tmsavg=$(grep -E 'TM-score=' "$f" | grep -m1 'average length' \
    | sed -E 's/.*TM-score[= ]+([0-9.]+).*/\1/' | tr -d '\r') || true

  # TMscore_refL = TM-score normalised to reference length (Structure_2)
  tms_ref="$tms2"

  # Skip rows with no usable data
  if [[ -z "$aligned" && -z "$tms1" && -z "$tmsavg" ]]; then
    return 0
  fi

  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$raw_basename" \
    "$accession" \
    "$region" \
    "$oligomer_size" \
    "$sample" \
    "${aligned:-}" \
    "${rmsd:-}" \
    "${seqid:-}" \
    "${tms_ref:-}" \
    "${tms1:-}" \
    "${tms2:-}" \
    "${tmsavg:-}"
}

# ---------------------------------------------------------------------------
# Process all reference sets into temp files, sort, append to CSVs
# ---------------------------------------------------------------------------
process_ref() {
  local -n _cands="$1"
  local out_csv="$2"

  local tmp
  tmp="$(mktemp)"
  for f in "${_cands[@]}"; do
    extract_one "$f" >> "$tmp"
  done
  sort -t',' -k12 -rn "$tmp" >> "$out_csv"
  rm -f "$tmp"
}

process_ref CANDS_7SC0   "$out_7SC0"
process_ref CANDS_9DN1   "$out_9DN1"
process_ref CANDS_8BJM   "$out_8BJM"
process_ref CANDS_7D60   "$out_7D60"
process_ref CANDS_7RSL   "$out_7RSL"
process_ref CANDS_AFQ    "$out_AFQ"
process_ref CANDS_AFComp "$out_AFComp"

echo "Wrote: $out_7SC0  $out_9DN1  $out_8BJM  $out_7D60  $out_7RSL  $out_AFQ  $out_AFComp"

# ---------------------------------------------------------------------------
# Build SummaryHits tables (filtered and unfiltered)
# ---------------------------------------------------------------------------
HITS_HEADER="accession,region,oligomer_size,sample,reference,TMscore_refL,TMscore_avgL"
echo "$HITS_HEADER" > "$out_hits"
echo "$HITS_HEADER" > "$out_hits_all"

tmp_hits="$(mktemp)"
tmp_hits_all="$(mktemp)"

build_hits() {
  local csv_file="$1"
  local ref_name="$2"

  # Input columns (1-indexed):
  #   1=basename 2=accession 3=region 4=oligomer_size 5=sample
  #   6=aligned_len 7=RMSD 8=Seq_ID
  #   9=TMscore_refL 10=TMscore_L1 11=TMscore_L2 12=TMscore_avgL
  awk -F',' -v ref="$ref_name" \
      -v filt="$tmp_hits" -v all="$tmp_hits_all" '
  NR == 1 { next }   # skip header
  {
    tmsavg = $12
    row = $2 "," $3 "," $4 "," $5 "," ref "," $9 "," tmsavg
    print row >> all
    if (tmsavg + 0 > 0.2) print row >> filt
  }' "$csv_file"
}

build_hits "$out_7SC0"   "7SC0"
build_hits "$out_9DN1"   "9DN1"
build_hits "$out_8BJM"   "8BJM"
build_hits "$out_7D60"   "7D60"
build_hits "$out_7RSL"   "7RSL"
build_hits "$out_AFQ"    "AF-Q03135"
build_hits "$out_AFComp" "AFComp"

sort -t',' -k7 -rn "$tmp_hits"     >> "$out_hits"
sort -t',' -k7 -rn "$tmp_hits_all" >> "$out_hits_all"
rm -f "$tmp_hits" "$tmp_hits_all"

echo "Wrote: $out_hits"
n_hits=$(( $(wc -l < "$out_hits") - 1 ))
echo "  Total hits (TMscore_avgL > 0.2): ${n_hits}"

echo "Wrote: $out_hits_all"
n_all=$(( $(wc -l < "$out_hits_all") - 1 ))
echo "  Total entries (all scores):      ${n_all}"

# ---------------------------------------------------------------------------
# Report any CIFs that were skipped due to missing inference output
# ---------------------------------------------------------------------------
SKIP_LOG="${ALIGN_DIR}/skipped.log"
if [[ -f "$SKIP_LOG" ]]; then
  n_skipped=$(wc -l < "$SKIP_LOG" | tr -d '[:space:]')
  echo ""
  echo "WARNING: ${n_skipped} CIF(s) were skipped due to missing inference output."
  echo "  See: ${SKIP_LOG}"
  echo "  Skipped entries:"
  awk -F'\t' '{printf "    task=%-6s  %s\n", $3, $4}' "$SKIP_LOG"
  echo ""
  echo "  To re-run failed inference tasks and fill in missing alignments:"
  echo "    1. Check logs/af3infer_*.err for the failed task IDs"
  echo "    2. Re-submit inference for those tasks"
  echo "    3. Re-run US_align for those CIFs (idempotent — completed alignments are skipped)"
else
  echo "  No skipped CIFs (Align_Out/skipped.log not present — all inference tasks succeeded)"
fi
