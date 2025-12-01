
#!/bin/bash
#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH -J "Make_CD_input"
#SBATCH --cpus-per-task=16
#SBATCH --partition=standard
#SBATCH --mail-user=rts3dd@virginia.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --account=cavorigins

###############################################################################
# Script Name: Make_CD_input
#
# Purpose:
#   1. Reads an input file containing accession numbers and start/end positions.
#   2. Retrieves sequences from BLAST databases (nr, env_nr, tsa_nr) using
#      blastdbcmd, expanding each range by ±50 amino acids.
#   3. Writes all sequences to a FASTA file.
#   4. Clusters sequences using CD-HIT at 95% identity.
#   5. Extracts representative accession IDs from the CD-HIT cluster file.
#
# Input:
#   - A tab-delimited file with columns: ACCESSION START END
#
# Output:
#   - FASTA file of retrieved sequences: sequences_##.fasta
#   - CD-HIT clustered FASTA: clustered_sequences_##.fasta
#   - Representative accession list: clstr_accessions_##.txt
#   - Log file: blastdbcmd_##.log
#
# Parameters:
#   - Expansion: ±50 amino acids around the specified range
#   - CD-HIT clustering: 95% identity (-c 0.95), word size (-n 5)
#   - Threads: 16
#
# Notes:
#   - If expanded start < 1, it is clamped to 1.
#   - If expanded end exceeds sequence length, BLAST truncates automatically.
###############################################################################

# Load BLAST module
ml load blast

# Check for input file argument
if [ $# -ne 1 ]; then
    echo "Usage: $0 <accession_file>"
    exit 1
fi

# Determine output file names
counter=1
OUTPUT_FASTA="sequences_$(printf %02d $counter).fasta"

while [ -f "$OUTPUT_FASTA" ]; do
    counter=$((counter + 1))
    OUTPUT_FASTA="sequences_$(printf %02d $counter).fasta"
done

ACCESSION_FILE="$1"
LOG_FILE="blastdbcmd_$(printf %02d $counter).log"

# Loop through each line of the input file
while IFS=$'\t' read -r ACCESSION START END; do
    echo "Processing $ACCESSION ($START-$END)" >> "$LOG_FILE"

    # Expand range by 50 amino acids on each side
    EXPANDED_START=$((START - 50))
    EXPANDED_END=$((END + 50))

    # Ensure start is at least 1
    if [ "$EXPANDED_START" -lt 1 ]; then
        EXPANDED_START=1
    fi

    # Retrieve sequence using blastdbcmd
    blastdbcmd -db "nr env_nr tsa_nr" \
               -entry "$ACCESSION" \
               -range "$EXPANDED_START"-"$EXPANDED_END" \
               -outfmt %f >> "$OUTPUT_FASTA" 2>> "$LOG_FILE"

done < "$ACCESSION_FILE"

echo "BLASTDBCMD completed. Output written to $OUTPUT_FASTA. Log written to $LOG_FILE."

# Load CD-HIT module
module purge
module load cd-hit

# Run CD-HIT clustering
cd-hit -i "$OUTPUT_FASTA" \
       -o "clustered_sequences_$(printf %02d $counter).fasta" \
       -c 0.95 -n 5 -T 16 -d 0

# Extract representative accessions from .clstr file
grep '\*' "clustered_sequences_$(printf %02d $counter).fasta.clstr" \
    | sed -E 's/.*>([^ ]+)\.\.\..*/\1/' \
    > "clstr_accessions_$(printf %02d $counter).txt"

echo "Representative accessions written to clstr_accessions_$(printf %02d $counter).txt"
