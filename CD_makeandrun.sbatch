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

ml load blast

# Check for input file argument
if [ $# -ne 1 ]; then
    echo "Usage: $0 <accession_file>"
    exit 1
fi



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
    
    blastdbcmd -db "nr env_nr tsa_nr" \
               -entry "$ACCESSION" \
               -range "$START"-"$END" \
               -outfmt %f >> "$OUTPUT_FASTA" 2>> "$LOG_FILE"
done < "$ACCESSION_FILE"

echo "BLASTDBCMD completed. Output written to $OUTPUT_FASTA. Log written to $LOG_FILE."


# Load CD-HIT module if required (uncomment if using module system)
module purge
module load cd-hit

# Run CD-HIT clustering
cd-hit -i "$OUTPUT_FASTA" -o "clustered_sequences_$(printf %02d $counter).fasta" -c 0.95 -n 5 -T 16 -d 0

# Extract representative accessions from .clstr file
grep '\*' "clustered_sequences_$(printf %02d $counter).fasta.clstr" | sed -E 's/.*>([^ ]+)\.\.\..*/\1/' > "clstr_accessions_$(printf %02d $counter).txt"
echo "Representative accessions written to clstr_accessions_$(printf %02d $counter).txt"
