#!/bin/bash
#SBATCH -A cavorigins           # your allocation account
#SBATCH -p standard               # partition
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH -t 12:00:00          # time
#SBATCH -o logs/alphafold_%A_%a.out
#SBATCH -e logs/alphafold_%A_%a.err
#SBATCH --array=0-4299%8  # Adjust based on number of input files

module purge
module load gcc/11.4.0 alphafold

# Optional: create output directory
mkdir -p AF3_Data
mkdir -p logs



# Get the list of input files
INPUT_FILES=(AF3_Seqs2/*.json)  # or *.fasta if using FASTA
INPUT_FILE=${INPUT_FILES[$SLURM_ARRAY_TASK_ID]}

python $EBROOTALPHAFOLD/app/alphafold/run_alphafold.py \
    --run_inference=false \
    --db_dir=$ALPHAFOLD_DATA_PATH \
    --jackhmmer_n_cpu=32 \
    --nhmmer_n_cpu=32 \
    --jax_compilation_cache_dir /scratch/rts3dd/jaxcache \
    --model_dir=/scratch/rts3dd/AF3 \
    --json_path="$INPUT_FILE" \
    --output_dir=$(pwd -P)/AF3_Data
