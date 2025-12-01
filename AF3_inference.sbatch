#!/bin/bash
#SBATCH -A cavorigins
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH -C "a6000|a100"  # compatible with A6000, A100
#SBATCH -t 3-00:00:00
#SBATCH --job-name=AF3_batch
#SBATCH --array=3000-3499%4  # Adjust based on number of input files
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

module purge
module load gcc/11.4.0 alphafold/3.0.0
INPUT_DIR=$(pwd -P)/AF3_Data
OUTPUT_DIR=$(pwd -P)/OUTPUT_MODELS
mkdir -p "$OUTPUT_DIR" logs

# Get the list of input files
#INPUT_FILES=($INPUT_DIR/*.json)  # or *.fasta if using FASTA
#INPUT_FILE=${INPUT_FILES[$SLURM_ARRAY_TASK_ID]}
INPUT_FILES=($(find "$INPUT_DIR" -type f -name "*_data.json" | sort))
INPUT_FILE=${INPUT_FILES[$SLURM_ARRAY_TASK_ID]}

python $EBROOTALPHAFOLD/app/alphafold/run_alphafold.py \
    --run_data_pipeline=false \
    --db_dir=$ALPHAFOLD_DATA_PATH \
    --jax_compilation_cache_dir /scratch/rts3dd/jaxcache \
    --model_dir=/scratch/rts3dd/AF3 \
    --json_path="$INPUT_FILE" \
    --output_dir="$OUTPUT_DIR"

