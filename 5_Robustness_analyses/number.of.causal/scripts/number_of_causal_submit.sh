#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --array=1-500
#SBATCH --time=10:00:00
#SBATCH --account=zps5164_cr_default
#SBATCH --partition=basic
#SBATCH --mail-user=mzp5919@psu.edu
#SBATCH --mail-type=ALL

# read hXX from command line
hXX=$1

# to check if the script reads it successfully
if [ -z "$hXX" ]; then
    echo "Error: hXX value not provided."
    echo "Usage: sbatch run.sh <hXX>"
    exit 1
fi

# to read the seed
seed=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${hXX}_500.list")

echo "Task ID: $SLURM_ARRAY_TASK_ID, Seed: $seed"

module load python/3.11.2
python number_of_causal.py "$hXX" "$seed" 1.0
python number_of_causal.py "$hXX" "$seed" 0.3