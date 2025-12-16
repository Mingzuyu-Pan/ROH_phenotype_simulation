#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20GB
#SBATCH --array=1-500
#SBATCH --time=6:00:00
#SBATCH --account=zps5164_cr_default
#SBATCH --partition=himem
#SBATCH --mail-user=mzp5919@psu.edu
#SBATCH --mail-type=ALL

# read hXX from command line
hXX=$1

# to check if we read it successfully
if [ -z "$hXX" ]; then
    echo "Error: hXX value not provided."
    echo "Usage: sbatch run.sh <hXX>"
    exit 1
fi

module load python/3.11.2
python run_sim_modified.py gravel_chr1_100Mbp_treeseq_del_only_"$hXX".slim 
sleep 30


