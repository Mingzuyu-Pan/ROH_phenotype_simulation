#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --array=1
#SBATCH --time=1:00:00
#SBATCH --account=zps5164_cr_default
#SBATCH --partition=basic
#SBATCH --mail-user=mzp5919@psu.edu
#SBATCH --mail-type=ALL


python del_mutation_calculate_per_replicate.py h00 
python del_mutation_calculate_per_replicate.py h05
python del_mutation_calculate_per_replicate.py h10
python del_mutation_calculate_per_replicate.py h15

python to_get_distribution.py h00 h05 h10 h15
