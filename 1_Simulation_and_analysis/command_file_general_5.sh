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

module load python/3.11.2
module use /storage/icds/RISE/sw8/modules/r
module load r/4.2.1-gcc-8.5.0

Rscript to_normalize_phenotype_score.R h00
Rscript to_normalize_phenotype_score.R h05
Rscript to_normalize_phenotype_score.R h10
Rscript to_normalize_phenotype_score.R h15

python process_data_frac_deleterious.py
python process_data_frac_neutral.py
python process_data_proportion_deleterious.py
python process_data_proportion_neutral.py
