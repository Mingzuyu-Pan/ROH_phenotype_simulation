#!/bin/bash
#SBATCH --job-name=concat_naive
#SBATCH --array=1-17
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1           
#SBATCH --mem=4GB                    
#SBATCH --time=1:00:00              
#SBATCH --account=zps5164_cr_default
#SBATCH --partition=basic
#SBATCH --mail-user=mzp5919@psu.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=concat_%A_pop%a.log
#SBATCH --error=concat_%A_pop%a.err

POP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" pop_list.txt)
echo "Working on pop: $POP"

PREFIX="ALL.chr"
SUFFIX=".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.biallelic_snps_${POP}.vcf.gz"

FILE_LIST=""
for i in {1..22}; do
    FILE_LIST="$FILE_LIST ${PREFIX}${i}${SUFFIX}"
done

OUT_FILE="WholeGenome_chr1-22_${POP}.vcf.gz"

/storage/group/zps5164/default/bin/bcftools concat  --naive -o $OUT_FILE $FILE_LIST
/storage/group/zps5164/default/bin/tabix -p vcf $OUT_FILE

echo "Finish working on pop: $POP"