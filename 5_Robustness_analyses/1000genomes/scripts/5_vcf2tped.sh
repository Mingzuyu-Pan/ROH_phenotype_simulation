#!/bin/bash
#SBATCH --job-name=vcf2tped
#SBATCH --array=1-17
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4 
#SBATCH --mem=8GB 
#SBATCH --time=1:00:00               
#SBATCH --account=zps5164_cr_default
#SBATCH --partition=basic
#SBATCH --mail-user=mzp5919@psu.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=convert_%A_pop%a.log
#SBATCH --error=convert_%A_pop%a.err


POP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" pop_list.txt)
echo "Working on pop: $POP"

INPUT_VCF="WholeGenome_chr1-22_${POP}.vcf.gz"
OUT_PREFIX="WholeGenome_chr1-22_${POP}"

/storage/group/zps5164/default/bin/plink  --vcf $INPUT_VCF --recode transpose --const-fid --maf 0.05 --out $OUT_PREFIX

echo "Compress tped now..."

gzip -f ${OUT_PREFIX}.tped

echo "Finish working on $POP"
