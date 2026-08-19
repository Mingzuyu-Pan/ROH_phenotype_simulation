#!/bin/bash
#SBATCH --job-name=filter_chrs2-22
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4GB
#SBATCH --array=2-22
#SBATCH --time=2:00:00
#SBATCH --account=zps5164_cr_default
#SBATCH --partition=basic
#SBATCH --mail-user=mzp5919@psu.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=filter_%A_chr%a.log 
#SBATCH --error=filter_%A_chr%a.err

CHR=$SLURM_ARRAY_TASK_ID
echo "Now processing Chromosome: $CHR"

INPUT="ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
OUTPUT="ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.biallelic_snps.vcf.gz"

/storage/group/zps5164/default/bin/bcftools view --threads 4 -m2 -M2 -v snps $INPUT -O z -o $OUTPUT