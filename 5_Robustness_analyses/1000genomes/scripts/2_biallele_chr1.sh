#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4GB
#SBATCH --array=1
#SBATCH --time=2:00:00
#SBATCH --account=zps5164_cr_default
#SBATCH --partition=basic
#SBATCH --mail-user=mzp5919@psu.edu
#SBATCH --mail-type=ALL

/storage/group/zps5164/default/bin/bcftools view --threads 4 -m2 -M2 -v snps ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -O z -o ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.biallelic_snps.vcf.gz