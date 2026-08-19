#!/bin/bash
#SBATCH --job-name=split
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --array=2-22
#SBATCH --time=6:00:00
#SBATCH --account=zps5164_cr_default
#SBATCH --partition=basic
#SBATCH --mail-user=mzp5919@psu.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=filter_%A_chr%a.log 
#SBATCH --error=filter_%A_chr%a.err

chr=$SLURM_ARRAY_TASK_ID
echo "Now processing Chromosome: $chr"

pop_map="pop_map_bcftools.txt"
prefix="ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.biallelic_snps"
input="${prefix}.vcf.gz"

temp_dir="temp_split_chr${chr}"
mkdir -p $temp_dir

/storage/group/zps5164/default/bin/bcftools +split $input -G $pop_map -O z -o $temp_dir/ 

for pop_file in $temp_dir/*.vcf.gz; do
    pop_name=$(basename "$pop_file" .vcf.gz)
    new_name="${prefix}_${pop_name}.vcf.gz"
    mv "$pop_file" "./$new_name"
done

rm -r $temp_dir