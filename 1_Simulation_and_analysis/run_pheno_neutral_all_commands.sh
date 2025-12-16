#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=3GB
#SBATCH --time=1-10:00:00
#SBATCH --array=1-500
#SBATCH --account=zps5164_cr_default
#SBATCH --partition=basic
#SBATCH --mail-user=mzp5919@psu.edu
#SBATCH --mail-type=ALL

hXX=$1

if [ -z "$hXX" ]; then
    echo "Error: hXX value not provided."
    echo "Usage: sbatch run.sh <hXX>"
    exit 1
fi

seed=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${hXX}_500.list")

module load python/3.11.2

# to get the file path
files_p1=("./${hXX}/${seed}"/*.p1.vcf.gz)
files_p2=("./${hXX}/${seed}"/*.p2.vcf.gz)
files_p3=("./${hXX}/${seed}"/*.p3.vcf.gz)


# to check if we find the file we want for this step
if [ "${#files_p1[@]}" -eq 0 ]; then
    echo "No .p1.vcf.gz files found in ./${hXX}/${seed}/"
    exit 1
fi

# to check if we find the file we want for this step
if [ "${#files_p2[@]}" -eq 0 ]; then
    echo "No .p2.vcf.gz files found in ./${hXX}/${seed}/"
    exit 1
fi

# to check if we find the file we want for this step
if [ "${#files_p3[@]}" -eq 0 ]; then
    echo "No .p3.vcf.gz files found in ./${hXX}/${seed}/"
    exit 1
fi

gunzip  "${files_p1[@]}"
gunzip  "${files_p2[@]}"
gunzip  "${files_p3[@]}"


python to_get_neutral_mutations_in_exon.py "$hXX" "$seed"
python sample_mutations_from_all_in_exon.py "$hXX" "$seed"

sampled_files_p1=("./${hXX}/${seed}"/*.p1.exon.in.smapled.vcf)
sampled_files_p2=("./${hXX}/${seed}"/*.p2.exon.in.smapled.vcf)
sampled_files_p3=("./${hXX}/${seed}"/*.p3.exon.in.smapled.vcf)

# to check if we find the file we want for this step
if [ "${#sampled_files_p1[@]}" -eq 0 ]; then
    echo "No .p1.exon.in.smapled.vcf found in ./${hXX}/${seed}/"
    exit 1
fi

if [ "${#sampled_files_p2[@]}" -eq 0 ]; then
    echo "No .p2.exon.in.smapled.vcf found in ./${hXX}/${seed}/"
    exit 1
fi

if [ "${#sampled_files_p3[@]}" -eq 0 ]; then
    echo "No .p3.exon.in.smapled.vcf found in ./${hXX}/${seed}/"
    exit 1
fi


# to check if there are more than one ALT values
if grep -v '^#' "${sampled_files_p1[@]}" | awk -F'\t' '$5 ~ /,/' | grep -q .; then
    echo "$hXX Seed $seed contains multiple ALT values in p1. Job ID: $SLURM_ARRAY_JOB_ID, Task ID: $SLURM_ARRAY_TASK_ID" | mail -s "VCF ALT Check Alert (p1)" mzp5919@psu.edu
    echo "$hXX Seed $seed contains multiple ALT values in p1. Job ID: $SLURM_ARRAY_JOB_ID, Task ID: $SLURM_ARRAY_TASK_ID"
else
    echo "$hXX Seed $seed does not contain multiple ALT values in p1."
fi

if grep -v '^#' "${sampled_files_p2[@]}" | awk -F'\t' '$5 ~ /,/' | grep -q .; then
    echo "$hXX Seed $seed contains multiple ALT values in p2. Job ID: $SLURM_ARRAY_JOB_ID, Task ID: $SLURM_ARRAY_TASK_ID" | mail -s "VCF ALT Check Alert (p2)" mzp5919@psu.edu
    echo "$hXX Seed $seed contains multiple ALT values in p2. Job ID: $SLURM_ARRAY_JOB_ID, Task ID: $SLURM_ARRAY_TASK_ID"
else
    echo "$hXX Seed $seed does not contain multiple ALT values in p2."
fi

if grep -v '^#' "${sampled_files_p3[@]}" | awk -F'\t' '$5 ~ /,/' | grep -q .; then
    echo "$hXX Seed $seed contains multiple ALT values in p3. Job ID: $SLURM_ARRAY_JOB_ID, Task ID: $SLURM_ARRAY_TASK_ID" | mail -s "VCF ALT Check Alert (p3)" mzp5919@psu.edu
    echo "$hXX Seed $seed contains multiple ALT values in p3. Job ID: $SLURM_ARRAY_JOB_ID, Task ID: $SLURM_ARRAY_TASK_ID"
else
    echo "$hXX Seed $seed does not contain multiple ALT values in p3."
fi

for tau in 70 80 90 100 110 120 130; do
    python get_feature_file_for_neutral.py "$hXX" "$seed" 0 "$tau"
    sleep 10
done

for tau in 70 80 90 100 110 120 130; do
	python neutral_pheno_modified.py "$hXX" "$seed" 0 "$tau"
	sleep 10
done

sleep 30

find "./${hXX}/${seed}" -maxdepth 1 -name "*.vcf" -type f -exec gzip -f {} \;


