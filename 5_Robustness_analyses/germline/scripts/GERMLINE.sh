#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=15GB
#SBATCH --array=1-500
#SBATCH --time=10:00:00
#SBATCH --account=zps5164_cr_default
#SBATCH --partition=himem
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

BCFTOOLS=/storage/group/zps5164/default/bin/bcftools
GERMLINE_DIR=/storage/group/zps5164/default/mzp5919/complex_traits_comparison/my_data_collection/out_of_africa/germline-1-5-3

for POP in p1 p2 p3; do
    PREFIX="./${hXX}/${seed}/gravel.del.only.${hXX}.${seed}.trees.${POP}"
    ${BCFTOOLS} convert --hapsample ${PREFIX} ${PREFIX}.vcf.gz
    gunzip -f ${PREFIX}.hap.gz
    awk 'BEGIN{srand(); split("A,G,C,T",a,",")} {r=int(rand()*2); if(r==0){$4="A";$5="T"}else{$4="G";$5="C"}; print}' ${PREFIX}.hap > ${PREFIX}.fixed.hap
    ${GERMLINE_DIR}/bin/impute_to_ped ${PREFIX}.fixed.hap ${PREFIX}.samples ${PREFIX}.fixed
    ${GERMLINE_DIR}/germline \
    -input ${PREFIX}.fixed.ped \
            ${PREFIX}.fixed.map \
    -homoz-only \
    -min_m 0.0000001 \
    -bits 128 \
    -output ${PREFIX}.fixed.germline.roh

    gzip -f ${PREFIX}.hap ${PREFIX}.fixed.hap ${PREFIX}.fixed.ped ${PREFIX}.fixed.map
done