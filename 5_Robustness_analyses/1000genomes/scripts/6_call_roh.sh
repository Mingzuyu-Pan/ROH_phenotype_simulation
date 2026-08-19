#!/bin/bash
#SBATCH --job-name=call_roh
#SBATCH --array=1-17
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB 
#SBATCH --time=10:00:00               
#SBATCH --account=zps5164_cr_default
#SBATCH --partition=himem
#SBATCH --mail-user=mzp5919@psu.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=call_roh_%A_pop%a.log
#SBATCH --error=call_roh_%A_pop%a.err


POP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" pop_list.txt)
echo "Working on pop: $POP"
PREFIX="WholeGenome_chr1-22_${POP}"
/storage/group/zps5164/default/bin/garlic --tped ${PREFIX}.tped.gz --tfam ${PREFIX}.tfam --cm --map chraut.map --build hg19 --size-bounds 0.25 1 --winsize 10 --auto-winsize --auto-overlap-frac --out ${PREFIX}_ROH --error 0.0001 --lod-cutoff 0 --tped-missing 9