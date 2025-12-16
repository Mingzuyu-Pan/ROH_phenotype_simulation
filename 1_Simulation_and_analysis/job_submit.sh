#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=2GB
#SBATCH --time=2-10:00:00
#SBATCH --account=zps5164_cr_default
#SBATCH --partition=basic
#SBATCH --mail-user=mzp5919@psu.edu
#SBATCH --mail-type=ALL

# Using squeue to check the status of a job
wait_for_job() {
    local job_id=$1
    echo "Waiting for $job_id to be done..."
    while [ -n "$(squeue -h -j "$job_id")" ]; do
        sleep 30
    done
    echo "Job $job_id is done."
}

# let multiple jobs wait
wait_for_jobs() {
    local dep_list=$1
    IFS=: read -r -a jobs <<< "$dep_list"
    for job in "${jobs[@]}"; do
        wait_for_job "$job"
    done
}

echo "Submit _1.sh..."
jid1_h00=$(sbatch --parsable command_file_general_1.sh h00)
jid1_h05=$(sbatch --parsable command_file_general_1.sh h05)
jid1_h10=$(sbatch --parsable command_file_general_1.sh h10)
jid1_h15=$(sbatch --parsable command_file_general_1.sh h15)

if [ -z "$jid1_h00" ] || [ -z "$jid1_h05" ] || [ -z "$jid1_h10" ] || [ -z "$jid1_h15" ]; then
    echo "Can't submit job _1.sh. Exiting."
    exit 1
fi

dep1="${jid1_h00}:${jid1_h05}:${jid1_h10}:${jid1_h15}"
echo "Job _1.sh is submitted: $dep1"
wait_for_jobs "$dep1"
echo "_1.sh stage completed."
）

echo "Submit _2.sh..."
jid2_h00=$(sbatch --parsable command_file_general_2.sh h00)
jid2_h05=$(sbatch --parsable command_file_general_2.sh h05)
jid2_h10=$(sbatch --parsable command_file_general_2.sh h10)
jid2_h15=$(sbatch --parsable command_file_general_2.sh h15)

if [ -z "$jid2_h00" ] || [ -z "$jid2_h05" ] || [ -z "$jid2_h10" ] || [ -z "$jid2_h15" ]; then
    echo "Can't submit job _2.sh. Exiting."
    exit 1
fi

dep2="${jid2_h00}:${jid2_h05}:${jid2_h10}:${jid2_h15}"
echo "Job _2.sh is submitted: $dep2"
wait_for_jobs "$dep2"
echo "_2.sh stage completed."


echo "Submit job _3.sh..."
jid3_h00=$(sbatch --parsable command_file_general_3.sh h00)
jid3_h05=$(sbatch --parsable command_file_general_3.sh h05)
jid3_h10=$(sbatch --parsable command_file_general_3.sh h10)
jid3_h15=$(sbatch --parsable command_file_general_3.sh h15)

if [ -z "$jid3_h00" ] || [ -z "$jid3_h05" ] || [ -z "$jid3_h10" ] || [ -z "$jid3_h15" ]; then
    echo "Can't submit job _3.sh. Exiting."
    exit 1
fi

dep3="${jid3_h00}:${jid3_h05}:${jid3_h10}:${jid3_h15}"
echo "Job _3.sh is submitted: $dep3"
wait_for_jobs "$dep3"
echo "_3.sh stage completed."

echo "Submit job _4.sh..."
jid4_pre=$(sbatch --parsable command_file_general_4.sh)
if [ -z "$jid4_pre" ]; then
    echo "Can't submit job _4.sh. Exiting."
    exit 1
fi
wait_for_job "$jid4_pre"
echo "_4.sh stage completed."

echo "Submit job pheno..."
jid_pheno_h00=$(sbatch --parsable run_pheno_neutral_all_commands_test.sh h00)
jid_pheno_h05=$(sbatch --parsable run_pheno_neutral_all_commands_test.sh h05)
jid_pheno_h10=$(sbatch --parsable run_pheno_neutral_all_commands_test.sh h10)
jid_pheno_h15=$(sbatch --parsable run_pheno_neutral_all_commands_test.sh h15)

if [ -z "$jid_pheno_h00" ] || [ -z "$jid_pheno_h05" ] || [ -z "$jid_pheno_h10" ] || [ -z "$jid_pheno_h15" ]; then
    echo "Can't submit job pheno. Exiting."
    exit 1
fi

dep_pheno="${jid_pheno_h00}:${jid_pheno_h05}:${jid_pheno_h10}:${jid_pheno_h15}"
echo "Job pheno is submitted: $dep_pheno"
wait_for_jobs "$dep_pheno"
echo "Pheno stage completed."

echo "Submit job _5.sh..."
jid5=$(sbatch --parsable command_file_general_5.sh)
if [ -z "$jid5" ]; then
    echo "Can't submit job _5.sh. Exiting."
    exit 1
fi
wait_for_job "$jid5"
echo "_5.sh stage completed."

echo "All jobs submitted successfully."
echo "All stages completed successfully."