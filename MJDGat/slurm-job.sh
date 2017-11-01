#!/bin/bash -l
#SBATCH -t 8:00:00  --ntasks=1
#SBATCH --mem 3400
#SBATCH --account=majorana
#SBATCH --workdir=
#SBATCH --output=./logs/slurm-%j.txt

echo "Job Start:"
date
echo "Node(s):  "$SLURM_JOB_NODELIST
echo "Job ID:  "$SLURM_JOB_ID

# Pass command here
echo "${@}"
${@}

echo "Job Complete:"
date