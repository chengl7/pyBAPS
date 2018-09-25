#!/bin/bash
#SBATCH --time=0-00:05:00    # 5 mins
#SBATCH -p short
##SBATCH -p debug
#SBATCH -o a.out
#SBATCH -N 4
##SBATCH -n 4
#SBATCH --mem=500    # 500MB of memory

# output goes into hello.out
# If you use srun for each command, the mem/cpu usage of each step
# can be seen individually with "slurm history"
#echo $SLURM_JOB_NODELIST
#srun hostname

server=$(echo "$SLURM_JOB_NODELIST" | cut -d ',' -f1)
server=$(echo "$server" | tr -d '[')
echo "list is ${SLURM_JOB_NODELIST} and server is ${server}"

srun script.sh "$server"
