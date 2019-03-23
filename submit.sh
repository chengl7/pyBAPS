#!/bin/bash
##SBATCH --time=0-03:25:00    # 5 mins
#SBATCH --time=0-00:05:00    # 5 mins
#SBATCH -p short
#SBATCH -p batch
#SBATCH -p debug
#SBATCH -o a.out
#SBATCH -e a.err
#SBATCH -N 8
##SBATCH -n 4
#SBATCH --exclusive
##SBATCH --mem=40000    # 500MB of memory

# output goes into hello.out
# If you use srun for each command, the mem/cpu usage of each step
# can be seen individually with "slurm history"
#echo $SLURM_JOB_NODELIST
#srun hostname

module load anaconda3

if [ $# -eq 0 ]; then
    echo "Usage: sbatch submit.sh X1K.npy test1K X10K.npy test10K"
    exit 1
fi
exit 1

server=$(echo "$SLURM_JOB_NODELIST" | cut -d ',' -f1)
server=$(echo "$server" | cut -d '-' -f1)
server=$(echo "$server" | tr -d '[')
echo "list is ${SLURM_JOB_NODELIST} and server is ${server}"
echo "n jobs is ${SLURM_JOB_NUM_NODES}"
echo " " 

start=`date +%s`
echo "input parameters: $@"
srun script.sh ${SLURM_JOB_NUM_NODES} "$server" "$@" 
end=`date +%s`

runtime=$((end-start))
echo "runtime: $runtime seconds"
