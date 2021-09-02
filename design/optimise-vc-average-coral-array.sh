#!/bin/bash
#SBATCH --partition             serial-normal
#SBATCH --ntasks                1
#SBATCH --nodes                 1
#SBATCH --tasks-per-node        1
#SBATCH --cpus-per-task         4
#SBATCH --array                 0-9
#SBATCH --time                  72:00:00
#SBATCH --mem                   62G
#SBATCH --job-name              a-GSA-D
#SBATCH --output                log/a-GSA-D.%A_%a.out
#SBATCH --error                 log/a-GSA-D.%A_%a.err
#SBATCH --mail-type             ALL
##SBATCH --mail-user            chonloklei@um.edu.mo

source /etc/profile
source /etc/profile.d/modules.sh
source /home/chonloklei/m  # Load miniconda

ulimit -s unlimited

# Load module
module purge

# Path and Python version checks
pwd
python --version
conda activate /home/chonloklei/action-potential-experimental-design/env  # Load miniconda venv
python --version
which python

# Set up
DESIGN="GSA-D"
MODEL="model-list"
PATH2MODEL="./${MODEL}.txt"
N=1
LOG="./log/optimise-vc-average"
mkdir -p $LOG

# We are using multiprocessing, so switch multi-threading off
# https://stackoverflow.com/a/43897781
# export OMP_NUM_THREADS=1

# Run
python -u optimise-vc-average.py -d ${DESIGN} -l ${PATH2MODEL} -n ${N} -r ${SLURM_ARRAY_TASK_ID} -p -1 --tmp --rand_x0 > ${LOG}/${DESIGN}_${MODEL}_n${N}_x${SLURM_ARRAY_TASK_ID}.log

echo "Done."
