#!/bin/bash
#SBATCH --partition             fhs-fast
#SBATCH --ntasks                1
#SBATCH --nodes                 1
#SBATCH --tasks-per-node        1
#SBATCH --cpus-per-task         48
#SBATCH --mem-per-cpu           4G
#SBATCH --time                  120:00:00
#SBATCH --job-name              OED-cc-GSA-E
#SBATCH --output                log/OED-cc-GSA-E.%j.out
#SBATCH --error                 log/OED-cc-GSA-E.%j.err
#SBATCH --mail-type             ALL
##SBATCH --mail-user            chonloklei@um.edu.mo

source /etc/profile
source /etc/profile.d/modules.sh
source /home/chonloklei/m  # Load miniconda

ulimit -s unlimited

# Load module
module purge
module load intel impi
# NOTE: Running multiple MPI jobs per node is currently only possible with IntelMPI and MVAPICH2 (i.e. OpenMPI does not work).
#       https://scitas-data.epfl.ch/confluence/display/DOC/Running+multiple+tasks+on+one+node

# Path and Python version checks
pwd
python --version
conda activate /home/chonloklei/action-potential-experimental-design/env  # Load miniconda venv
python --version
which python

# Set up
DESIGN="GSA-E"
MODEL="ohara-2011"
PATH2MODEL="../mmt/${MODEL}.mmt"
N=10
LOG="./log/optimise-cc-single"
mkdir -p $LOG

# We are using multiprocessing, so switch multi-threading off
# https://stackoverflow.com/a/43897781
# export OMP_NUM_THREADS=1

# Run

#for ((x=0; x<2; x++));
#do
#	# NOTE: We are not doing MPI here, so it takes only 1 task, but needs more than 1 CPU to do multiprocessing.
#	#       https://login.scg.stanford.edu/faqs/cores/#nodes-vs-tasks-vs-cpus-vs-cores
#    echo "Logging to ${LOG}/${DESIGN}_${MODEL}_n${N}_x${x}.log"
#	srun --exclusive --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --mem=45G python -u optimise-vc-single.py -d ${DESIGN} -l ${PATH2MODEL} -n ${N} -r ${x} -p ${SLURM_CPUS_PER_TASK} --tmp > ${LOG}/${DESIGN}_${MODEL}_n${N}_x${x}.log &
#done
#wait

python -u optimise-cc-single.py -d ${DESIGN} -l ${PATH2MODEL} -n ${N} > ${LOG}/${DESIGN}_${MODEL}_n${N}.log

echo "Done."
