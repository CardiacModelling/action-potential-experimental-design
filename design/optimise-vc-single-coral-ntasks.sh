#!/bin/bash
#SBATCH --partition             fhs-fast
#SBATCH --ntasks                2
#SBATCH --nodes                 1
#SBATCH --tasks-per-node        2
#SBATCH --cpus-per-task         24
#SBATCH --time                  100:00:00
#SBATCH --mem                   60G
#SBATCH --job-name              OED-GSA-A
#SBATCH --output                log/OED-GSA-A-.%j.out
#SBATCH --error                 log/OED-GSA-A-.%j.err
#SBATCH --mail-type		ALL
##SBATCH --mail-user		chonloklei@um.edu.mo

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
python -u optimise-vc-single.py -d LSA-E -l ../mmt/ohara-2011.mmt -n 10
DESIGN="GSA-A"
MODEL="ohara-2011"
PATH2MODEL="../mmt/${MODEL}.mmt"
N=5
LOG="./log/optimise-vc-single"
mkdir -p $LOG

# We are using multiprocessing, so switch multi-threading off
# https://stackoverflow.com/a/43897781
export OMP_NUM_THREADS=1

# Run
for ((x=0; x<2; x++));
do
	echo "${CELL_LIST[x]}"
	# NOTE: We are not doing MPI here, so it takes only 1 task, but needs more than 1 CPU to do multiprocessing.
	#       https://login.scg.stanford.edu/faqs/cores/#nodes-vs-tasks-vs-cpus-vs-cores
	srun --exclusive --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --mem=30G python -u optimise-vc-single.py -d ${DESIGN} -l ${PATH2MODEL} -n ${N} -r ${x} -p ${SLURM_CPUS_PER_TASK} --tmp > ${LOG}/${DESIGN}-${MODEL}.log &
done

wait
