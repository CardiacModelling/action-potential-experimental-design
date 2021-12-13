#!/bin/bash
#SBATCH --partition             fhs-fast
#SBATCH --ntasks                1
#SBATCH --nodes                 1
#SBATCH --tasks-per-node        1
#SBATCH --cpus-per-task         48
#SBATCH --mem-per-cpu           4G
#SBATCH --time                  120:00:00
#SBATCH --job-name              x-cri-eva
#SBATCH --output                log/x-cri-eva.%j.out
#SBATCH --error                 log/x-cri-eva.%j.err
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

# We are using multiprocessing, so switch multi-threading off
# https://stackoverflow.com/a/43897781
# export OMP_NUM_THREADS=1

# Run

python -u cross-criteria-evaluate.py

echo "Done."
