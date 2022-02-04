#!/bin/bash
#SBATCH --job-name              oed-mcmc
#SBATCH --partition             serial-normal
#SBATCH --nodes                 1
#SBATCH --tasks-per-node        4
#SBATCH --array                 0-5
#SBATCH --time                  72:00:00
#SBATCH --mem                   30G
#SBATCH --output                log/oed-mcmc-ohara-cc.%A_%a.out
#SBATCH --error                 log/oed-mcmc-ohara-cc.%A_%a.err
#SBATCH --mail-type		        ALL
##SBATCH --mail-user		    chonloklei@um.edu.mo

source /etc/profile
source /etc/profile.d/modules.sh
source /home/chonloklei/m  # Load miniconda

ulimit -s unlimited

# Path and Python version checks
pwd
python --version
conda activate /home/chonloklei/action-potential-experimental-design/env  # Load miniconda venv
python --version
which python

# Set up
PATH2FILES="./practicality-input/"
FILE="${PATH2FILES}true_ohara-fit_ohara-row_models-col_measures-cc.txt"
declare -A m

read_matrix() {
    local i=0
    local line
    local j
    local J=0
    while read -r line; do
        j=0
        # split on spaces
        for v in `echo $line`; do
            # Not doing GSA for now!
            if [ $j -lt 3 ]
            then
                V=$(printf '%03d' "$v")
                m[$((i*J+j))]="$V"
                j=$((j+1))
            fi
        done
        # Get the size of column
        if [ $i -lt 1 ]
        then
            J=$j
        fi
        i=$((i+1))
    done
}

read_matrix < $FILE

ID=${m[${SLURM_ARRAY_TASK_ID}]}

# Run
F=${PATH2FILES}id_${ID}.py
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "ID: $ID"
echo "Running MCMC with $F"
python -u practicality-mcmc-cc.py $F
