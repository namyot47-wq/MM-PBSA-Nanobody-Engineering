#!/bin/bash
#SBATCH -p public
#SBATCH -q public
#SBATCH -G a100:1
#SBATCH -c 4
#SBATCH -t 0-4
#SBATCH -J nb_equil
#SBATCH -o slurm_logs/equil_%j.out
#SBATCH -e slurm_logs/equil_%j.err

set -eo pipefail

module load mamba/latest
source activate amber_pipeline
module load amber/24        
unset PYTHONPATH

echo "CONDA_DEFAULT_ENV=$CONDA_DEFAULT_ENV"
which python3 cpptraj sander pmemd.cuda
python3 -c "import yaml; print('yaml OK')"

cd /scratch/namyot/BartelleLab/MM-PBSA-Nanobody-Engineering/amber_pipeline
python -m pipeline.run_pipeline config.yaml --step equil
