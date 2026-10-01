#!/bin/bash
#SBATCH -p public
#SBATCH -q public
#SBATCH -c 4
#SBATCH -t 0-8
#SBATCH -J nb_mmpbsa
#SBATCH -o slurm_logs/mmpbsa_%j.out
#SBATCH -e slurm_logs/mmpbsa_%j.err

set -eo pipefail

module load mamba/latest
source activate amber_pipeline
module load amber/24        
unset PYTHONPATH

echo "CONDA_DEFAULT_ENV=$CONDA_DEFAULT_ENV"
which python3 cpptraj sander pmemd.cuda
python3 -c "import yaml; print('yaml OK')"

cd /scratch/namyot/BartelleLab/MM-PBSA-Nanobody-Engineering/amber_pipeline
python -m pipeline.run_pipeline config.yaml --step mmpbsa