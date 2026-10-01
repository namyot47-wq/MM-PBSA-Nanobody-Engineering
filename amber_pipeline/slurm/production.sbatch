#!/bin/bash
#SBATCH -p public
#SBATCH -q public
#SBATCH -G a100:1
#SBATCH -c 4
#SBATCH -t 1-0
#SBATCH -J nb_prod
#SBATCH -o slurm_logs/prod_%j.out
#SBATCH -e slurm_logs/prod_%j.err

set -eo pipefail

module load mamba/latest
source activate amber_pipeline

echo "CONDA_DEFAULT_ENV=$CONDA_DEFAULT_ENV"
which python3 cpptraj sander
python3 -c "import yaml; print('yaml OK')"

cd /scratch/namyot/BartelleLab/MM-PBSA-Nanobody-Engineering/amber_pipeline
python -m pipeline.run_pipeline config.yaml --step production