#!/bin/bash
#SBATCH -p public
#SBATCH -q public
#SBATCH -G a100:1
#SBATCH -c 4
#SBATCH -t 0-4
#SBATCH -J nb_equil
#SBATCH -o slurm_logs/equil_%j.out
#SBATCH -e slurm_logs/equil_%j.err

set -euo pipefail
mkdir -p /scratch/namyot/BartelleLab/MM-PBSA-Nanobody-Engineering/amber_pipeline/slurm_logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate amber_pipeline
export PATH="/home/namyot/.conda/envs/amber_pipeline/bin:$PATH"
unset PYTHONPATH

echo "CONDA_DEFAULT_ENV=$CONDA_DEFAULT_ENV"
which python3
python3 -c "import sys; print(sys.executable)"
python3 -c "import yaml; print('yaml OK')"
which cpptraj
which sander

cd /scratch/namyot/BartelleLab/MM-PBSA-Nanobody-Engineering/amber_pipeline
python -m pipeline.run_pipeline config.yaml --step equil
