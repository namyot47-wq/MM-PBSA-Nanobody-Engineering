#!/bin/bash
#SBATCH -p public
#SBATCH -q public
#SBATCH -c 4
#SBATCH -t 0-8
#SBATCH -J nb_mmpbsa
#SBATCH -o slurm_logs/mmpbsa_%j.out
#SBATCH -e slurm_logs/mmpbsa_%j.err

module load mamba/latest
module load amber/24
source activate amber_pipeline

cd /scratch/namyot/BartelleLab/MM-PBSA-Nanobody-Engineering/amber_pipeline
python -m pipeline.run_pipeline config.yaml --step mmpbsa