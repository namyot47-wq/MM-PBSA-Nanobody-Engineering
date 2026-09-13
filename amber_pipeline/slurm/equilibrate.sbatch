#!/bin/bash
#SBATCH -p public
#SBATCH -q public
#SBATCH -G a100:1
#SBATCH -c 4
#SBATCH -t 0-4
#SBATCH -J nb_equil
#SBATCH -o slurm_logs/equil_%j.out
#SBATCH -e slurm_logs/equil_%j.err

module load amber/24
source activate amber_pipeline

cd /scratch/namyot/BartelleLab/MM-PBSA-Nanobody-Engineering/amber_pipeline
python -m pipeline.run_pipeline config.yaml --step equil
