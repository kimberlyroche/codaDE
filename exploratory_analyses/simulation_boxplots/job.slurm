#!/bin/bash
#SBATCH -J condition1
#SBATCH --mem=16GB
#SBATCH --get-user-env
#SBATCH --time=02:00:00
#

## Submit as: sbatch --array=1-20 job.slurm

module add R/3.6.1-gcb03
module add gcc/7.1.0-fasrc01

cd /data/mukherjeelab/roche/codaDE/exploratory/simulation_boxplots
srun Rscript generate_data.R --cond=1
