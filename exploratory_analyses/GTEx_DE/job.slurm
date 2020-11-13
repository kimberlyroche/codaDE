#!/bin/bash
#SBATCH -J DE_GTEx
#SBATCH --mem=24GB
#SBATCH --get-user-env
#SBATCH --time=01:00:00
#

## Submit as: sbatch --array=1-50 job.slurm

module add R/3.6.1-gcb03
module add gcc/7.1.0-fasrc01

cd /data/mukherjeelab/roche/codaDE/exploratory/GTEx_DE
srun Rscript run.R
