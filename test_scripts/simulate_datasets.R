#!/bin/bash
#SBATCH -J sim_data
#SBATCH --mem=16GB
#SBATCH --get-user-env
#SBATCH --time=12:00:00
#

module add R/3.6.1-gcb03
module add gcc/7.1.0-fasrc01

cd /data/mukherjeelab/roche/codaDE

srun Rscript simulate_datasets.R
