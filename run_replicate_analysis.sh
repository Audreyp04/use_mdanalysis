#!/bin/bash
#SBATCH -n 1
#SBATCH -A bevanlab
#SBATCH -p normal_q
#SBATCH -t 6:00:00
#SBATCH --job-name=hex1-analysis
#SBATCH --mail-user=audreyp04@vt.edu
#SBATCH --mail-type=all
#SBATCH -o %x-slurm-%j.out

python run_analysis.py