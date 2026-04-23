#!/bin/bash
#SBATCH -n 1
#SBATCH -A bevanlab
#SBATCH -p normal_q
#SBATCH -t 2:00:00
#SBATCH --job-name=mindist
#SBATCH --mail-user=audreyp04@vt.edu
#SBATCH --mail-type=all
#SBATCH -o %x-slurm-%j.out

module load gromacs/2025.2

echo 1 5|gmx mindist -f ../traj/cat_pbc_nowat.xtc -s nowat.tpr -n ../index.ndx -dt 10