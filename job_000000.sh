#!/bin/bash

#SBATCH --job-name=pc_PYMC000000
#SBATCH -p general
#SBATCH -o pc_PYMC_%j.txt
#SBATCH -e pc_PYMC_%j.err
#SBATCH --nodes=None
#SBATCH --ntasks-per-node=None
#SBATCH --cpus-per-task=None
#SBATCH --time=0-00:00:00
