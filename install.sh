#!/bin/bash
#SBATCH --job-name='install'
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=log.ins.log
#SBATCH --mem=20G
#SBATCH -t 14-00:00:00

conda env create -n SpecImmune -f environment.yml
