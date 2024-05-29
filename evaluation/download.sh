#!/bin/bash
#SBATCH --job-name=clr
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=20-10:10:00
#SBATCH --output=download.log


python3 get_real_data.py
