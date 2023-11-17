#!/bin/bash
#
#SBATCH --job-name=gNOMO2           # Job name, will show up in squeue output
#SBATCH --ntasks=<cores>            # Number of cores
#SBATCH --time=0-48:00:00           # Runtime in DAYS-HH:MM:SS format
#SBATCH --mem=250000                # RAM 
#SBATCH --nodes=<nodes>             # Number of nodes
#SBATCH --output=gNOMO2_%j.out      # File to which standard out will be written
#SBATCH --error=gNOMO2_%j.err       # File to which standard err will be written
#SBATCH --mail-type=ALL             # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=<email>         # Email to which notifications will be sent SBATCH --qos=medium

cd /path/to/gNOMO2
source $(conda info --root)/etc/profile.d/conda.sh
conda activate snakemake
snakemake --cores <cores> --use-conda 
conda deactivate
