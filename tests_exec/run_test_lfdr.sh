#!/usr/bin/bash

#SBATCH --job-name=TEST_clustmut

#SBATCH --output=.slurm_%j.out

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=1

#SBATCH --time=00:60:00

#SBATCH --mem=10000

Rscript ../exec/clustmut_distance.R -i data/ \
                --glob "*.tsv" \
                --recursive \
                -o mela_au_fdr \
                -N 1 \
                -Vvu 
