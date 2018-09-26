#!/usr/bin/bash
#
#SBATCH --job-name=cluster_call

#SBATCH --output=.slurm_%j.out

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=1

#SBATCH --time=00:30:00

#SBATCH --mem=400

Rscript ../exec/clustmut.R -i data/ \
                --glob "*.tsv" \
                --recursive \
                --mode distance \
                -o mela_au_sFDR_500 \
                -N 1 \
                --fdr_method "FDR" \
                -d 500 \
                -Vvu 
