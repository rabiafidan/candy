#!/usr/bin/env bash
#SBATCH --partition=cpu                # write partition name(cpu|hmem|int)
#SBATCH --job-name=snakemake                # write your job name
#SBATCH --ntasks=1                   # number of task
#SBATCH --cpus-per-task=1        # Number of CPU cores per task
#SBATCH --nodes=1
#SBATCH --mem=2G
#SBATCH --time 3-00:00:00  #ask this to be extended
#SBATCH --output=%j.txt
#SBATCH --error=%j.err

source ~/.bashrc
conda activate snakemake
snakemake --latency-wait 300 --rerun-incomplete -p --cluster "sbatch --ntasks 1 --cpus-per-task {threads} --partition cpu --job-name {rule} --time 5:00:00 -e {params.err} -o {params.out} --mem {resources.mem_mb}" --jobs 500 --keep-going #--use-singularity
#snakemake --latency-wait 300 --rerun-incomplete -p --cluster "sbatch --ntasks=1 --cpus-per-task={threads} --partition=cpu --job-name={rule} --time=5:00:00 -e={params.err} -o={params.out} --mem={resources.mem_mb}" --jobs 500 --keep-going 