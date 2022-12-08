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


##### uncomment as you wish ######

# single sample cluster
#snakemake Mutect2/temp/multi/RFH001_norm.vcf.gz --latency-wait 100 --rerun-incomplete -p \
#--cluster "sbatch --ntasks 1 --cpus-per-task {threads} --partition cpu --job-name {rule} --time 5:00:00 -e logs/{rule}/{params.err} -o logs/{rule}/{params.out} --mem {resources.mem_mb} --parsable" \
#--jobs 10 --keep-going --use-envmodules --cluster-status ./status-sacct.sh

#single sample local
#snakemake Mutect2/multi/RFH001.vcf.gz --latency-wait 100 --rerun-incomplete -p --use-envmodules -c1

#all cluster
snakemake -n --latency-wait 100 --rerun-incomplete -p \
--cluster "sbatch --ntasks 1 --cpus-per-task {threads} --partition cpu --job-name {rule} --time 5:00:00 -e logs/{rule}/{params.err} -o logs/{rule}/{params.out} --mem {resources.mem_mb} --parsable" \
--jobs 500 --keep-going --use-envmodules --cluster-status ./status-sacct.sh

#singe samle dryrun
#snakemake Mutect2/multi/RFH001.vcf.gz --rerun-incomplete -p --use-envmodules --dry-run

#all dryrun
#snakemake --rerun-incomplete -p --use-envmodules --dry-run