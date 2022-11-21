#!/usr/bin/env bash
#SBATCH --partition=cpu                # write partition name(cpu|hmem|int)
#SBATCH --job-name=sarek-test                # write your job name
#SBATCH --ntasks=1                   # number of task
#SBATCH --cpus-per-task=1        # Number of CPU cores per task
#SBATCH --nodes=1
#SBATCH --time 3-00:00:00
#SBATCH --output=%j.txt
#SBATCH --error=%j.err

ml purge
ml Nextflow/21.04.0-Java-11
ml Singularity/3.4.2
module load Graphviz/2.38.0-foss-2016b

export NXF_SINGULARITY_CACHEDIR=/camp/project/tracerX/working/PIPELINES/nf-core/sarek/2.7/singularity-images/

nextflow run /camp/project/tracerX/working/PIPELINES/nf-core/sarek/2.7/workflow \
    -profile crick \
    --step mapping \
    --genome GRCh38 \
    --input /nemo/project/proj-tracerX/working/SRM/RABIA/samplesheet.tsv \
    --target_bed /nemo/project/proj-tracerX/working/PIPELINES/nf-core/sarek/2.7/custom_references/bed_files/Panel_v6_hg38.bed \
    --outdir /nemo/project/proj-tracerX/working/SRM/RABIA \
    -c /camp/project/tracerX/working/PIPELINES/nf-core/sarek/2.7/custom_configs/custom_crick_sarek_hp.config \
    -resume
    # -c /nemo/project/proj-tracerX/working/SRM/RABIA/nf.config



#  no target bed needed in WGS
##SBATCH --memory=0
# 1 node has max 32 cpus