#!/usr/bin/env bash
#SBATCH --partition=cpu                # write partition name(cpu|hmem|int)
#SBATCH --job-name=sarek1                # write your job name
#SBATCH --ntasks=1                   # number of task
#SBATCH --cpus-per-task=1        # Number of CPU cores per task
#SBATCH --nodes=1
#SBATCH --time 3-00:00:00
#SBATCH --output=%j.txt
#SBATCH --error=%j.err

ml purge
ml Nextflow/22.04.0
ml Singularity/3.4.2

export NXF_SINGULARITY_CACHEDIR=/camp/project/tracerX/working/PIPELINES/nf-core/sarek/2.7/singularity-images/

##########   ALIGNMENT  ###############
#nextflow run /camp/project/tracerX/working/PIPELINES/nf-core/sarek/2.7/workflow \
#    -profile crick \
#    --step mapping \
#    --genome GRCh38 \
#    --input /nemo/project/proj-tracerX/working/SRM/RABIA/SAREK_1/samples_1.tsv \
#    --target_bed /nemo/project/proj-tracerX/working/PIPELINES/nf-core/sarek/2.7/custom_references/bed_files/Panel_v6_hg38.bed \
#    --outdir /nemo/project/proj-tracerX/working/SRM/RABIA/SAREK_1 \
#    -c /camp/project/tracerX/working/PIPELINES/nf-core/sarek/2.7/custom_configs/custom_crick_sarek_hp.config \
#    -resume \
#    -dsl1
#    # -c /nemo/project/proj-tracerX/working/SRM/RABIA/nf.config


###### ASCAT###########
nextflow run /camp/project/tracerX/working/PIPELINES/nf-core/sarek/2.7/workflow \
    -profile crick \
    --step variant_calling \
    --tools ASCAT \
    --genome GRCh38 \
    --input /nemo/project/proj-tracerX/working/SRM/RABIA/SAREK_1/Preprocessing/TSV/recalibrated.tsv \
    --target_bed /nemo/project/proj-tracerX/working/PIPELINES/nf-core/sarek/2.7/custom_references/bed_files/Panel_v6_hg38.bed \
    --outdir /nemo/project/proj-tracerX/working/SRM/RABIA/SAREK_1 \
    -c /camp/project/tracerX/working/PIPELINES/nf-core/sarek/2.7/custom_configs/custom_crick_sarek_hp.config \
    -resume \
    -dsl1


#  no target bed needed in WGS
##SBATCH --memory=0
# 1 node has max 32 cpus
