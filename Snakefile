configfile: "config.yaml"
include: 'params.snake'
import pandas as pd
import sys
import logging

logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.StreamHandler()])

recal=pd.read_csv(RECAL,sep='\t',names=['PatientName','Sex','SampleType','SampleName','bam','bai'])
ALL_TUMOURS=[x for x in recal.loc[recal['SampleType']==1,'SampleName']]
patients=list(set(recal['PatientName']))
recal["rep"]=recal.apply(lambda x: "multi" if len(recal.loc[(recal['PatientName']==x.PatientName) & (recal['SampleType']==1),:])>1 else "single" , axis=1)

normals=[]
tumours=[]
for i in patients:
    tab=recal.loc[(recal['PatientName']==i) & (recal['SampleType']==0) ,'SampleName'] #germline
    tab2=recal.loc[(recal['PatientName']==i) & (recal['SampleType']==1) ,'SampleName'] #tumour
    tumours.append([x for x in tab2])
    normals.append([x for x in tab])

localrules: all

def get_normal_name(wildcards):
    pat=recal.loc[recal['SampleName']==wildcards.tum,'PatientName']
    nor=recal.loc[(recal['PatientName']==pat.iloc[0]) & (recal['SampleType']==0) ,'SampleName']
    return nor[0]
def get_normal_bam(wildcards):
    pat=recal.loc[recal['SampleName']==wildcards.tum,'PatientName']
    norbam=recal.loc[(recal['PatientName']==pat.iloc[0]) & (recal['SampleType']==0) ,'bam']
    return list(norbam)
def get_normal_bai(wildcards):
    pat=recal.loc[recal['SampleName']==wildcards.tum,'PatientName']
    norbai=recal.loc[(recal['PatientName']==pat.iloc[0]) & (recal['SampleType']==0) ,'bai']
    return list(norbai)
def get_tumour_bam(wildcards):
    tumbam=recal.loc[recal['SampleName']==wildcards.tum ,'bam']
    return list(tumbam)
def get_tumour_bai(wildcards):
    tumbai=recal.loc[recal['SampleName']==wildcards.tum ,'bai']
    return list(tumbai)

def get_patient_tumour_vcf(wildcards):
    idx=patients.index(wildcards.patient)
    pat_tum= tumours[idx]
    return [f"Mutect2/temp/fil_{a}.vcf.gz" for a in pat_tum]

onstart:
    if set([len(x)] for x in normals) != {1}:
        logging.error("Each patient must have only 1 normal!")
        sys.exit(1)

rule all:
    input:
        "Mutect2/temp/test1_valid_positions"
        #['Mutect2/or_filtered/{tum}.vcf.gz'.format(tum=tums)for tums in ALL_TUMOURS]
        #['Mutect2/unfiltered/{tum}.vcf.gz'.format(tum=tums)for tums in ALL_TUMOURS]

rule Mutect2_call:
    """
    Call variants and count reads for orientation filtering
    #Step2 of mutect2_ROB_isLegacy_targetedOrWGS_hg38.sh
    """
    input:
        get_normal_bai,
        get_tumour_bai,
        tb=get_tumour_bam,
        nb=get_normal_bam,
        ref=REF_GEN,
        gr=GR,
        intv=INT,
        pon=PON
    output:
        v='Mutect2/unfiltered/{tum}.vcf.gz',
        n='Mutect2/noise/{tum}.tar.gz'
    params:
        nn=get_normal_name,
        err=lambda wildcards: "logs/{rule}/{wildcards.tum}.err",
	    out=lambda wildcards: "logs/{rule}/{wildcards.tum}.out"	
    envmodules :
        "GATK/4.1.8.1-GCCcore-9.3.0-Java-1.8"
    threads: 8
    resources:
        mem_mb=16000
    shell:
        "java -jar $EBROOTGATK/gatk-package-4.1.8.1-local.jar Mutect2 \
        --native-pair-hmm-threads {threads} \
        -R {input.ref} \
        -L {input.intv} \
        -I {input.tb} \
        -I {input.nb} \
        -normal {params.nn} \
        -germline-resource {input.gr} \
        -pon {input.pon} \
        --f1r2-tar-gz {output.n} \
        -O {output.v}"

rule learn_model:
    """
    Infer noise priors
    #Step3 of mutect2_ROB_isLegacy_targetedOrWGS_hg38.sh
    """
    input:
        'Mutect2/noise/{tum}.tar.gz'
    output:
        'Mutect2/noise/{tum}-read-orientation-model.tar.gz'
    params:
        err=lambda wildcards: "logs/{rule}/{wildcards.tum}.err",
	    out=lambda wildcards: "logs/{rule}/{wildcards.tum}.out"	
    envmodules :
        "GATK/4.1.8.1-GCCcore-9.3.0-Java-1.8"
    threads: 1
    resources:
        mem_mb=4000
    shell:
        "java -jar $EBROOTGATK/gatk-package-4.1.8.1-local.jar LearnReadOrientationModel \
        -I {input} \
        -O {output}"

   
rule add_flags:
    """
    Add filtering flags to vcf.
    #Step4 of mutect2_ROB_isLegacy_targetedOrWGS_hg38.sh
    """
    input:
        r=REF_GEN,
        v='Mutect2/unfiltered/{tum}.vcf.gz',
        priors='Mutect2/noise/{tum}-read-orientation-model.tar.gz'
    output:
        v='Mutect2/or_filtered/{tum}.vcf.gz'
    params:
        err=lambda wildcards: "logs/{rule}/{wildcards.tum}.err",
	    out=lambda wildcards: "logs/{rule}/{wildcards.tum}.out"	
    envmodules :
        "GATK/4.1.8.1-GCCcore-9.3.0-Java-1.8"
    threads: 1
    resources:
        mem_mb=4000
    shell:
        """
        java -jar $EBROOTGATK/gatk-package-4.1.8.1-local.jar FilterMutectCalls \
        -R {input.r} \
        -V {input.v} \
        --ob-priors {input.priors} \
        -O {output.v}
        """

rule vcf_index:
    """
    Index vcf
    #Step4 of mutect2_ROB_isLegacy_targetedOrWGS_hg38.sh
    """
    input:
        'Mutect2/or_filtered/{tum}.vcf.gz'
    output:
        'Mutect2/or_filtered/{tum}.vcf.gz.tbi'
    params:
        err=lambda wildcards: "logs/{rule}/{wildcards.tum}.err",
	    out=lambda wildcards: "logs/{rule}/{wildcards.tum}.out"	
    envmodules :
        "BCFtools/1.12-GCC-10.2.0"
    threads: 1
    resources:
        mem_mb=4000
    shell:
        "tabix -p vcf {input}"


rule vcf_normalise:
    """
    Normalize vcf, fix a bug in vcf tag, and index the vcf
    #Step1 of filterAndMergeSingleSample_noGLrescue.sh
    """
    input:
        'Mutect2/or_filtered/{tum}.vcf.gz.tbi',
        v='Mutect2/or_filtered/{tum}.vcf.gz'
    output:
        v=temp('Mutect2/temp/norm_{tum}.vcf.gz'),
        i=temp('Mutect2/temp/norm_{tum}.vcf.gz.tbi')
    params:
        err=lambda wildcards: "logs/{rule}/{wildcards.tum}.err",
	    out=lambda wildcards: "logs/{rule}/{wildcards.tum}.out"	
    envmodules :
        "BCFtools/1.12-GCC-10.2.0"
    threads: 1
    resources:
        mem_mb=10000
    conda:
        "snakemake"
    shell:
        """
        #echo "##INFO=<ID=AS_FilterStatus,Number=A" | sed "s/##INFO=<ID=AS_FilterStatus,Number=A/##INFO=<ID=AS_FilterStatus,Number=1/" && touch {output.v} && touch {output.i}
        /nemo/lab/turajlics/home/users/fidanr/bcftools/bin/bcftools norm -m -any {input.v} | sed "s/##INFO=<ID=AS_FilterStatus,Number=A/##INFO=<ID=AS_FilterStatus,Number=1/" | pbgzip -c -t 0 > {output.v}
	    /nemo/lab/turajlics/home/users/fidanr/bcftools/bin/bcftools tabix {output.v}
        """

rule qual_filter:
    """
    Filter vcf based on quality metrics
    #Step2 of filterAndMergeSingleSample_noGLrescue.sh
    #Variants are kept if AD> <threshold> & ROQ>= <threshold> & FILTER is PASS or clustered_events

    For now, it uses a local program path. We can add it to the conda env when v1.16 is installable by anaconda
    """
    input:
        'Mutect2/temp/norm_{tum}.vcf.gz.tbi',
        v='Mutect2/temp/norm_{tum}.vcf.gz'
    output:
        v=temp('Mutect2/temp/fil_{tum}.vcf.gz'),
        i=temp('Mutect2/temp/fil_{tum}.vcf.gz.tbi')
    params:
        AD=AD,
        ROQ=ROQ,
        err=lambda wildcards: "logs/{rule}/{wildcards.tum}.err",
	    out=lambda wildcards: "logs/{rule}/{wildcards.tum}.out"	
    envmodules :
        "BCFtools/1.12-GCC-10.2.0"
    threads: 1
    resources:
        mem_mb=10000
    shell:
        """
        /nemo/lab/turajlics/home/users/fidanr/bcftools/bin/bcftools view \
        -i '(ROQ>{params.ROQ} & AD[0:1]>{params.AD}) & (FILTER="PASS" | FILTER="clustered_events" )' -Oz -o {output.v} {input.v}
        /nemo/lab/turajlics/home/users/fidanr/bcftools/bin/bcftools tabix {output.v}
        """


rule patient_positions:
    """
    Union of all valid postions (allele aware) of samples coming from a patient.
    Step #4 of filterAndMergeSingleSample_noGLrescue.sh
    """
    input:
        expand('Mutect2/temp/fil_{tum}.vcf.gz', tum=get_patient_tumours())
    output:
        "Mutect2/temp/{patient}_valid_positions"
    params:
        err=lambda wildcards: "logs/{rule}/{wildcards.patient}.err",
	    out=lambda wildcards: "logs/{rule}/{wildcards.patient}.out"
    threads: 1
    resources:
        mem_mb=10000
    shell:
        """
        zcat {input} | grep -v "^#" | cut -f1-5 | sort -V | uniq > {output}
        echo $(wc -l {output} | awk '{print $1}') "valid variants found"
        """