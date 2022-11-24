configfile: "config.yaml"
include: 'params.snake'
import pandas as pd
import sys
import logging


recal=pd.read_csv(RECAL,sep='\t',names=['PatientName','Sex','SampleType','SampleName','bam','bai'])
recal["rep"]=recal.apply(lambda x: "multi" if len(recal.loc[(recal['PatientName']==x.PatientName) & (recal['SampleType']==1),:])>1 else "single" , axis=1)
ALL_TUMOURS=[x for x in recal.loc[recal['SampleType']==1,'SampleName']]
all_single_tumour=[x for x in recal.loc[(recal['SampleType']==1)&(recal['rep']=='single'),'SampleName']]
all_multi_tumour=[x for x in recal.loc[(recal['SampleType']==1)&(recal['rep']=='multi'),'SampleName']]
patients=list(set(recal['PatientName']))
multi_patients=list(set(recal.loc[recal['rep']=='multi',:]['PatientName']))
single_patients=list(set(recal.loc[recal['rep']=='single',:]['PatientName']))

normals=[]
tumours=[]
nb=[]
tb=[]
for i in patients:
    tab=recal.loc[(recal['PatientName']==i) & (recal['SampleType']==0) ,['SampleName',"bam"]] #germline
    tab2=recal.loc[(recal['PatientName']==i) & (recal['SampleType']==1) ,['SampleName',"bam"]] #tumour
    tumours.append([x for x in tab2["SampleName"]])
    normals.append([x for x in tab["SampleName"]])
    tb.append([x for x in tab2["bam"]])
    nb.append([x for x in tab["bam"]])

def get_normal_name(wildcards):
    pat=recal.loc[recal['SampleName']==wildcards.tum,'PatientName']
    nor=recal.loc[(recal['PatientName']==pat.iloc[0]) & (recal['SampleType']==0) ,'SampleName']
    return nor.iloc[0]
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


onstart:
    if set([len(x) for x in normals]) != {1}:
        print("Each patient must have only 1 normal!")
        sys.exit(1)
    shell('mkdir -p logs/Mutect2_call logs/learn_model logs/add_flags logs/vcf_index logs/vcf_normalise logs/qual_filter logs/patient_positions')

localrules: all 

rule all:
    input:
        #['Mutect2/temp/single/fil_{tum}.vcf.gz'.format(tum=tums)for tums in all_single_tumour],
        #['Mutect2/temp/multi/fil_{tum}.vcf.gz'.format(tum=tums)for tums in all_multi_tumour],
        ['Mutect2/temp/multi/{patient}_valid_pos.vcf'.format(patient=pats)for pats in multi_patients],
        ['Mutect2/temp/single/fil_{tum}.vcf.gz.tbi'.format(tum=tums)for tums in all_single_tumour]
        

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
        temp("Mutect2/unfiltered/{tum}.vcf.gz.stats"),
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
        'Mutect2/temp/or_flagged_{tum}.vcf.gz.filteringStats.tsv',
        v='Mutect2/temp/or_flagged_{tum}.vcf.gz'
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
        'Mutect2/temp/or_flagged_{tum}.vcf.gz'
    output:
        'Mutect2/temp/or_flagged_{tum}.vcf.gz.tbi'
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
        'Mutect2/temp/or_flagged_{tum}.vcf.gz.tbi',
        v='Mutect2/temp/or_flagged_{tum}.vcf.gz'
    output:
        v='Mutect2/temp/norm_{tum}.vcf.gz',
        i='Mutect2/temp/norm_{tum}.vcf.gz.tbi'
    params:
        err=lambda wildcards: "logs/{rule}/{wildcards.tum}.err",
	    out=lambda wildcards: "logs/{rule}/{wildcards.tum}.out"	
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
        v='Mutect2/temp/{rep}/fil_{tum}.vcf.gz',
        i='Mutect2/temp/{rep}/fil_{tum}.vcf.gz.tbi'
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
        -i '(ROQ>{params.ROQ} & AD[1:1]>{params.AD}) & (FILTER="PASS" | FILTER="clustered_events" )' -Oz -o {output.v} {input.v}
        /nemo/lab/turajlics/home/users/fidanr/bcftools/bin/bcftools tabix {output.v}
        """

### MULTI SAMPLE ###
rule patient_positions:
    """
    Union of all valid postions (allele aware) of samples coming from a patient.
    Step #3 of filterAndMergeSingleSample_noGLrescue.sh
    """
    input:
       v=lambda wildcards: ["Mutect2/temp/multi/fil_"+tum+".vcf.gz" for tum in tumours[patients.index(wildcards.patient)]],
       i=lambda wildcards: ["Mutect2/temp/multi/fil_"+tum+".vcf.gz.tbi" for tum in tumours[patients.index(wildcards.patient)]]
    output:
       "Mutect2/temp/multi/{patient}_valid_pos.txt"
    params:
        err=lambda wildcards: "logs/{rule}/{wildcards.patient}.err",
	    out=lambda wildcards: "logs/{rule}/{wildcards.patient}.out"
    threads: 1
    resources:
        mem_mb=10000
    shell:
        """
        zcat {input.v} | grep -v "^#" | cut -f1-5 | sort -V | uniq > {output}
        """
   
rule retrieve_from_vcf:
    """
    Filter unfiltered VCF with these positions
    Step #4 of filterAndMergeSingleSample_noGLrescue.sh
    """
    input:
        pos=lambda wildcards: "Mutect2/temp/multi/" + patients[tumours.index(wildcards.tum)]+ "_valid_pos.txt",
        v="Mutect2/temp/norm_{tum}.vcf.gz",
        i="Mutect2/temp/norm_{tum}.vcf.gz.tbi"
    output:
       v="Mutect2/temp/multi/{tum}_premerge.vcf.gz",
       i="Mutect2/temp/multi/{tum}_premerge.vcf.gz.tbi"
    params:
        err=lambda wildcards: "logs/{rule}/{wildcards.tum}.err",
	    out=lambda wildcards: "logs/{rule}/{wildcards.tum}.out"
    threads: 1
    resources:
        mem_mb=10000
    envmodules:
        "HTSlib/1.14-GCC-11.2.0"
    shell:
        """
	    zcat {input.v} | grep "^#" | bgzip >{output.v}
	    for x in $(seq 1 1000 $(wc -l {input.pos} | awk '{{print $1}}'))
	    do
		y=$((x+999))
		zcat {input.v} | grep -v "^#" | grep --line-buffered -f <(sed -n "${{x}},${{y}}p" {input.pos}) - | bgzip >> {output.v}
	    done
	    /nemo/lab/turajlics/home/users/fidanr/bcftools/bin/bcftools tabix {output.v}
        """


rule merge_vcf:
    """
    Merge filtered vcf files
    Step #5 of filterAndMergeSingleSample_noGLrescue.sh
    """
    input:
       v=lambda wildcards: ["Mutect2/temp/multi/"+tum+"_premerge.vcf.gz" for tum in tumours[patients.index(wildcards.patient)]],
       i=lambda wildcards: ["Mutect2/temp/multi/"+tum+"_premerge.vcf.gz.tbi" for tum in tumours[patients.index(wildcards.patient)]]
    output:
       "Mutect2/temp/multi/{patient}_merged.vcf.gz.tbi",
       v="Mutect2/temp/multi/{patient}_merged.vcf.gz"
    params:
        err=lambda wildcards: "logs/{rule}/{wildcards.patient}.err",
	    out=lambda wildcards: "logs/{rule}/{wildcards.patient}.out"
    threads: 1
    resources:
        mem_mb=10000
    envmodules:
        "HTSlib/1.14-GCC-11.2.0"
    shell:
        """
        /nemo/lab/turajlics/home/users/fidanr/bcftools/bin/bcftools merge -0 -Oz -o {output.v} --force-samples {input.v}
	    /nemo/lab/turajlics/home/users/fidanr/bcftools/bin/bcftools tabix {output.v}
        """


rule remove_germline:
    """
    Remove extra germline columns in vcfs
    Step #6 of filterAndMergeSingleSample_noGLrescue.sh
    """
    input:
       "Mutect2/temp/multi/{patient}_merged.vcf.gz.tbi",
       v="Mutect2/temp/multi/{patient}_merged.vcf.gz"
    output:
       "Mutect2/temp/multi/{patient}_merged_GL_removed_sorted.vcf.gz.tbi",
       s=temp("{patient}_to_keep.txt"),
       v1=temp("Mutect2/temp/multi/{patient}_merged_GL_removed.vcf.gz"),
       v2="Mutect2/temp/multi/{patient}_merged_GL_removed_sorted.vcf.gz"
    params:
        err=lambda wildcards: "logs/{rule}/{wildcards.patient}.err",
	    out=lambda wildcards: "logs/{rule}/{wildcards.patient}.out",
        GL=lambda wildcards: normals[patients.index(wildcards.patient)][0]
    threads: 1      
    resources:
        mem_mb=10000
    envmodules:
        "HTSlib/1.14-GCC-11.2.0"
    shell:
        """
        echo {params.GL} > {output.s}
        /nemo/lab/turajlics/home/users/fidanr/bcftools/bin/bcftools query -l {input.v} | grep -v {params.GL} >> {output.s}
        /nemo/lab/turajlics/home/users/fidanr/bcftools/bin/bcftools view -S {output.s} {input.v} -Oz -o {output.v1} 
        /nemo/lab/turajlics/home/users/fidanr/bcftools/bin/bcftools sort {output.v1} -Oz -o {output.v2}
        tabix -p vcf {output.v2}
        """


rule normalise_merged:
    """
    Normalise the merged vcfs
    """
    input:
       "Mutect2/temp/multi/{patient}_merged_GL_removed_sorted.vcf.gz.tbi",
       v="Mutect2/temp/multi/{patient}_merged_GL_removed_sorted.vcf.gz"
    output:
       "Mutect2/temp/multi/{patient}_norm.vcf.gz"
    params:
        err=lambda wildcards: "logs/{rule}/{wildcards.patient}.err",
	    out=lambda wildcards: "logs/{rule}/{wildcards.patient}.out"
    threads: 1      
    resources:
        mem_mb=10000
    envmodules:
        "HTSlib/1.14-GCC-11.2.0"  
    shell:
        """
        /nemo/lab/turajlics/home/users/fidanr/bcftools/bin/bcftools norm -m -any {input.v} | bgzip > {output}
        """


rule rescue_mission:
    """
    Read the missing information vcf files from corresponding bam files and fill in that information in vcfs.
    """
    input:
        v="Mutect2/temp/multi/{patient}_norm.vcf.gz",
        nb=lambda wildcards: nb[patients.index(wildcards.patient)][0],
        tb=lambda wildcards: tb[patients.index(wildcards.patient)]
    output:
       v="Mutect2/multi/{patient}.vcf.gz",
       bamlist=temp("Mutect2/temp/multi/{patient}_bamlist.txt"),
       obamlist=temp("Mutect2/temp/multi/{patient}_ordered_bamlist.txt"),
       s=temp("Mutect2/temp/multi/{patient}_samples.txt")
    params:
        err=lambda wildcards: "logs/{rule}/{wildcards.patient}.err",
	    out=lambda wildcards: "logs/{rule}/{wildcards.patient}.out",
        res=RES
    threads: 1      
    resources:
        mem_mb=10000
    envmodules:
        "Python/3.6.6-foss",
        "Pysam/0.15.1-foss-2018b-Python-3.6.6"    
    shell:
        """
        # Make a list of bams in VCF order
        echo {input.nb} > {output.bamlist}
        for i in {input.tb}; do echo $i >> {output.bamlist}; done
        /nemo/lab/turajlics/home/users/fidanr/bcftools/bin/bcftools query -l {input.v} >> {output.s}
        while read i; do grep $i {output.bamlist} >> {output.obamlist}; done< {output.s}
        #Genotype (SNVs, MNVs and indels)
        python3.6 {params.res} {input.v} {output.obamlist} {output.v}
        """