configfile: "config.yaml"
include: 'params.snake'
import pandas as pd
import sys
import logging


onstart:
    if set([len(x) for x in normals]) != {1}:
        print("Each patient must have only 1 normal!")
        sys.exit(1)
    shell('mkdir -p logs/germline_SNP logs/VAF_filter logs/fill_tags_VAF logs/vcf2maf logs/Mutect2_call logs/learn_model logs/add_flags logs/vcf_index logs/vcf_normalise logs/qual_filter logs/patient_positions logs/single_sample logs/retrieve_from_vcf logs/merge_vcf logs/remove_germline logs/normalise_merged logs/rescue_mission')

localrules: all 

rule all:
    input:
        #['Mutect2/temp/single/fil_{tum}.vcf.gz'.format(tum=tums)for tums in all_single_tumour],
        ['Mutect2/temp/multi/fil_{tum}.vcf.gz'.format(tum=tums)for tums in all_multi_tumour],
        ['Mutect2/multi/{patient}.vcf.gz'.format(patient=pats)for pats in multi_patients],
        ['Mutect2/single/{patient}.vcf.gz'.format(patient=pats)for pats in single_patients],
        ['GL_variant/'+x+'.vcf' for x in recal["SampleName"]]
        #["MAF/multi/{tum}.maf".format(tum=tums)for tums in all_multi_tumour],
        #["MAF/single/{tum}.maf".format(tum=tums)for tums in all_single_tumour]

        

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
        "Mutect2/unfiltered/{tum}.vcf.gz.stats",
        v='Mutect2/unfiltered/{tum}.vcf.gz',
        n='Mutect2/noise/{tum}.tar.gz'
    params:
        nn=get_normal_name,
        err=lambda wildcards: wildcards.tum+".err",
        out=lambda wildcards: wildcards.tum+".out"    
    envmodules :
        "GATK/4.1.8.1-GCCcore-9.3.0-Java-1.8"
    threads: 8
    resources:
        mem_mb=16000
    shell:
        """
        java -jar $EBROOTGATK/gatk-package-4.1.8.1-local.jar Mutect2 \
        --native-pair-hmm-threads {threads} \
        -R {input.ref} \
        -L {input.intv} \
        -I {input.tb} \
        -I {input.nb} \
        -normal {params.nn} \
        -germline-resource {input.gr} \
        -pon {input.pon} \
        --f1r2-tar-gz {output.n} \
        -O {output.v}
        """

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
        err=lambda wildcards: wildcards.tum+".err",
        out=lambda wildcards: wildcards.tum+".out"    
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
        "Mutect2/unfiltered/{tum}.vcf.gz.stats",
        r=REF_GEN,
        v='Mutect2/unfiltered/{tum}.vcf.gz',
        priors='Mutect2/noise/{tum}-read-orientation-model.tar.gz'
    output:
        'Mutect2/temp/or_flagged_{tum}.vcf.gz.filteringStats.tsv',
        v='Mutect2/temp/or_flagged_{tum}.vcf.gz'
    params:
        err=lambda wildcards: wildcards.tum+".err",
        out=lambda wildcards: wildcards.tum+".out"    
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
        err=lambda wildcards: wildcards.tum+".err",
        out=lambda wildcards: wildcards.tum+".out"    
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
        err=lambda wildcards: wildcards.tum+ ".err",
        out=lambda wildcards: wildcards.tum+ ".out"    
    threads: 1
    resources:
        mem_mb=10000
    envmodules:
        "HTSlib/1.15.1-GCC-11.3.0" 
    shell:
        """
        /nemo/lab/turajlics/home/users/fidanr/bcftools/bin/bcftools norm -m -any {input.v} | sed "s/##INFO=<ID=AS_FilterStatus,Number=A/##INFO=<ID=AS_FilterStatus,Number=1/" | bgzip > {output.v}
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
        err=lambda wildcards: wildcards.tum+ ".err",
        out=lambda wildcards: wildcards.tum+ ".out"    
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
### SINGLE SAMPLE ###
rule single_sample:
    """
    For single sample tumours, this is the final filtered vcf. Copy this to a permanent result diretory.
    """
    input:
        v=lambda wildcards: "Mutect2/temp/single/fil_"+ tumours[patients.index(wildcards.patient)][0]+".vcf.gz",
        i=lambda wildcards: "Mutect2/temp/single/fil_"+ tumours[patients.index(wildcards.patient)][0]+".vcf.gz.tbi"
    output:
        v='Mutect2/single/{patient}.vcf.gz',
        i='Mutect2/single/{patient}.vcf.gz.tbi'
    params:
        err=lambda wildcards: wildcards.patient+ ".err",
        out=lambda wildcards: wildcards.patient+ ".out"    
    threads: 1
    resources:
        mem_mb=10000
    shell:
        """
        cp {input.v} {output.v}
        cp {input.i} {output.i}
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
        err=lambda wildcards: wildcards.patient+ ".err",
        out=lambda wildcards: wildcards.patient+ ".out"
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
        pos=lambda wildcards: "Mutect2/temp/multi/" + recal.loc[recal['SampleName']==wildcards.tum,'PatientName'].iloc[0]+ "_valid_pos.txt",
        v="Mutect2/temp/norm_{tum}.vcf.gz",
        i="Mutect2/temp/norm_{tum}.vcf.gz.tbi"
    output:
       v="Mutect2/temp/multi/{tum}_premerge.vcf.gz",
       i="Mutect2/temp/multi/{tum}_premerge.vcf.gz.tbi"
    params:
        err=lambda wildcards: wildcards.tum+ ".err",
        out=lambda wildcards: wildcards.tum+ ".out"
    threads: 1
    resources:
        mem_mb=10000
    envmodules:
        "HTSlib/1.15.1-GCC-11.3.0"
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
        err=lambda wildcards: wildcards.patient+ ".err",
        out=lambda wildcards: wildcards.patient+ ".out"
    threads: 1
    resources:
        mem_mb=10000
    envmodules:
        "HTSlib/1.15.1-GCC-11.3.0"
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
       s="{patient}_to_keep.txt",
       v1="Mutect2/temp/multi/{patient}_merged_GL_removed.vcf.gz",
       v2="Mutect2/temp/multi/{patient}_merged_GL_removed_sorted.vcf.gz"
    params:
        err=lambda wildcards: wildcards.patient+ ".err",
        out=lambda wildcards: wildcards.patient+ ".out",
        GL=lambda wildcards: normals[patients.index(wildcards.patient)][0]
    threads: 1      
    resources:
        mem_mb=10000
    envmodules:
        "HTSlib/1.15.1-GCC-11.3.0"
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
        err=lambda wildcards: wildcards.patient+".err",
        out=lambda wildcards: wildcards.patient+".out"
    threads: 1      
    resources:
        mem_mb=8000
    envmodules:
        "HTSlib/1.15.1-GCC-11.3.0"  
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
       bamlist="Mutect2/temp/multi/{patient}_bamlist.txt",
       obamlist="Mutect2/temp/multi/{patient}_ordered_bamlist.txt",
       s="Mutect2/temp/multi/{patient}_samples.txt"
    params:
        err=lambda wildcards: wildcards.patient+ ".err",
        out=lambda wildcards: wildcards.patient+ ".out",
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


rule fill_tags_VAF:
    """
    Add VAF tags to the VCF files. 
    """
    input:
        v="Mutect2/multi/{patient}.vcf.gz"
    output:
        v="Mutect2/multi/tag_VAF/{patient}.vcf.gz",
        i="Mutect2/multi/tag_VAF/{patient}.vcf.gz.tbi"
    params:
        err=lambda wildcards: wildcards.patient+ ".err",
        out=lambda wildcards: wildcards.patient+ ".out"
    threads: 1      
    resources:
        mem_mb=10000   
    shell:
        """
        bcftools +fill-tags {input.v} \
        --threads {threads} \
        -Oz -o {output.v} \
        -- -t FORMAT/VAF
        bcftools index --force --tbi --threads {threads} {output.v}
        """


rule VAF_filter:
    """
    Filter based on VAF. 
    """
    input:
        v="Mutect2/multi/tag_VAF/{patient}.vcf.gz",
        i="Mutect2/multi/tag_VAF/{patient}.vcf.gz.tbi"
    output:
        v="Mutect2/multi/tag_VAF/filtered_{patient}.vcf.gz"
    params:
        err=lambda wildcards: wildcards.patient+ ".err",
        out=lambda wildcards: wildcards.patient+ ".out"
    threads: 1      
    resources:
        mem_mb=10000   
    shell:
        """
        bcftools view --include 'VAF>0.1' {input} -Ov -o {output}
        """


################ SAREK ANNOTATION AT /nemo/project/proj-tracerX/working/SRM/RABIA/Mutect2/multi/tag_VAF/annotate.sh ############

rule vcf2maf:
    input:
        ref=REF_GEN,
        zvcf= lambda wildcards: '/nemo/project/proj-tracerX/working/SRM/RABIA/VEP/'+wildcards.rep+'/Annotation/'+wildcards.rep+'/VEP/filtered_' +recal.loc[recal['SampleName']==wildcards.tum,'PatientName'].iloc[0]+'_VEP.ann.vcf.gz'
    output:
        "MAF/{rep}/{tum}.maf"
    params:
        err=lambda wildcards: wildcards.tum+ ".err",
        out=lambda wildcards: wildcards.tum+ ".out",
        nn=get_normal_name,
        vcf= lambda wildcards, input: input.zvcf[:-3]
    threads: 1
    resources:
        mem_mb=5000
    envmodules:
        "HTSlib/1.15.1-GCC-11.3.0",
        "SAMtools/1.16.1-GCC-11.3.0"
    shell:
        """
        gunzip {input.zvcf}
        perl /camp/project/tracerX/working/PIPELINES/Analysis_Modules/vcf2maf/v1.0/vcf2maf.pl \
			--input-vcf {params.vcf} \
			--output-maf {output} \
			--tumor-id {wildcards.tum} \
			--normal-id {params.nn} \
			--inhibit-vep \
			--ref-fasta {input.ref}
        """


rule germline_SNP:
    """
    Call germline SNPs from all tumour and normal samples
    """
    input:
        lambda wildcards: recal.loc[(recal['SampleName']==wildcards.samp) ,'bai'],
        b= lambda wildcards: recal.loc[(recal['SampleName']==wildcards.samp) ,'bam'],
        ref=REF_GEN,
        intv=INT
    output:
        mpu=temp("GL_variant/{samp}.mpu"),
        mpui=temp("GL_variant/{samp}.mpu.tbi"),
        v="GL_variant/{samp}.vcf"
    params:
        err=lambda wildcards: wildcards.samp+ ".err",
        out=lambda wildcards: wildcards.samp+ ".out"
    threads: 1
    resources:
        mem_mb=5000
    envmodules:
        "HTSlib/1.15.1-GCC-11.3.0" 
    shell:
        """
        /nemo/lab/turajlics/home/users/fidanr/bcftools/bin/bcftools mpileup -f {input.ref} --annotate FORMAT/AD --max-depth 4000 {input.b} -Oz -o {output.mpu}
        /nemo/lab/turajlics/home/users/fidanr/bcftools/bin/bcftools index --tbi {output.mpu}
        bcftools call -mv -R {input.intv} {output.mpu} -Ov -o {output.v}
        """


rule concordanceQC:
    """
    Check if blood and tumour samples belong to the same patient.
    /camp/project/tracerX/working/PIPELINES/Analysis_Modules/qc_concordance/count_discordant.sh 
    """
    input:

    output:

    params:
        err=lambda wildcards: wildcards.patient+ ".err",
        out=lambda wildcards: wildcards.patient+ ".out"
    threads: 1      
    resources:
        mem_mb=10000   
    shell:
        """
        """