# -*- mode: Snakemake -*-
#
# Illumina quality control rules

from sunbeamlib import qc

rule all_qc:
    """Runs trimmomatic and fastqc on all input files."""
    input:
        TARGET_QC

ruleorder: adapter_removal_paired > adapter_removal_unpaired

rule sample_intake:
    input:
        lambda wildcards: Samples[wildcards.sample][wildcards.rp]
    output:
        str(QC_FP/'00_samples'/'{sample}_{rp}.fastq.gz')
    params:
        suffix = Cfg['qc']['seq_id_ending']
    conda:
        "../../envs/qc.yml"
    script:
        "../../scripts/qc/sample_intake.py"

rule adapter_removal_unpaired:
    input:
        str(QC_FP/'00_samples'/'{sample}_1.fastq.gz')
    params:
        tmp = str(QC_FP/'01_cutadapt'/'{sample}_1.fastq'),
    log: str(QC_FP/'log'/'cutadapt'/'{sample}.log') 
    output:
        str(QC_FP/'01_cutadapt'/'{sample}_1.fastq.gz')
    threads:
        Cfg['qc']['threads']
    conda:
        "../../envs/qc.yml"
    script:
        "../../scripts/qc/adapter_removal_unpaired.py"
    
rule adapter_removal_paired:
    input:
        r1 = str(QC_FP/'00_samples'/'{sample}_1.fastq.gz'),
        r2 = str(QC_FP/'00_samples'/'{sample}_2.fastq.gz')
    params:
        r1 = str(QC_FP/'01_cutadapt'/'{sample}_1.fastq.gz'),
        r2 = str(QC_FP/'01_cutadapt'/'{sample}_2.fastq.gz')
    log: str(QC_FP/'log'/'cutadapt'/'{sample}.log') 
    output:
        gr1 = str(QC_FP/'01_cutadapt'/'{sample}_1.fastq.gz'),
        gr2 = str(QC_FP/'01_cutadapt'/'{sample}_2.fastq.gz')
    threads:
        Cfg['qc']['threads']
    conda:
        "../../envs/qc.yml"
    script:
        "../../scripts/qc/adapter_removal_paired.py"

ruleorder: trimmomatic_paired > trimmomatic_unpaired
        
rule trimmomatic_unpaired:
    input:
        str(QC_FP/'01_cutadapt'/'{sample}_1.fastq.gz')
    output:
        str(QC_FP/'02_trimmomatic'/'{sample}_1.fastq.gz')
    log: str(QC_FP/'log'/'trimmomatic'/'{sample}.out')
    params:
        sw_start = Cfg['qc']['slidingwindow'][0],
        sw_end = Cfg['qc']['slidingwindow'][1]
    threads:
        Cfg['qc']['threads']
    conda:
        "../../envs/qc.yml"
    shell:
        """
        trimmomatic \
        SE -threads {threads} -phred33 \
        {input} {output} \
        ILLUMINACLIP:{Cfg[qc][adapter_template]}:2:30:10:8:true \
        LEADING:{Cfg[qc][leading]} \
        TRAILING:{Cfg[qc][trailing]} \
        SLIDINGWINDOW:{params.sw_start}:{params.sw_end} \
        MINLEN:{Cfg[qc][minlen]} \
        > >(tee -a {log}) 2> >(tee -a {log} >&2)
        """

            
rule trimmomatic_paired:
    input:
        r1 = str(QC_FP/'01_cutadapt'/'{sample}_1.fastq.gz'),
        r2 = str(QC_FP/'01_cutadapt'/'{sample}_2.fastq.gz')
    output:
        pair_r1 = str(QC_FP/'02_trimmomatic'/'{sample}_1.fastq.gz'),
        pair_r2 = str(QC_FP/'02_trimmomatic'/'{sample}_2.fastq.gz'),
        unpair_r1 = temp(str(QC_FP/'02_trimmomatic'/'unpaired'/'{sample}_1_unpaired.fastq.gz')),
        unpair_r2 = temp(str(QC_FP/'02_trimmomatic'/'unpaired'/'{sample}_2_unpaired.fastq.gz'))
    log: str(QC_FP/'log'/'trimmomatic'/'{sample}.out'),
    params:
        sw_start = Cfg['qc']['slidingwindow'][0],
        sw_end = Cfg['qc']['slidingwindow'][1]
    threads:
        Cfg['qc']['threads']
    conda:
        "../../envs/qc.yml"
    shell:
        """
        trimmomatic \
        PE -threads {threads} -phred33 \
        {input.r1} {input.r2} \
        {output.pair_r1} {output.unpair_r1} \
        {output.pair_r2} {output.unpair_r2} \
        ILLUMINACLIP:{Cfg[qc][adapter_template]}:2:30:10:8:true \
        LEADING:{Cfg[qc][leading]} \
        TRAILING:{Cfg[qc][trailing]} \
        SLIDINGWINDOW:{params.sw_start}:{params.sw_end} \
        MINLEN:{Cfg[qc][minlen]} \
        > >(tee -a {log}) 2> >(tee -a {log} >&2)
        """

rule fastqc:
    input:
        reads = expand(
            str(QC_FP/"02_trimmomatic"/"{{sample}}_{rp}.fastq.gz"),
            rp=Pairs)
    output:
        expand(
            str(QC_FP/'reports'/'{{sample}}_{rp}_fastqc/fastqc_data.txt'),
            rp=Pairs)
    params:
        outdir = str(QC_FP/'reports')
    conda:
        "../../envs/qc.yml"
    shell:
        "fastqc -o {params.outdir} {input.reads} -extract"

rule find_low_complexity:
    input:
        expand(
            str(QC_FP/'02_trimmomatic'/'{{sample}}_{rp}.fastq.gz'),
            rp=Pairs)
    output:
        str(QC_FP/'log'/'komplexity'/'{sample}.filtered_ids')
    conda:
        "../../envs/qc.yml"
    shell:
        """
        for rp in {input}; do
          gzip -dc $rp | kz | \
          awk '{{ if ($4<{Cfg[qc][kz_threshold]}) print $1 }}' >> {output}
        done
        """

rule remove_low_complexity:
    input:
        reads = str(QC_FP/'02_trimmomatic'/'{sample}_{rp}.fastq.gz'),
        ids = str(QC_FP/'log'/'komplexity'/'{sample}.filtered_ids')
    output:
        str(QC_FP/'03_komplexity'/'{sample}_{rp}.fastq.gz')
    conda:
        "../../envs/qc.yml"
    shell:
        """
        gzip -dc {input.reads} | rbt fastq-filter {input.ids} |\
        gzip > {output}
        """

rule qc_final:
    input:
        str(QC_FP/'03_komplexity'/'{sample}_{rp}.fastq.gz')
    output:
        str(QC_FP/'cleaned'/'{sample}_{rp}.fastq.gz')
    shell:
        """cp {input} {output}"""

rule clean_qc:
    input:
        expand(
            str(QC_FP/'cleaned'/'{sample}_{rp}.fastq.gz'),
            sample=Samples.keys(), rp=Pairs)
    params:
        cutadapt_fp = str(QC_FP/'01_cutadapt'),
        trimmomatic_fp = str(QC_FP/'02_trimmomatic'),
        komplexity_fp = str(QC_FP/'03_komplexity')
    output:
        touch(".qc_cleaned")
    shell:
        """
        rm -r {params.cutadapt_fp} && \
        rm -r {params.trimmomatic_fp} && \
        rm -r {params.komplexity_fp}
        """
