# -*- mode: Snakemake -*-
#
# Illumina quality control rules

from sunbeamlib import qc


rule all_qc:
    """Runs trimmomatic and fastqc on all input files."""
    input:
        TARGET_QC,


ruleorder: adapter_removal_paired > adapter_removal_unpaired


rule sample_intake:
    input:
        lambda wildcards: Samples[wildcards.sample][wildcards.rp],
    output:
        QC_FP / "00_samples" / "{sample}_{rp}.fastq.gz",
    log:
        LOG_FP / "sample_intake_{sample}_{rp}.log",
    script:
        "../scripts/sample_intake.py"


rule adapter_removal_unpaired:
    input:
        QC_FP / "00_samples" / "{sample}_1.fastq.gz",
    output:
        QC_FP / "01_cutadapt" / "{sample}_1.fastq.gz",
    log:
        LOG_FP / "adapter_removal_unpaired_{sample}.log",
    benchmark:
        BENCHMARK_FP / "adapter_removal_unpaired_{sample}.tsv"
    params:
        str(QC_FP / "01_cutadapt" / "{sample}_1.fastq"),
    threads: 4
    conda:
        "../envs/cutadapt.yml"
    script:
        "../scripts/adapter_removal_unpaired.py"


rule adapter_removal_paired:
    input:
        r1=QC_FP / "00_samples" / "{sample}_1.fastq.gz",
        r2=QC_FP / "00_samples" / "{sample}_2.fastq.gz",
    output:
        r1=QC_FP / "01_cutadapt" / "{sample}_1.fastq.gz",
        r2=QC_FP / "01_cutadapt" / "{sample}_2.fastq.gz",
    log:
        LOG_FP / "adapter_removal_paired_{sample}.log",
    benchmark:
        BENCHMARK_FP / "adapter_removal_paired_{sample}.tsv"
    params:
        r1=str(QC_FP / "01_cutadapt" / "{sample}_1.fastq"),
        r2=str(QC_FP / "01_cutadapt" / "{sample}_2.fastq"),
    threads: 4
    conda:
        "../envs/cutadapt.yml"
    script:
        "../scripts/adapter_removal_paired.py"


ruleorder: trimmomatic_paired > trimmomatic_unpaired


rule trimmomatic_unpaired:
    input:
        QC_FP / "01_cutadapt" / "{sample}_1.fastq.gz",
    output:
        QC_FP / "02_trimmomatic" / "{sample}_1.fastq.gz",
    log:
        LOG_FP / "trimmomatic_{sample}.log",
    benchmark:
        BENCHMARK_FP / "trimmomatic_unpaired_{sample}.tsv"
    params:
        sw_start=Cfg["qc"]["slidingwindow"][0],
        sw_end=Cfg["qc"]["slidingwindow"][1],
    threads: 4
    conda:
        "../envs/qc.yml"
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
        > >(tee {log}) 2> >(tee {log} >&2)
        """


rule trimmomatic_paired:
    input:
        r1=QC_FP / "01_cutadapt" / "{sample}_1.fastq.gz",
        r2=QC_FP / "01_cutadapt" / "{sample}_2.fastq.gz",
    output:
        pair_r1=QC_FP / "02_trimmomatic" / "{sample}_1.fastq.gz",
        pair_r2=QC_FP / "02_trimmomatic" / "{sample}_2.fastq.gz",
        unpair_r1=temp(
            QC_FP / "02_trimmomatic" / "unpaired" / "{sample}_1_unpaired.fastq.gz"
        ),
        unpair_r2=temp(
            QC_FP / "02_trimmomatic" / "unpaired" / "{sample}_2_unpaired.fastq.gz"
        ),
    log:
        LOG_FP / "trimmomatic_{sample}.log",
    benchmark:
        BENCHMARK_FP / "trimmomatic_paired_{sample}.tsv"
    params:
        sw_start=Cfg["qc"]["slidingwindow"][0],
        sw_end=Cfg["qc"]["slidingwindow"][1],
    threads: 4
    conda:
        "../envs/qc.yml"
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
        > >(tee {log}) 2> >(tee {log} >&2)
        """


rule fastqc:
    input:
        reads=expand(QC_FP / "02_trimmomatic" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
    output:
        expand(QC_FP / "reports" / "{{sample}}_{rp}_fastqc/fastqc_data.txt", rp=Pairs),
    log:
        LOG_FP / "fastqc_{sample}.log",
    benchmark:
        BENCHMARK_FP / "fastqc_{sample}.tsv"
    params:
        outdir=QC_FP / "reports",
    conda:
        "../envs/qc.yml"
    shell:
        "fastqc -o {params.outdir} {input.reads} -extract 2>&1 | tee {log}"


rule fastqc_report:
    """ make fastqc reports """
    input:
        files=expand(
            QC_FP / "reports" / "{sample}_{rp}_fastqc/fastqc_data.txt",
            sample=Samples.keys(),
            rp=Pairs,
        ),
    output:
        QC_FP / "reports" / "fastqc_quality.tsv",
    log:
        LOG_FP / "fastqc_report.log",
    benchmark:
        BENCHMARK_FP / "fastqc_report.tsv"
    conda:
        "../envs/reports.yml"
    script:
        "../scripts/fastqc_report.py"


# REMOVE THIS ONCE IT'S ON PYPI
rule install_heyfastq:
    output:
        touch(QC_FP / ".heyfastq_installed"),
    conda:
        "../envs/reports.yml"
    shell:
        """
        git clone https://github.com/kylebittinger/heyfastq.git
        pip install heyfastq/
        """


rule remove_low_complexity:
    input:
        reads=expand(QC_FP / "02_trimmomatic" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
        installed=QC_FP / ".heyfastq_installed",
    output:
        reads=expand(QC_FP / "03_komplexity" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
    log:
        LOG_FP / "remove_low_complexity_{sample}.log",
    benchmark:
        BENCHMARK_FP / "remove_low_complexity_{sample}.tsv"
    conda:
        "../envs/reports.yml"
    params:
        min_kscore=Cfg["qc"]["kz_threshold"],
        ir1=QC_FP / "03_komplexity" / "{{sample}}_1.fastq",
        ir2=QC_FP / "03_komplexity" / "{sample}_2.fastq",
        or1=QC_FP / "03_komplexity" / "{sample}_1.fastq",
        or2=QC_FP / "03_komplexity" / "{sample}_2.fastq",
    shell:
        """
        gzip -d {input.reads}
        heyfastq filter-kscore --input {params.ir1} {params.ir2} --output {params.or1} {params.or2} --min-kscore {params.min_kscore}
        gzip {input.reads} {output.reads}
        """


rule qc_final:
    input:
        QC_FP / "03_komplexity" / "{sample}_{rp}.fastq.gz",
    output:
        QC_FP / "cleaned" / "{sample}_{rp}.fastq.gz",
    shell:
        """cp {input} {output}"""


rule clean_qc:
    input:
        expand(
            QC_FP / "cleaned" / "{sample}_{rp}.fastq.gz",
            sample=Samples.keys(),
            rp=Pairs,
        ),
    params:
        cutadapt_fp=QC_FP / "01_cutadapt",
        trimmomatic_fp=QC_FP / "02_trimmomatic",
        komplexity_fp=QC_FP / "03_komplexity",
    output:
        touch(QC_FP / ".qc_cleaned"),
    shell:
        """
        rm -r {params.cutadapt_fp} || true
        rm -r {params.trimmomatic_fp} || true
        rm -r {params.komplexity_fp} || true
        """
