# -*- mode: Snakemake -*-
#
# Illumina quality control rules
localrules:
    all_qc,
    sample_intake,
    qc_final,
    clean_qc,


rule all_qc:
    """Runs trimmomatic and fastqc on all input files."""
    input:
        TARGET_QC,


rule sample_intake:
    input:
        lambda wildcards: Samples[wildcards.sample][wildcards.rp],
    output:
        QC_FP / "00_samples" / "{sample}_{rp}.fastq.gz",
    log:
        LOG_FP / "sample_intake_{sample}_{rp}.log",
    script:
        "../scripts/sample_intake.py"


rule adapter_removal:
    input:
        reads=expand(QC_FP / "00_samples" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
    output:
        reads=expand(QC_FP / "01_adapters" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
        fail=QC_FP / "01_adapters" / "{sample}_adapter_removal_failed.fastq.gz",
        json=QC_FP / "reports" / "01_{sample}_adapter_removal.json",
    log:
        LOG_FP / "adapter_removal_{sample}.log",
    benchmark:
        BENCHMARK_FP / "adapter_removal_{sample}.tsv"
    params:
        adapter=Cfg["qc"]["adapter_template"],
    resources:
        runtime=lambda wc, input: max(MIN_RUNTIME, input.size_mb / 5),
    threads: 4
    conda:
        "../envs/qc.yml"
    container:
        get_docker_str("qc")
    script:
        "../scripts/adapter_removal.py"


#rule 


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
    resources:
        mem_mb=lambda wc, input: max(MIN_MEM_MB, (input.size_mb / 2000) * MIN_MEM_MB),
        runtime=lambda wc: max(MIN_RUNTIME, 240),
    threads: 4
    conda:
        "../envs/qc.yml"
    container:
        get_docker_str("qc")
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
        > {log} 2>&1
        """

rule remove_low_complexity:
    input:
        reads=expand(
            QC_FP / "02_trimmomatic" / "{{sample}}_{rp}.fastq.gz",
            rp=Pairs,
        ),
    output:
        reads=expand(
            QC_FP / "cleaned" / "{{sample}}_{rp}.fastq.gz",
            rp=Pairs,
        ),
        counter=QC_FP / "reports" / "03_{sample}_count.txt",
    params:
        kmer_size=Cfg["qc"]["kmer_size"],
        min_kscore=Cfg["qc"]["kz_threshold"],
    log:
        LOG_FP / "remove_low_complexity_{sample}.log",
    benchmark:
        BENCHMARK_FP / "remove_low_complexity_{sample}.tsv"
    resources:
        mem_mb=lambda wc, input: max(MIN_MEM_MB, 2 * input.size_mb),
        runtime=lambda wc: max(MIN_RUNTIME, 120),
    conda:
        "../envs/qc.yml"
    container:
        get_docker_str("qc")
    script:
        "../scripts/remove_low_complexity.py"


rule fastqc:
    input:
        reads=expand(QC_FP / "cleaned" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
    output:
        expand(QC_FP / "reports" / "{{sample}}_{rp}_fastqc/fastqc_data.txt", rp=Pairs),
    log:
        LOG_FP / "fastqc_{sample}.log",
    benchmark:
        BENCHMARK_FP / "fastqc_{sample}.tsv"
    resources:
        runtime=lambda wc: max(MIN_RUNTIME, 120),
    conda:
        "../envs/qc.yml"
    container:
        get_docker_str("qc")
    shell:
        """
        sample_dir=$(dirname {output[0]})
        outdir=$(dirname $sample_dir)

        fastqc -o $outdir {input.reads} -extract > {log} 2>&1
        """


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
    script:
        "../scripts/fastqc_report.py"


rule clean_qc:
    input:
        expand(
            QC_FP / "cleaned" / "{sample}_{rp}.fastq.gz",
            sample=Samples.keys(),
            rp=Pairs,
        ),
    output:
        touch(QC_FP / ".qc_cleaned"),
    shell:
        """
        qc_dir=$(dirname {output[0]})

        rm -r $qc_dir/01_adapters || true
        rm -r $qc_dir/02_trimmomatic || true
        rm -r $qc_dir/03_complexity || true
        """
