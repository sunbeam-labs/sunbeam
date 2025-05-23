# -*- mode: Snakemake -*-
#
# Illumina quality control rules

import os
from sunbeam import get_docker_str


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


rule fastp:
    """Run fastp on the input files."""
    input:
        r1=QC_FP / "00_samples" / "{sample}_1.fastq.gz",
        r2=QC_FP / "00_samples" / "{sample}_2.fastq.gz",
    output:
        r1=QC_FP / "01_fastp" / "{sample}_1.fastq.gz",
        r2=QC_FP / "01_fastp" / "{sample}_2.fastq.gz",
        ur1=QC_FP / "01_fastp" / "{sample}_unpaired_1.fastq.gz",
        ur2=QC_FP / "01_fastp" / "{sample}_unpaired_2.fastq.gz",
        json=QC_FP / "01_fastp" / "{sample}.json",
        html=QC_FP / "01_fastp" / "{sample}.html",
    params:
        adapter_template=Cfg["qc"]["adapter_template"],
        leading=Cfg["qc"]["leading"],
        trailing=Cfg["qc"]["trailing"],
        sw_start=Cfg["qc"]["slidingwindow"][0],
        sw_end=Cfg["qc"]["slidingwindow"][1],
        minlen=Cfg["qc"]["minlen"],
    log:
        LOG_FP / "fastp_{sample}_{rp}.log",
    benchmark:
        BENCHMARK_FP / "fastp_{sample}_{rp}.tsv"
    resources:
        mem_mb=lambda wc, input: max(MIN_MEM_MB, 2 * input.size_mb),
        runtime=lambda wc: max(MIN_RUNTIME, 120),
    threads: 4
    conda:
        "../envs/fastp.yml"
    container:
        get_docker_str("fastp")
    shell:
        """
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            --unpaired1 {output.ur1} \
            --unpaired2 {output.ur2} \
            --adapter_fasta {params.adapter_template} \
            --thread {threads} \
            --cut_front --cut_front_window_size {params.leading} \
            --cut_tail --cut_tail_window_size {params.trailing} \
            --cut_right --cut_right_window_size {params.sw_start} \
            --cut_mean_quality {params.sw_end} \
            --length_required {params.minlen} \
            --json {output.json} --html {output.html} 2>&1 | tee {log}
        """


rule fastqc:
    input:
        reads=expand(QC_FP / "01_fastp" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
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

        fastqc -o $outdir {input.reads} -extract 2>&1 | tee {log}
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


rule find_low_complexity:
    input:
        expand(QC_FP / "01_fastp" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
    output:
        QC_FP / "log" / "komplexity" / "{sample}.filtered_ids",
    log:
        LOG_FP / "find_low_complexity_{sample}.log",
    benchmark:
        BENCHMARK_FP / "find_low_complexity_{sample}.tsv"
    conda:
        "../envs/komplexity.yml"
    container:
        get_docker_str("komplexity")
    shell:
        """
        for rp in {input}; do
          gzip -dc $rp | kz | \
          awk '{{ if ($4<{Cfg[qc][kz_threshold]}) print $1 }}' >> {output}
        done
        """


rule remove_low_complexity:
    input:
        reads=QC_FP / "01_fastp" / "{sample}_{rp}.fastq.gz",
        ids=QC_FP / "log" / "komplexity" / "{sample}.filtered_ids",
    output:
        QC_FP / "03_komplexity" / "{sample}_{rp}.fastq.gz",
    log:
        LOG_FP / "remove_low_complexity_{sample}_{rp}.log",
    benchmark:
        BENCHMARK_FP / "remove_low_complexity_{sample}_{rp}.tsv"
    resources:
        mem_mb=lambda wc, input: max(MIN_MEM_MB, 2 * input.size_mb),
        runtime=lambda wc: max(MIN_RUNTIME, 120),
    script:
        "../scripts/remove_low_complexity.py"


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
    output:
        touch(QC_FP / ".qc_cleaned"),
    shell:
        """
        cleaned_dir=$(dirname {input[0]})
        qc_dir=$(dirname $cleaned_dir)

        rm -r $qc_dir/01_cutadapt || true
        rm -r $qc_dir/02_trimmomatic || true
        rm -r $qc_dir/03_komplexity || true
        """
