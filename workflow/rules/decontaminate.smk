rule all_decontam:
    input:
        TARGET_DECONTAM,


rule build_host_index:
    input:
        Cfg["qc"]["host_fp"] / "{host}.fasta",
    output:
        Cfg["qc"]["host_fp"] / "{host}.fasta.amb",
    log:
        LOG_FP / "build_host_index_{host}.log",
    benchmark:
        BENCHMARK_FP / "build_host_index_{host}.tsv"
    params:
        host="{host}",
        index_fp=Cfg["qc"]["host_fp"],
    conda:
        "../envs/qc.yml"
    shell:
        "cd {Cfg[qc][host_fp]} && bwa index {input} 2>&1 | tee {log}"


rule align_to_host:
    input:
        reads=expand(QC_FP / "cleaned" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
        index=Cfg["qc"]["host_fp"] / "{host}.fasta.amb",
        host=Cfg["qc"]["host_fp"] / "{host}.fasta",
    output:
        QC_FP / "decontam" / "intermediates" / "{host}" / "{sample}.sam",
    log:
        LOG_FP / "align_to_host_{host}_{sample}.log",
    benchmark:
        BENCHMARK_FP / "align_to_host_{host}_{sample}.tsv"
    threads: 4
    conda:
        "../envs/qc.yml"
    shell:
        """
        bwa mem -M -t {threads} {input.host} \
        {input.reads} -o {output} 2>&1 | tee {log}
        """


rule b_align_to_host:
    input:
        QC_FP / "decontam" / "intermediates" / "{host}" / "{sample}.sam",
    output:
        QC_FP / "decontam" / "intermediates" / "{host}" / "{sample}.bam",
    log:
        LOG_FP / "b_align_to_host_{host}_{sample}.log",
    benchmark:
        BENCHMARK_FP / "b_align_to_host_{host}_{sample}.tsv"
    threads: 4
    conda:
        "../envs/qc.yml"
    shell:
        """
        samtools view -bSF4 {input} -o {output} 2>&1 | tee {log}
        """


rule get_mapped_reads:
    input:
        QC_FP / "decontam" / "intermediates" / "{host}" / "{sample}.bam",
    output:
        ids=QC_FP / "decontam" / "intermediates" / "{host}" / "{sample}.ids",
    log:
        LOG_FP / "get_mapped_reads_{host}_{sample}.log",
    benchmark:
        BENCHMARK_FP / "get_mapped_reads_{host}_{sample}.tsv"
    params:
        pct_id=Cfg["qc"]["pct_id"],
        frac=Cfg["qc"]["frac"],
    conda:
        "../envs/reports.yml"
    script:
        "../scripts/get_mapped_reads.py"


rule aggregate_reads:
    input:
        expand(
            QC_FP / "decontam" / "intermediates" / "{host}" / "{{sample}}.ids",
            host=HostGenomes.keys(),
        ),
    output:
        temp(QC_FP / "decontam" / "intermediates" / "{sample}_hostreads.ids"),
    log:
        LOG_FP / "aggregate_reads_{sample}.log",
    shell:
        """
        arr=({input})
        if (( ${{#arr[@]}} == 0 )); then
            echo 'No intermediate ids...' > {log}
            touch {output}
        else
            echo 'Aggregating reads...' > {log}
            cat {input} > {output}
        fi
        """


rule filter_reads:
    input:
        hostreads=QC_FP / "decontam" / "intermediates" / "{sample}_hostreads.ids",
        reads=QC_FP / "cleaned" / "{sample}_{rp}.fastq.gz",
        hostids=expand(
            QC_FP / "decontam" / "intermediates" / "{host}" / "{{sample}}.ids",
            host=HostGenomes.keys(),
        ),
    output:
        reads=QC_FP / "decontam" / "{sample}_{rp}.fastq.gz",
        log=QC_FP / "log" / "decontam" / "{sample}_{rp}.txt",
    log:
        LOG_FP / "filter_reads_{sample}_{rp}.log",
    benchmark:
        BENCHMARK_FP / "filter_reads_{sample}_{rp}.tsv"
    conda:
        "../envs/rbt.yml"
    script:
        "../scripts/filter_reads.py"


rule preprocess_report:
    """Combines the information from multiple preprocessing steps"""
    input:
        trim_files=expand(
            LOG_FP / "trimmomatic_{sample}.log",
            sample=sorted(Samples.keys()),
        ),
        decontam_files=expand(
            QC_FP / "log" / "decontam" / "{sample}_1.txt",
            sample=sorted(Samples.keys()),
        ),
        komplexity_files=expand(
            QC_FP / "log" / "komplexity" / "{sample}.filtered_ids",
            sample=sorted(Samples.keys()),
        ),
    output:
        QC_FP / "reports" / "preprocess_summary.tsv",
    log:
        LOG_FP / "preprocess_report.log",
    benchmark:
        BENCHMARK_FP / "preprocess_report.tsv"
    conda:
        "../envs/reports.yml"
    script:
        "../scripts/preprocess_report.py"
