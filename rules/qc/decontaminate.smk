rule all_decontam:
    input:
        TARGET_DECONTAM,


ruleorder: build_host_index > build_genome_index


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
        "../../envs/qc.yml"
    shell:
        "cd {Cfg[qc][host_fp]} && bwa index {input}"


rule align_to_host:
    input:
        reads=expand(QC_FP / "cleaned" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
        index=Cfg["qc"]["host_fp"] / "{host}.fasta.amb",
    output:
        temp(QC_FP / "decontam" / "intermediates" / "{host}" / "{sample}.bam"),
    log:
        LOG_FP / "align_to_host_{host}_{sample}.log",
    benchmark:
        BENCHMARK_FP / "align_to_host_{host}_{sample}.tsv"
    params:
        sam=temp(QC_FP / "decontam" / "intermediates" / "{host}" / "{sample}.sam"),
        index_fp=Cfg["qc"]["host_fp"],
    threads: 4  # Should be overridden by profile's set-threads (https://github.com/snakemake/snakemake/issues/1983)
    conda:
        "../../envs/qc.yml"
    shell:
        """
        bwa mem -M -t {threads} \
        {params.index_fp}/{wildcards.host}.fasta \
        {input.reads} -o {params.sam} && \
        samtools view -bSF4 {params.sam} > {output} && \
        rm {params.sam}
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
        "../../envs/qc.yml"
    script:
        "../../scripts/qc/get_mapped_reads.py"


rule aggregate_reads:
    input:
        expand(
            QC_FP / "decontam" / "intermediates" / "{host}" / "{{sample}}.ids",
            host=HostGenomes.keys(),
        ),
    output:
        temp(QC_FP / "decontam" / "intermediates" / "{sample}_hostreads.ids"),
    script:
        "../../scripts/qc/aggregate_reads.py"


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
        QC_FP / "log" / "decontam" / "{sample}_{rp}.txt",
    benchmark:
        BENCHMARK_FP / "filter_reads_{sample}_{rp}.tsv"
    conda:
        "../../envs/qc.yml"
    script:
        "../../scripts/qc/filter_reads.py"


rule preprocess_report:
    """Combines the information from multiple preprocessing steps"""
    input:
        trim_files=expand(
            QC_FP / "log" / "trimmomatic" / "{sample}.out",
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
        "../../envs/reports.yml"
    script:
        "../../scripts/qc/preprocess_report.py"
