localrules:
    all_decontam,
    aggregate_reads,
    clean_decontam,


rule all_decontam:
    input:
        TARGET_DECONTAM,


rule build_host_index:
    input:
        Cfg["qc"]["host_fp"] / "{host}.fasta",
    output:
        [
            Cfg["qc"]["host_fp"] / ("{host}.fasta." + ext)
            for ext in ["amb", "ann", "bwt", "pac", "sa"]
        ],
    log:
        LOG_FP / "build_host_index_{host}.log",
    benchmark:
        BENCHMARK_FP / "build_host_index_{host}.tsv"
    conda:
        "../envs/qc.yml"
    container:
        get_docker_str("qc")
    shell:
        """
        bwa index {input} 2>&1 | tee {log}
        """


rule align_to_host:
    input:
        reads=expand(QC_FP / "cleaned" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
        host=Cfg["qc"]["host_fp"] / "{host}.fasta",
        ids=[
            Cfg["qc"]["host_fp"] / ("{host}.fasta." + ext)
            for ext in ["amb", "ann", "bwt", "pac", "sa"]
        ],
    output:
        temp(QC_FP / "decontam" / "intermediates" / "{host}" / "{sample}.sam"),
    log:
        LOG_FP / "align_to_host_{host}_{sample}.log",
    benchmark:
        BENCHMARK_FP / "align_to_host_{host}_{sample}.tsv"
    resources:
        mem_mb=lambda wc: max(MIN_MEM_MB, 16000),
        runtime=lambda wc: max(MIN_RUNTIME, 240),
    threads: 4
    conda:
        "../envs/qc.yml"
    container:
        get_docker_str("qc")
    shell:
        """
        bwa mem -M -t {threads} {input.host} \
        {input.reads} -o {output} 2>&1 | tee {log}
        """


rule get_mapped_reads:
    input:
        QC_FP / "decontam" / "intermediates" / "{host}" / "{sample}.sam",
    output:
        ids=QC_FP / "decontam" / "intermediates" / "{host}" / "{sample}.ids",
    log:
        LOG_FP / "get_mapped_reads_{host}_{sample}.log",
    benchmark:
        BENCHMARK_FP / "get_mapped_reads_{host}_{sample}.tsv"
    params:
        pct_id=Cfg["qc"]["pct_id"],
        frac=Cfg["qc"]["frac"],
    script:
        "../scripts/get_mapped_reads.py"


rule aggregate_reads:
    input:
        expand(
            QC_FP / "decontam" / "intermediates" / "{host}" / "{{sample}}.ids",
            host=HostGenomes,
        ),
    output:
        temp(QC_FP / "decontam" / "intermediates" / "{sample}_hostreads.ids"),
    log:
        LOG_FP / "aggregate_reads_{sample}.log",
    params:
        host_fp=Cfg["qc"]["host_fp"],
    script:
        "../scripts/aggregate_reads.py"


rule filter_reads:
    input:
        hostreads=QC_FP / "decontam" / "intermediates" / "{sample}_hostreads.ids",
        reads=QC_FP / "cleaned" / "{sample}_{rp}.fastq.gz",
        hostids=expand(
            QC_FP / "decontam" / "intermediates" / "{host}" / "{{sample}}.ids",
            host=HostGenomes,
        ),
    output:
        reads=QC_FP / "decontam" / "{sample}_{rp}.fastq.gz",
        log=QC_FP / "log" / "decontam" / "{sample}_{rp}.txt",
    log:
        LOG_FP / "filter_reads_{sample}_{rp}.log",
    benchmark:
        BENCHMARK_FP / "filter_reads_{sample}_{rp}.tsv"
    resources:
        mem_mb=lambda wc: max(MIN_MEM_MB, 24000),
        runtime=lambda wc: max(MIN_RUNTIME, 240),
    script:
        "../scripts/filter_reads.py"


rule preprocess_report:
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
    script:
        "../scripts/preprocess_report.py"


rule clean_decontam:
    input:
        expand(
            QC_FP / "decontam" / "{sample}_{rp}.fastq.gz",
            sample=Samples.keys(),
            rp=Pairs,
        ),
        QC_FP / ".qc_cleaned",
    output:
        touch(QC_FP / ".decontam_cleaned"),
    shell:
        """
        cleaned_dir=$(dirname {input[0]})
        qc_dir=$(dirname $cleaned_dir)

        rm -r $qc_dir/cleaned || true
        rm -r $qc_dir/decontam/intermediates || true
        """
