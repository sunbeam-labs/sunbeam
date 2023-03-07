rule all_decontam:
    input:
        TARGET_DECONTAM,


if Cfg["all"]["paired_end"]:

    ruleorder: sam_convert_paired > sam_convert_unpaired

else:

    ruleorder: sam_convert_unpaired > sam_convert_paired


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
    threads: 4
    conda:
        "../envs/qc.yml"
    shell:
        """
        bwa mem -M -t {threads} {input.host} \
        {input.reads} -o {output} 2>&1 | tee {log}
        """


rule get_unmapped_reads:
    input:
        QC_FP / "decontam" / "intermediates" / "{host}" / "{sample}.sam",
    output:
        temp(QC_FP / "decontam" / "{host}" / "unmapped_{sample}.sam"),
    log:
        LOG_FP / "get_unmapped_reads_{host}_{sample}.log",
    benchmark:
        BENCHMARK_FP / "get_unmapped_reads_{host}_{sample}.tsv"
    threads: 4
    conda:
        "../envs/qc.yml"
    shell:
        """
        samtools view -f4 {input} -o {output} 2>&1 | tee {log}
        """


rule sam_convert_unpaired:
    input:
        QC_FP / "decontam" / "{host}" / "unmapped_{sample}.sam",
    output:
        temp(QC_FP / "decontam" / "{host}" / "unmapped_{sample}_1.fastq"),
    log:
        LOG_FP / "sam_convert_{host}_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sam_convert_{host}_{sample}.tsv"
    script:
        "../scripts/sam_convert_unpaired.py"


rule sam_convert_paired:
    input:
        QC_FP / "decontam" / "{host}" / "unmapped_{sample}.sam",
    output:
        rp_1=temp(QC_FP / "decontam" / "{host}" / "unmapped_{sample}_1.fastq"),
        rp_2=temp(QC_FP / "decontam" / "{host}" / "unmapped_{sample}_2.fastq"),
    log:
        LOG_FP / "sam_convert_{host}_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sam_convert_{host}_{sample}.tsv"
    script:
        "../scripts/sam_convert_paired.py"


rule filter_unmapped_reads:
    input:
        unmapped_reads=expand(
            QC_FP / "decontam" / "{host}" / "unmapped_{{sample}}_{{rp}}.fastq",
            host=HostGenomes.keys(),
        ),
        reads=QC_FP / "cleaned" / "{sample}_{rp}.fastq.gz",
    output:
        reads=QC_FP / "decontam" / "{sample}_{rp}.fastq.gz",
        log=QC_FP / "log" / "decontam" / "{sample}_{rp}.txt",
    log:
        LOG_FP / "filter_unmapped_reads_{sample}_{rp}.log",
    benchmark:
        BENCHMARK_FP / "filter_unmapped_reads_{sample}_{rp}.tsv"
    script:
        "../scripts/filter_unmapped_reads.py"


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
    conda:
        "../envs/reports.yml"
    script:
        "../scripts/preprocess_report.py"
