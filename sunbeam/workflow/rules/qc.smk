# -*- mode: Snakemake -*-


localrules:
    all_qc,
    sample_intake,
    clean_qc,


rule all_qc:
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


def fastp_in(wildcards, input):
    if len(input.reads) > 1:
        return f"--in1 {input.reads[0]} --in2 {input.reads[1]}"
    else:
        return f"--in1 {input.reads[0]}"


def fastp_out(wildcards, output):
    if len(output.reads) > 1:
        return f"--out1 {output.reads[0]} --out2 {output.reads[1]}"
    else:
        return f"--out1 {output.reads[0]}"


rule adapter_removal:
    input:
        reads=expand(QC_FP / "00_samples" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
    output:
        reads=expand(QC_FP / "01_adapters" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
        json=QC_FP / "reports" / "01_{sample}_adapter_removal.json",
        html=QC_FP / "reports" / "01_{sample}_adapter_removal.html",
    log:
        LOG_FP / "adapter_removal_{sample}.log",
    benchmark:
        BENCHMARK_FP / "adapter_removal_{sample}.tsv"
    conda:
        "../envs/qc.yml"
    container:
        get_docker_str("qc")
    threads: Cfg["qc"].get("threads", os.cpu_count())
    resources:
        runtime=lambda wc, input: max(MIN_RUNTIME, input.size_mb / 5),
    params:
        inargs=fastp_in,
        outargs=fastp_out,
        adapter=Cfg["qc"]["adapter_fp"],
    shell:
        """
        fastp {params.inargs} {params.outargs} \
          --adapter_fasta {params.adapter} \
          --disable_quality_filtering \
          --disable_length_filtering \
          --json {output.json} \
          --html {output.html} \
          --thread {threads}
        """


rule trim_quality:
    input:
        reads=expand(QC_FP / "01_adapters" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
    output:
        reads=expand(QC_FP / "02_quality" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
        report=QC_FP / "reports" / "02_quality_{sample}.json",
    log:
        LOG_FP / "trim_quality_{sample}.log",
    benchmark:
        BENCHMARK_FP / "trim_quality_{sample}.tsv"
    conda:
        "../envs/qc.yml"
    threads: Cfg["qc"].get("threads", os.cpu_count())
    resources:
        mem_mb=lambda wc, input: max(MIN_MEM_MB, (input.size_mb / 2000) * MIN_MEM_MB),
        runtime=lambda wc, input: max(MIN_RUNTIME, input.size_mb / 5),
    params:
        window=Cfg["qc"]["slidingwindow"],
        start_threshold=Cfg["qc"]["leading"],
        end_threshold=Cfg["qc"]["trailing"],
        min_length=Cfg["qc"]["minlen"],
    shell:
        """
        heyfastq trim-qual \
          --input {input.reads} \
          --output {output.reads} \
          --report {output.report} \
          --window-width {params.window[0]} \
          --window-threshold {params.window[1]} \
          --start-threshold {params.start_threshold} \
          --end-threshold {params.end_threshold} \
          --min-length {params.min_length} \
          --threads {threads}
        """


rule remove_low_complexity:
    input:
        reads=expand(
            QC_FP / "02_quality" / "{{sample}}_{rp}.fastq.gz",
            rp=Pairs,
        ),
    output:
        reads=expand(
            QC_FP / "cleaned" / "{{sample}}_{rp}.fastq.gz",
            rp=Pairs,
        ),
        report=QC_FP / "reports" / "03_complexity_{sample}.json",
    log:
        LOG_FP / "remove_low_complexity_{sample}.log",
    benchmark:
        BENCHMARK_FP / "remove_low_complexity_{sample}.tsv"
    conda:
        "../envs/qc.yml"
    container:
        get_docker_str("qc")
    threads: Cfg["qc"].get("threads", os.cpu_count())
    resources:
        mem_mb=lambda wc, input: max(MIN_MEM_MB, 2 * input.size_mb),
        runtime=lambda wc: max(MIN_RUNTIME, 120),
    params:
        kmer_size=Cfg["qc"]["kmer_size"],
        min_kscore=Cfg["qc"]["kz_threshold"],
    shell:
        """
        heyfastq filter-kscore \
          --input {input.reads} \
          --output {output.reads} \
          --report {output.report} \
          --kmer-size {params.kmer_size} \
          --min-kscore {params.min_kscore} \
          --threads {threads}
        """


rule fastqc:
    input:
        reads=expand(QC_FP / "cleaned" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
    output:
        expand(QC_FP / "reports" / "{{sample}}_{rp}_fastqc/fastqc_data.txt", rp=Pairs),
    log:
        LOG_FP / "fastqc_{sample}.log",
    benchmark:
        BENCHMARK_FP / "fastqc_{sample}.tsv"
    conda:
        "../envs/qc.yml"
    container:
        get_docker_str("qc")
    resources:
        runtime=lambda wc: max(MIN_RUNTIME, 120),
    shell:
        """
        sample_dir=$(dirname {output[0]})
        outdir=$(dirname $sample_dir)

        fastqc -o $outdir {input.reads} -extract > {log} 2>&1
        """


rule fastqc_report:
    input:
        reports=expand(
            QC_FP / "reports" / "{sample}_{rp}_fastqc/fastqc_data.txt",
            sample=Samples.keys(),
            rp=Pairs,
        ),
    output:
        report=QC_FP / "reports" / "fastqc_quality.tsv",
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
