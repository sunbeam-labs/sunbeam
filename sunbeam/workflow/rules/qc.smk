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


ruleorder: adapter_removal_paired > adapter_removal_unpaired


rule adapter_removal_unpaired:
    input:
        QC_FP / "00_samples" / "{sample}_1.fastq.gz",
    output:
        r=QC_FP / "01_cutadapt" / "{sample}_1.fastq.gz",
        ngz=temp(QC_FP / "01_cutadapt" / "{sample}_1.fastq"),
    log:
        LOG_FP / "adapter_removal_unpaired_{sample}.log",
    benchmark:
        BENCHMARK_FP / "adapter_removal_unpaired_{sample}.tsv"
    resources:
        runtime=lambda wc, input: max(MIN_RUNTIME, input.size_mb / 5),
    threads: 4
    conda:
        "../envs/cutadapt.yml"
    container:
        get_docker_str("cutadapt")
    script:
        "../scripts/adapter_removal_unpaired.py"


rule adapter_removal_paired:
    input:
        r1=QC_FP / "00_samples" / "{sample}_1.fastq.gz",
        r2=QC_FP / "00_samples" / "{sample}_2.fastq.gz",
    output:
        r1=QC_FP / "01_cutadapt" / "{sample}_1.fastq.gz",
        r2=QC_FP / "01_cutadapt" / "{sample}_2.fastq.gz",
        ngz1=temp(QC_FP / "01_cutadapt" / "{sample}_1.fastq"),
        ngz2=temp(QC_FP / "01_cutadapt" / "{sample}_2.fastq"),
    log:
        LOG_FP / "adapter_removal_paired_{sample}.log",
    benchmark:
        BENCHMARK_FP / "adapter_removal_paired_{sample}.tsv"
    resources:
        runtime=lambda wc, input: max(MIN_RUNTIME, input.size_mb / 10),
    threads: 4
    conda:
        "../envs/cutadapt.yml"
    container:
        get_docker_str("cutadapt")
    script:
        "../scripts/adapter_removal_paired.py"


ruleorder: trimmomatic_paired > trimmomatic_unpaired


# Check that the adapter template file exists
if os.environ.get("SUNBEAM_NO_ADAPTER", None):
    assert os.path.exists(Cfg["qc"]["adapter_template"])
    assert os.path.isfile(Cfg["qc"]["adapter_template"])
    assert os.stat(Cfg["qc"]["adapter_template"]).st_size > 0


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
    resources:
        mem_mb=lambda wc, input: max(MIN_MEM_MB, (input.size_mb / 1000) * MIN_MEM_MB),
        runtime=lambda wc: max(MIN_RUNTIME, 240),
    threads: 4
    conda:
        "../envs/qc.yml"
    container:
        get_docker_str("qc")
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
        2>&1 | tee {log}
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
        2>&1 | tee {log}
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
        expand(QC_FP / "02_trimmomatic" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
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
        reads=QC_FP / "02_trimmomatic" / "{sample}_{rp}.fastq.gz",
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
