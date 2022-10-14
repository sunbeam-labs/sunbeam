rule all_decontam:
    input:
        expand(
            QC_FP/'decontam'/'{sample}_{rp}.fastq.gz',
            sample=Samples.keys(), rp=Pairs)

ruleorder: build_host_index > build_genome_index
        
rule build_host_index:
    input:
        Cfg['qc']['host_fp']/'{host}.fasta'
    output:
        Cfg['qc']['host_fp']/'{host}.fasta.amb'
    params:
        host='{host}',
        index_fp=Cfg['qc']['host_fp']
    conda:
        "../../envs/qc.yml"
    shell:
        "cd {Cfg[qc][host_fp]} && bwa index {input}"

rule align_to_host:
    input:
        reads = expand(
            QC_FP/'cleaned'/'{{sample}}_{rp}.fastq.gz',
            rp = Pairs),
        index = Cfg['qc']['host_fp']/'{host}.fasta.amb'
    output:
        temp(QC_FP/'decontam'/'intermediates'/'{host}'/'{sample}.bam')
    threads:
        Cfg['qc']['threads']
    params:
        sam = temp(QC_FP/'decontam'/'intermediates'/'{host}'/'{sample}.sam'),
        index_fp = Cfg['qc']['host_fp']
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
        QC_FP/'decontam'/'intermediates'/'{host}'/'{sample}.bam'
    output:
        ids = QC_FP/'decontam'/'intermediates'/'{host}'/'{sample}.ids',
    params:
        pct_id =  Cfg['qc']['pct_id'],
        frac = Cfg['qc']['frac']
    conda:
        "../../envs/qc.yml"
    script:
        "../../scripts/qc/get_mapped_reads.py"

rule aggregate_reads:
    input:
        expand(
            QC_FP/'decontam'/'intermediates'/'{host}'/'{{sample}}.ids',
            host=HostGenomes.keys())
    output:
        temp(QC_FP/'decontam'/'intermediates'/'{sample}_hostreads.ids'),
    script:
        "../../scripts/qc/aggregate_reads.py"

rule filter_reads:
    input:
        hostreads = QC_FP/'decontam'/'intermediates'/'{sample}_hostreads.ids',
        reads = QC_FP/'cleaned'/'{sample}_{rp}.fastq.gz',
        hostids = expand(QC_FP/'decontam'/'intermediates'/'{host}'/'{{sample}}.ids', host=HostGenomes.keys())
    output:
        reads = QC_FP/'decontam'/'{sample}_{rp}.fastq.gz',
        log = QC_FP/'log'/'decontam'/'{sample}_{rp}.txt'
    conda:
        "../../envs/qc.yml"
    script:
        "../../scripts/qc/filter_reads.py"

rule preprocess_report:
    """Combines the information from multiple preprocessing steps"""
    input:
        trim_files = expand(
            QC_FP/'log'/'trimmomatic'/'{sample}.out',
            sample=sorted(Samples.keys())),
        decontam_files = expand(
            QC_FP/'log'/'decontam'/'{sample}_1.txt',
            sample=sorted(Samples.keys())),
        komplexity_files = expand(
            QC_FP/'log'/'komplexity'/'{sample}.filtered_ids',
            sample=sorted(Samples.keys())),
    output:
        QC_FP/'reports'/'preprocess_summary.tsv'
    conda:
        "../../envs/reports.yml"
    script:
        "../../scripts/reports/preprocess_report.py"