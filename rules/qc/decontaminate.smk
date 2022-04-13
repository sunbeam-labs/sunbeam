rule all_decontam:
    input:
        expand(
            str(QC_FP/'decontam'/'{sample}_{rp}.fastq.gz'),
            sample=Samples.keys(), rp=Pairs)

ruleorder: build_host_index > build_genome_index
        
rule build_host_index:
    input:
        str(Cfg['qc']['host_fp']/'{host}.fasta')
    output:
        str(Cfg['qc']['host_fp']/'{host}.fasta.amb')
    params:
        host='{host}',
        index_fp=str(Cfg['qc']['host_fp'])
    conda:
        "../../envs/qc.yml"
    shell:
        "cd {Cfg[qc][host_fp]} && bwa index {input}"

rule align_to_host:
    input:
        reads = expand(
            str(QC_FP/'cleaned'/'{{sample}}_{rp}.fastq.gz'),
            rp = Pairs),
        index = str(Cfg['qc']['host_fp']/'{host}.fasta.amb')
    output:
        temp(str(QC_FP/'decontam'/'intermediates'/'{host}'/'{sample}.bam'))
    threads:
        Cfg['qc']['threads']
    params:
        sam = temp(str(QC_FP/'decontam'/'intermediates'/'{host}'/'{sample}.sam')),
        index_fp = str(Cfg['qc']['host_fp'])
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
        str(QC_FP/'decontam'/'intermediates'/'{host}'/'{sample}.bam')
    output:
        ids = str(QC_FP/'decontam'/'intermediates'/'{host}'/'{sample}.ids'),
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
            str(QC_FP/'decontam'/'intermediates'/'{host}'/'{{sample}}.ids'),
            host=HostGenomes.keys())
    output:
        temp(str(QC_FP/'decontam'/'intermediates'/'{sample}_hostreads.ids')),
    script:
        "../../scripts/qc/aggregate_reads.py"

rule filter_reads:
    input:
        hostreads = str(QC_FP/'decontam'/'intermediates'/'{sample}_hostreads.ids'),
        reads = str(QC_FP/'cleaned'/'{sample}_{rp}.fastq.gz'),
        hostids = expand(str(QC_FP/'decontam'/'intermediates'/'{host}'/'{{sample}}.ids'), host=HostGenomes.keys())
    output:
        reads = str(QC_FP/'decontam'/'{sample}_{rp}.fastq.gz'),
        log = str(QC_FP/'log'/'decontam'/'{sample}_{rp}.txt')
    conda:
        "../../envs/qc.yml"
    script:
        "../../scripts/qc/filter_reads.py"
