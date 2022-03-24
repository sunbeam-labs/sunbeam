# -*- mode: Snakemake -*-
#
# Map reads to contigs and calculate per base coverage
#
# Requires Minimap2 and samtools.

rule all_coverage:
    input:
        str(ASSEMBLY_FP/'contigs_coverage.txt')

rule _contigs_mapping:
    input:
        expand(str(ASSEMBLY_FP/'contigs'/'coverage'/'{sample}.depth'),
               sample=Samples.keys())        

rule _all_coverage:
    input:
        expand(str(ASSEMBLY_FP/'contigs'/'coverage'/'{sample}.csv'),
               sample=Samples.keys())

rule contigs_mapping:
    input:
        contig = str(ASSEMBLY_FP/'contigs'/'{sample}-contigs.fa'),
        reads = expand(str(QC_FP/'decontam'/'{{sample}}_{rp}.fastq.gz'),rp = Pairs)
    output:
        bam = str(ASSEMBLY_FP/'contigs'/'minimap2'/'{sample}.sorted.bam'),
        bai = str(ASSEMBLY_FP/'contigs'/'minimap2'/'{sample}.sorted.bam.bai')
    params:
        temp = temp(str(ASSEMBLY_FP/'contigs'/'minimap2'/'{sample}.sorted.tmp'))
    threads: 
        Cfg['mapping']['threads']
    conda:
        "../../envs/minimap2_samtools.yml"
    shell:
        """
        minimap2 -ax sr -t {threads} {input.contig} {input.reads} | \
            samtools sort -@ {threads} -o {output.bam} -T {params.temp} - 
        samtools index {output.bam} {output.bai}
        """

rule mapping_depth:
    input:
        bam = str(ASSEMBLY_FP/'contigs'/'minimap2'/'{sample}.sorted.bam'),
        bai = str(ASSEMBLY_FP/'contigs'/'minimap2'/'{sample}.sorted.bam.bai')
    output:
        depth = str(ASSEMBLY_FP/'contigs'/'coverage'/'{sample}.depth')
    conda:
        "../../envs/samtools.yml"
    shell:
        """
        samtools depth -aa {input.bam} > {output.depth}
        """


rule get_coverage:
    input:
        str(ASSEMBLY_FP/'contigs'/'coverage'/'{sample}.depth')
    output:
        str(ASSEMBLY_FP/'contigs'/'coverage'/'{sample}.csv')
    conda:
        "../../envs/numpy.yml"
    script:
        "../../scripts/assembly/get_coverage.py"


rule summarize_coverage:
    input:
        expand(str(ASSEMBLY_FP/'contigs'/'coverage'/'{sample}.csv'), 
               sample = Samples.keys())
    output:
        str(ASSEMBLY_FP/'contigs_coverage.txt')
    shell:
        "(head -n 1 {input[0]}; tail -q -n +2 {input}) > {output}"

