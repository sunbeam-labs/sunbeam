# -*- mode: Snakemake -*-
#
# Map reads to contigs and calculate per base coverage
#
# Requires Minimap2 and samtools.


rule all_coverage:
    input:
        ASSEMBLY_FP / "contigs_coverage.txt",


rule _contigs_mapping:
    input:
        expand(
            ASSEMBLY_FP / "contigs" / "coverage" / "{sample}.depth",
            sample=Samples.keys(),
        ),


rule _all_coverage:
    input:
        expand(
            ASSEMBLY_FP / "contigs" / "coverage" / "{sample}.csv",
            sample=Samples.keys(),
        ),


rule contigs_mapping:
    input:
        contig=ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa",
        reads=expand(QC_FP / "decontam" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
    output:
        bam=ASSEMBLY_FP / "contigs" / "minimap2" / "{sample}.sorted.bam",
        bai=ASSEMBLY_FP / "contigs" / "minimap2" / "{sample}.sorted.bam.bai",
    benchmark:
        BENCHMARK_FP / "contigs_mapping_{sample}.tsv"
    params:
        temp=temp(ASSEMBLY_FP / "contigs" / "minimap2" / "{sample}.sorted.tmp"),
    threads: 4  # Should be overridden by profile's set-threads (https://github.com/snakemake/snakemake/issues/1983)
    conda:
        "../../envs/assembly.yml"
    shell:
        """
        minimap2 -ax sr -t {threads} {input.contig} {input.reads} | \
            samtools sort -@ {threads} -o {output.bam} -T {params.temp} - 
        samtools index {output.bam} {output.bai}
        """


rule mapping_depth:
    input:
        bam=ASSEMBLY_FP / "contigs" / "minimap2" / "{sample}.sorted.bam",
        bai=ASSEMBLY_FP / "contigs" / "minimap2" / "{sample}.sorted.bam.bai",
    output:
        depth=ASSEMBLY_FP / "contigs" / "coverage" / "{sample}.depth",
    benchmark:
        BENCHMARK_FP / "mapping_depth_{sample}.tsv"
    conda:
        "../../envs/assembly.yml"
    shell:
        """
        samtools depth -aa {input.bam} > {output.depth}
        """


rule get_coverage:
    input:
        ASSEMBLY_FP / "contigs" / "coverage" / "{sample}.depth",
    output:
        ASSEMBLY_FP / "contigs" / "coverage" / "{sample}.csv",
    benchmark:
        BENCHMARK_FP / "get_coverage_{sample}.tsv"
    conda:
        "../../envs/assembly.yml"
    script:
        "../../scripts/assembly/get_coverage.py"


rule summarize_coverage:
    input:
        expand(
            ASSEMBLY_FP / "contigs" / "coverage" / "{sample}.csv",
            sample=Samples.keys(),
        ),
    output:
        ASSEMBLY_FP / "contigs_coverage.txt",
    shell:
        "(head -n 1 {input[0]}; tail -q -n +2 {input}) > {output}"
