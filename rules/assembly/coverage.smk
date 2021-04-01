# -*- mode: Snakemake -*-
#
# Map reads to contigs and calculate per base coverage
#
# Requires Minimap2 and samtools.


import csv
import numpy


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
    shell:
        """
        samtools depth -aa {input.bam} > {output.depth}
        """

def _get_coverage(input_fp, sample, output_fp):
    """
    Summarize stats for coverage data for each sample 
    """

    with open(input_fp) as f:
        reader = csv.reader(f, delimiter='\t')    
        data = {}
        for row in reader:
            if not data.get(row[0]):
                data[row[0]] = []
            data[row[0]].append(int(row[2]))
    
    # summarize stats for all segments present and append to output
    output_rows = []
    for segment in data.keys():
        sumval     = numpy.sum(data[segment])
        minval     = numpy.min(data[segment])
        maxval     = numpy.max(data[segment])
        mean       = numpy.mean(data[segment])
        median     = numpy.median(data[segment])
        stddev     = numpy.std(data[segment])
        gen_cov    = len(list(filter(lambda x: x!=0, data[segment])))
        gen_length = len(data[segment])
        row = [sample, segment, sumval, minval, maxval, mean, median, stddev, gen_cov, gen_length]
        output_rows.append(row)

    # write out stats per segment per sample
    fields = ['sample','contig', 'sum', 'min', 'max', 'mean', 'median', 'stddev', 'coverage', 'length']
    with open(output_fp, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        writer.writerows(output_rows)


rule get_coverage:
    input:
        str(ASSEMBLY_FP/'contigs'/'coverage'/'{sample}.depth')
    output:
        str(ASSEMBLY_FP/'contigs'/'coverage'/'{sample}.csv')
    run:
        _get_coverage(input[0], wildcards.sample, output[0])


rule summarize_coverage:
    input:
        expand(str(ASSEMBLY_FP/'contigs'/'coverage'/'{sample}.csv'), 
               sample = Samples.keys())
    output:
        str(ASSEMBLY_FP/'contigs_coverage.txt')
    shell:
        "(head -n 1 {input[0]}; tail -q -n +2 {input}) > {output}"

