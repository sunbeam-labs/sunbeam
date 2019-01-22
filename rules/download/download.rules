# -*- mode: Snakemake -*-
#
# Rules to download SRA data

ruleorder: download_paired > download_unpaired

rule download_paired:
    output:
        r1 = str(DOWNLOAD_FP/'{sample}_1.fastq.gz'),
        r2 = str(DOWNLOAD_FP/'{sample}_2.fastq.gz')
    params:
        accession = '{sample}',
        outdir = str(DOWNLOAD_FP)
    threads:
        Cfg['download']['threads']
    shadow: "shallow"
    shell:
        """
        mkdir -p {params.outdir}
        grabseqs sra -r 3 -t {threads} -o {params.outdir} {params.accession}
        """

rule download_unpaired:
    output:
        str(DOWNLOAD_FP/'{sample}.fastq.gz')
    params:
        accession = '{sample}',
        outdir = str(DOWNLOAD_FP)
    threads:
        Cfg['download']['threads']	
    shadow: "shallow"
    shell:
        """
        mkdir -p {params.outdir}
        grabseqs sra -r 3 -t {threads} -o {params.outdir} {params.accession}
        """
