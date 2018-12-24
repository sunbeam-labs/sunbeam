# -*- mode: Snakemake -*-
#
# Contig annotation:
# 	Rules for BLASTing against databases
#
# See Readme.md

import csv

from Bio import SearchIO
from xml.etree.ElementTree import ParseError

rule run_blastn:
    """Run BLASTn against a given database and write the results to XML."""    
    input:
        contigs=str(ASSEMBLY_FP/'contigs'/'{sample}-contigs.fa')
    output:
        str(ANNOTATION_FP/'blastn'/'{db}'/'{contigs}'/'{sample}.xml')
    params:
        db=lambda wildcard: Blastdbs['nucl'][wildcard.db] 
    threads: 
        Cfg['blast']['threads']
    shell:
        """
        blastn \
        -query {input.contigs} \
        -db {params.db} \
        -outfmt 5 \
        -num_threads {threads} \
        -evalue 1e-10 \
        -max_target_seqs 5000 \
        -out {output} \
        """

rule run_blastp:
    """Run BLASTp on translated genes against a target db and write to XML."""
    input:
        genes=str(ANNOTATION_FP/'genes'/'{orf_finder}'/'{sample}_genes_prot.fa'),
        db=lambda wildcard: Blastdbs['prot'][wildcard.db]
    output:
        str(ANNOTATION_FP/'blastp'/'{db}'/'{orf_finder}'/'{sample}.xml')
    threads:
        Cfg['blast']['threads']
    shell:
        """
        blastp \
        -query {input.genes} \
        -db {input.db} \
        -outfmt 5 \
        -num_threads {threads} \
        -evalue 1e-10 \
        -max_target_seqs 2475 \
        -out {output} \
        """

rule run_blastx:
    """Run BLASTx on untranslated genes against a target db and write to XML."""
    input:
        genes=str(ANNOTATION_FP/'genes'/'{orf_finder}'/'{sample}_genes_nucl.fa'),
        db=lambda wildcard: Blastdbs['prot'][wildcard.db]
    output:
        str(ANNOTATION_FP/'blastx'/'{db}'/'{orf_finder}'/'{sample}.xml')
    threads:
        Cfg['blast']['threads']
    shell:
        """
        blastx \
        -query {input.genes} \
        -db {input.db} \
        -outfmt 5 \
        -num_threads {threads} \
        -evalue 1e-10 \
        -max_target_seqs 2475 \
        -out {output} \
        """
        
rule blast_report:
    """Create a summary of results from a BLAST call."""
    input:
        expand(
            str(ANNOTATION_FP/'{{blast_prog}}'/'{{db}}'/'{{query}}'/'{sample}.xml'),
            sample=Samples.keys())
    output:
        str(ANNOTATION_FP/'{blast_prog}'/'{db}'/'{query}'/'report.tsv')
    run:
        with open(output[0], 'w') as out:
            writer = csv.DictWriter(
	        out,
	        fieldnames=['sample','query','hit'],
                delimiter='\t')
            writer.writeheader()
            list(writer.writerow(result) for result in blast_summary(input))

rule _test_blastpx:
    input:
        expand(str(ANNOTATION_FP/'{blastpx}'/'card'/'prodigal'/'{sample}.xml'), 
               blastpx=['blastx','blastp'], sample=Samples.keys())
    
rule _test_blastpx_report:
    input:
        expand(str(ANNOTATION_FP/'{blastpx}'/'card'/'prodigal'/'report.tsv'),
        blastpx=['blastx','blastp'])

rule clean_xml:
    input:
        expand(str(ANNOTATION_FP/'summary'/'{sample}.tsv'), sample=Samples.keys())
    params:
        blastn_fp = str(ANNOTATION_FP/'blastn'),
        blastp_fp = str(ANNOTATION_FP/'blastp'),
        blastx_fp = str(ANNOTATION_FP/'blastx')
    output:
        touch(".xml_cleaned")
    shell:
        """
        if [ -d {params.blastn_fp} ]; then rm -r {params.blastn_fp}; fi && \
        if [ -d {params.blastp_fp} ]; then rm -r {params.blastp_fp}; fi && \
        if [ -d {params.blastx_fp} ]; then rm -r {params.blastx_fp}; fi
        """
