# -*- mode: Snakemake -*-
#
# Contig annotation:
# 	Rules for BLASTing against databases
#
# See Readme.md

import csv

from Bio import SearchIO
from xml.etree.ElementTree import ParseError

#rule make_blast_db_prot:
#    """Use makeblastdb to create any necessary db indeces that don't exist."""
#    input:
#        prot = lambda wildcards: Blastdbs['prot'][wildcard.db]
#    output:
#        expand(str(GENES_DIR/'{{gene}}.fasta.{index}'),index=['psq','pin','phr'])
#    conda:
#        "../../envs/annotation.yml"
#    shell:
#        """
#        makeblastdb -in {input} -dbtype prot
#        """
    
#rule make_blast_db_nucl:
#    """Use makeblastdb to create any necessary db indeces that don't exist."""
#    input:
#        prot = lambda wildcards: Blastdbs['nucl'][wildcard.db]
#    output:
#        expand(str(GENES_DIR/'{{gene}}.fasta.{index}'),index=['psq','pin','phr'])
#    conda:
#        "../../envs/annotation.yml"
#    shell:
#        """
#        makeblastdb -in {input} -dbtype nucl
#        """

rule build_diamond_db:
    """Use diamond makedb to create any necessary db indeces that don't exist."""
    input:
        [Blastdbs['prot'][db] for db in Blastdbs['prot']]
    output:
        [Blastdbs['prot'][db] + '.dmnd' for db in Blastdbs['prot']]
    conda:
        "../../envs/annotation.yml"
    shell:
        """
        diamond makedb --in {input} -d {input}
        """

rule run_blastn:
    """Run BLASTn against a given database and write the results to blast tabular format."""    
    input:
        contigs=str(ASSEMBLY_FP/'contigs'/'{sample}-contigs.fa')
    output:
        str(ANNOTATION_FP/'blastn'/'{db}'/'{contigs}'/'{sample}.btf')
    params:
        db=lambda wildcard: Blastdbs['nucl'][wildcard.db] 
    threads: 
        Cfg['blast']['threads']
    conda:
        "../../envs/annotation.yml"
    shell:
        """
        blastn \
        -query {input.contigs} \
        -db {params.db} \
        -outfmt 7 \
        -num_threads {threads} \
        -evalue 1e-10 \
        -max_target_seqs 5000 \
        -out {output} \
        """

rule run_blastp:
    """Run BLASTp on translated genes against a target db and write to blast tabular format."""
    input:
        genes=str(ANNOTATION_FP/'genes'/'{orf_finder}'/'{sample}_genes_prot.fa'),
        db=lambda wildcard: Blastdbs['prot'][wildcard.db]
    output:
        str(ANNOTATION_FP/'blastp'/'{db}'/'{orf_finder}'/'{sample}.btf')
    threads:
        Cfg['blast']['threads']
    conda:
        "../../envs/annotation.yml"
    shell:
        """
        blastp \
        -query {input.genes} \
        -db {input.db} \
        -outfmt 7 \
        -num_threads {threads} \
        -evalue 1e-10 \
        -max_target_seqs 2475 \
        -out {output} \
        """

rule run_diamond_blastp:
    """Run diamond blastp on translated genes against a target db and write to blast tabular format."""
    input:
        genes=str(ANNOTATION_FP/'genes'/'{orf_finder}'/'{sample}_genes_prot.fa'),
        db=lambda wildcard: Blastdbs['prot'][wildcard.db],
        indeces=rules.build_diamond_db.output
    output:
        str(ANNOTATION_FP/'blastp'/'{db}'/'{orf_finder}'/'{sample}.btf_IGNORE_RULE')
    threads:
        Cfg['blast']['threads']
    conda:
        "../../envs/annotation.yml"
    shell:
        """
        diamond blastp \
        -q {input.genes} \
        --db {input.db} \
        --outfmt 6 \
        --threads {threads} \
        --evalue 1e-10 \
        --max-target-seqs 2475 \
        --out {output} \
        || if [ $? == 1 ]; then echo "Caught empty query error from diamond" && touch {output}; fi
        """

rule run_blastx:
    """Run BLASTx on untranslated genes against a target db and write to blast tabular format."""
    input:
        genes=str(ANNOTATION_FP/'genes'/'{orf_finder}'/'{sample}_genes_nucl.fa'),
        db=lambda wildcard: Blastdbs['prot'][wildcard.db]
    output:
        str(ANNOTATION_FP/'blastx'/'{db}'/'{orf_finder}'/'{sample}.btf_IGNORE_RULE')
    threads:
        Cfg['blast']['threads']
    conda:
        "../../envs/annotation.yml"
    shell:
        """
        blastx \
        -query {input.genes} \
        -db {input.db} \
        -outfmt 7 \
        -num_threads {threads} \
        -evalue 1e-10 \
        -max_target_seqs 2475 \
        -out {output} \
        && cat {output}
        """

rule run_diamond_blastx:
    """Run diamond blastx on untranslated genes against a target db and write to blast tabular format."""
    input:
        genes=str(ANNOTATION_FP/'genes'/'{orf_finder}'/'{sample}_genes_nucl.fa'),
        db=lambda wildcard: Blastdbs['prot'][wildcard.db],
        indeces=rules.build_diamond_db.output
    output:
        str(ANNOTATION_FP/'blastx'/'{db}'/'{orf_finder}'/'{sample}.btf')
    threads:
        Cfg['blast']['threads']
    conda:
        "../../envs/annotation.yml"
    shell:
        """
        diamond blastx \
        -q {input.genes} \
        --db {input.db} \
        --outfmt 6 \
        --threads {threads} \
        --evalue 1e-10 \
        --max-target-seqs 2475 \
        --out {output} \
        || if [ $? == 1 ]; then echo "Caught empty query error from diamond" && touch {output}; fi
        """
        
rule blast_report:
    """Create a summary of results from a BLAST call."""
    input:
        expand(
            str(ANNOTATION_FP/'{{blast_prog}}'/'{{db}}'/'{{query}}'/'{sample}.btf'),
            sample=Samples.keys())
    output:
        str(ANNOTATION_FP/'{blast_prog}'/'{db}'/'{query}'/'report.tsv')
    conda:
        "../../envs/annotation.yml"
    script:
        "../../scripts/annotation/blast_report.py"

rule _test_blastpx:
    input:
        expand(str(ANNOTATION_FP/'{blastpx}'/'card'/'prodigal'/'{sample}.btf'), 
               blastpx=['blastx','blastp'], sample=Samples.keys())
    
rule _test_blastpx_report:
    input:
        expand(str(ANNOTATION_FP/'{blastpx}'/'card'/'prodigal'/'report.tsv'),
        blastpx=['blastx','blastp'])

#rule clean_xml:
#    input:
#        expand(str(ANNOTATION_FP/'summary'/'{sample}.tsv'), sample=Samples.keys())
#    params:
#        blastn_fp = str(ANNOTATION_FP/'blastn'),
#        blastp_fp = str(ANNOTATION_FP/'blastp'),
#        blastx_fp = str(ANNOTATION_FP/'blastx')
#    output:
#        touch(".xml_cleaned")
#    shell:
#        """
#        if [ -d {params.blastn_fp} ]; then rm -r {params.blastn_fp}; fi && \
#        if [ -d {params.blastp_fp} ]; then rm -r {params.blastp_fp}; fi && \
#        if [ -d {params.blastx_fp} ]; then rm -r {params.blastx_fp}; fi
#        """
