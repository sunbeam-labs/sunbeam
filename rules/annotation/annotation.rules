# -*- mode: Snakemake -*-
#
# Contig annotation.
#
# See Readme.md

import csv
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from sunbeamlib import circular


rule all_annotate:
    input:
        TARGET_ANNOTATE

rule aggregate_results:
    input:
        contigs=str(ASSEMBLY_FP/'contigs'/'{sample}-contigs.fa'),
        contig_results=expand(
            str(ANNOTATION_FP/'blastn'/'{db}'/'contig'/'{{sample}}.xml'),
            db=Blastdbs['nucl']),
        gene_results=expand(
            str(ANNOTATION_FP/'{blastpx}'/'{db}'/'{orf_finder}'/'{{sample}}.xml'),
            blastpx=['blastp','blastx'],
            db=Blastdbs['prot'],
            orf_finder=['prodigal'])
    output:
        str(ANNOTATION_FP/'summary'/'{sample}.tsv')
    params:
        dbs=list(Blastdbs['nucl'].keys()) + list(Blastdbs['prot'].keys())
    run:
        contigs = {r.id: r.seq for r in SeqIO.parse(input.contigs, 'fasta')}
        # Separate each set of result files by the database it was blasted against
        contig_results = {
            db: blast_contig_summary(f for f in input.contig_results if db in f)
            for db in Blastdbs['nucl']
        }
        # We only care about the number of hits for the protein database
        gene_hits = {
            db: blast_hits(f for f in input.gene_results if db in f)
            for db in Blastdbs['prot']
        }
        with open(output[0], 'w') as out:
            writer = csv.DictWriter(
                out,
                fieldnames=['sample', 'contig', 'length', 'circular'] + params.dbs,
                delimiter='\t')
            writer.writeheader()
            for contig, contig_seq in contigs.items():
                is_circular = circular(
                    contig_seq,
                    Cfg['annotation']['circular_kmin'],
                    Cfg['annotation']['circular_kmax'],
                    Cfg['annotation']['circular_min_len'])
                results = {
                    'sample':wildcards.sample,
                    'contig':contig,
                    'length':len(contig_seq),
                    'circular':is_circular
                }
                for db in Blastdbs['nucl']:                
                    results[db] = contig_results[db].get(contig, "NA")
                # Report the number of hits of each contig/gene for each prot. db
                # Genes are reported from prodigal as contig_1,contig_2, etc so
                # aggregate all hits for all of a contig's genes together using sum
                for db in Blastdbs['prot']:
                    results[db] = sum(
                        gene_hits[db][gene] for gene in gene_hits[db] if contig in gene)
                writer.writerow(results)

rule aggregate_all:
    input:
        expand(
            str(ANNOTATION_FP/'summary'/'{sample}.tsv'),
            sample=Samples.keys())
    output:
        str(ANNOTATION_FP/'all_samples.tsv')
    run:
        with open(output[0], 'w') as out:
            out.writelines(open(input[0]).readlines())
            for infile in input[1:]:
                out.writelines(l for i,l in enumerate(open(infile)) if i > 0)
