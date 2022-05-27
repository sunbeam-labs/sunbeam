import csv
from Bio import SeqIO
from sunbeamlib import circular
from sunbeamlib.reports import blast_contig_summary, blast_hits

contigs = {r.id: r.seq for r in SeqIO.parse(snakemake.input.contigs, 'fasta')}
# Separate each set of result files by the database it was blasted against
contig_results = {
    db: blast_contig_summary(f for f in snakemake.input.contig_results if db in f)
    for db in snakemake.params.nucl
}
# We only care about the number of hits for the protein database
gene_hits = {
    db: blast_hits(f for f in snakemake.input.gene_results if db in f)
    for db in snakemake.params.prot
}

with open(snakemake.output[0], 'w') as out:
    writer = csv.DictWriter(
        out,
        fieldnames=['sample', 'contig', 'length', 'circular'] + snakemake.params.dbs,
        delimiter='\t')
    writer.writeheader()
    for contig, contig_seq in contigs.items():
        is_circular = circular(
            contig_seq,
            snakemake.config['annotation']['circular_kmin'],
            snakemake.config['annotation']['circular_kmax'],
            snakemake.config['annotation']['circular_min_len'])
        results = {
            'sample':snakemake.wildcards.sample,
            'contig':contig,
            'length':len(contig_seq),
            'circular':is_circular
        }
        for db in snakemake.params.nucl:                
            results[db] = contig_results[db].get(contig, "NA")
        # Report the number of hits of each contig/gene for each prot. db
        # Genes are reported from prodigal as contig_1,contig_2, etc so
        # aggregate all hits for all of a contig's genes together using sum
        for db in snakemake.params.prot:
            results[db] = sum(
                gene_hits[db][gene] for gene in gene_hits[db] if contig in gene)
        writer.writerow(results)