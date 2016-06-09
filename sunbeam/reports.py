import warnings
warnings.filterwarnings('ignore', '.*experimental.*')
import csv
from pathlib import Path
from collections import Counter, defaultdict

from Bio import SeqIO

from Bio import SearchIO
from Bio.SeqRecord import SeqRecord

from .parsers import MetaGeneAnnotation

def blast_summary(blast_files, blast_format="blast-xml"):
    """Summarize BLAST results from an set of BLAST output files."""
    for infile in blast_files:
        sample = Path(infile).stem
        try:
            for result in SearchIO.parse(infile, blast_format):
                if len(result.hits) > 0:
                    yield {
                        'sample': sample,
                        'query': result.id,
                        'hit':result.hits[0].id
                    }
        except ParseError:
            print("Skipping empty/malformed %s" % infile)
            continue

def blast_contig_summary(xml_files):
    return {r['query']: r['hit'] for r in blast_summary(xml_files)}

def blast_hits(blast_xml_files):
    return Counter(c['query'] for c in blast_summary(blast_xml_files))
    
def mga_summary(orf_fp, contigs_fp, sample_id):
    """Summarizes results from MetaGene Annotator."""
    genes = []
    gene_no = 0
    contigs = SeqIO.parse(contigs_fp, 'fasta')
    contig_seqs = {r.description: r.seq for r in contigs}
    annotations = MetaGeneAnnotation.parse(open(orf_fp))
    for anno in annotations:
        contig_seq = contig_seqs[anno.id]
        # Gather putative genes and create SeqRecords for them
        for pgene in anno.genes:
            gene_seq = contig_seq[pgene.start:pgene.end]
            if pgene.strand == '-':
                gene_seq = gene_seq.reverse_complement()
            gene = SeqRecord(
                gene_seq,
                id=anno.id,
                description="{}:{}".format(sample_id, gene_no))
            genes.append(gene)
            gene_no += 1
    return genes


        
    
