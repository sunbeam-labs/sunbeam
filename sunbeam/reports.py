import csv

from Bio import SeqIO
from Bio import SearchIO

from .parsers import MetaGeneAnnotation

def blast_summary(blast_xml_files, out_handle):
    """Create a summary of BLAST results from an set of XML files."""
    writer = csv.DictWriter(
	out_handle,
	fieldnames=['sample','contig','hit'], delimiter='\t')
    writer.writeheader()
    for infile in blast_xml_files:
        sample = Path(infile).stem
        try:
            results = [
                {
                    'sample': sample,
                    'query': result.id,
                    'hit':result.hits[0].id
                }
                for result in SearchIO.parse(infile, 'blast-xml')
                if len(result.hits) > 0
            ]
            writer.writerows(results)
        except ParseError:
            print("Skipping empty/malformed %s" % infile)
            continue

def mga_summary(orf_fp, contigs_fp, sample_id):
    """Summarizes results from MetaGene Annotator."""
    genes = []
    gene_no = 0
    contigs = SeqIO.parse(contigs_fp, 'fasta')
    contig_seqs = {r.description: r.seq for r in contigs}
    annotations = MetaGeneAnnotation.parse(open(orf_fp))
    for anno in annotations:
        contig_seq = contigs[anno.id]
        # Gather putative genes and create SeqRecords for them
        for pgene in anno.genes:
            gene_seq = contig_seq[pgene.start:pgene.end]
            if pgene.strand == '-':
                gene_seq = gene_seq.reverse_complement()
            gene = SeqRecord(
                gene_seq,
                id="{}:{}".format(sample_id, gene_no),
                description=anno.id)
            genes.append(gene)
            gene_no += 1
    return genes


        
    
