import warnings
warnings.filterwarnings('ignore', '.*experimental.*')
from pathlib import Path
from collections import Counter

from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord

from xml.etree.ElementTree import ParseError

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
    


        
    
