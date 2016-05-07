from collections import namedtuple

class ContigAnnotation(object):

    def __init__(self, contig_id, model):
        self.id = contig_id
        self.model = model
        self.genes = []

    def __repr__(self):
        return self.id

    @staticmethod
    def parse(lines):
        contig_id = ''
        model = ''
        contig = None
        in_record = True
        for line in lines:
            if line.startswith("# gc ="):
                # I don't care about this stuff
                continue

            elif line.startswith("# self:"):
                model = line[8]
                continue

            elif line.startswith("# "):

                if contig is not None:
                    # We've reached the next record, so yield the current contig
                    # and start a new one
                    yield contig
                    contig = None
                    
                contig_id = line.strip("# \n")
                continue

            else:
                # No comment region, now parsing gene calls
                if contig is None:
                    contig = ContigAnnotation(contig_id, model)
                gene = PutativeGene(*line.strip("\n").split("\t"))
                contig.genes.append(gene)
    

class PutativeGene(object):

    completeness = {
        '11':'both',
        '01':'stop',
        '10':'start',
        '00':'neither'}

    def __init__(self, geneid, start, end, strand, frame, complete, score, model,
                 rbs_start, rbs_end, rbs_score):
        self.geneid = geneid
        self.start = int(start)
        self.end = int(end)
        self.length = self.end - self.start
        assert strand in ['-','+']
        self.strand = strand
        self.frame = int(frame)
        self.sscodons = self.completeness[complete]
        self.score = float(score)
        assert model in ['s','b','a','p']
        self.model = model
        self.rbs_start = int(rbs_start) if rbs_start != '-' else None
        self.rbs_end = int(rbs_end) if rbs_end != '-' else None
        self.rbs_score = float(rbs_score) if rbs_score != '-' else None

    def __repr__(self):
        return "{}bp '{}'-sense gene".format(self.length, self.strand)

