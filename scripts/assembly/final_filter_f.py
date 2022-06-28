from io import TextIOWrapper

def parse_fasta(f: TextIOWrapper) -> list:
    desc = ""
    seq = ""
    for l in f:
        l = l.strip()
        if l[0] == ">":
            yield desc, seq
            desc = l
            seq = ""
        else:
            seq += l
    yield desc, seq

def filter_seqs(seqs: list, min: int) -> list:
    for desc, seq in seqs:
        if len(seq) >= min:
            yield desc, seq

def write_fasta(f: TextIOWrapper, seqs: list):
    seqs.sort(key=lambda t: -len(t[1])) # Forces everything into mem at once but doesn't seem like there's a good way around that
    for desc, seq in seqs:
        f.write(f"{desc}\n{seq}\n")