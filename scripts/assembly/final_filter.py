from io import TextIOWrapper
import shutil
#from collections import OrderedDict

#seqs = OrderedDict()

#with open(snakemake.input[0]) as f:
#    label = ""
#    while True:
#        l = f.readline().rstrip('\n')
#        if not l:
#            break
#        elif l[0] == ">":
#            label = l
#            seqs[label] = ""
#        else:
#            seqs[label] += l

#seqs = OrderedDict(sorted(seqs.items(), key=lambda t: -len(t[1]))) # descending order sort by length

# Remove seqs from end until they are greater than minlength
#if seqs:
#    while len(seqs[next(reversed(seqs))]) < int(snakemake.params.len):
#        del seqs[next(reversed(seqs))]

#with open(f"{snakemake.input[0]}.{snakemake.params.len}f", 'w') as f:
#    for k, v in seqs.items():
#        f.write(f"{k}\n")
#        f.write(f"{v}\n")

#shutil.copyfile(f"{snakemake.input[0]}.{snakemake.params.len}f", snakemake.output[0])

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

with open(snakemake.input[0]) as f:
    with open(f"{snakemake.input[0]}.{snakemake.params.len}f", "w") as g:
        write_fasta(g, list(filter_seqs(parse_fasta(f), snakemake.params.len)))

shutil.copyfile(f"{snakemake.input[0]}.{snakemake.params.len}f", snakemake.output[0])