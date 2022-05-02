import shutil
from collections import OrderedDict

seqs = OrderedDict()

with open(snakemake.input[0]) as f:
    label = ""
    while True:
        l = f.readline().rstrip('\n')
        if not l:
            break
        elif l[0] == ">":
            label = l
            seqs[label] = ""
        else:
            seqs[label] += l

seqs = OrderedDict(sorted(seqs.items(), key=lambda t: -len(t[1]))) # descending order sort by length

# Remove seqs from end until they are greater than minlength
if seqs:
    while len(seqs[next(reversed(seqs))]) < int(snakemake.params.len):
        del seqs[next(reversed(seqs))]

with open(f"{snakemake.input[0]}.{snakemake.params.len}f", 'w') as f:
    for k, v in seqs.items():
        f.write(f"{k}\n")
        f.write(f"{v}\n")

shutil.copyfile(f"{snakemake.input[0]}.{snakemake.params.len}f", snakemake.output[0])