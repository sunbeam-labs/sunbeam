# Chunyu Zhao 20161202
# https://www.ncbi.nlm.nih.gov/genome/?term=Bacteroides+fragilis
# https://www.ncbi.nlm.nih.gov/genome/167


import sys
import os
import random
import gzip

# So that "random" values in this testing will actually be arbitrary but still
# deterministic:
random.seed(sum(bytes("sunbeam", "ascii")))


def reverse_complement(dna):
    dnadict = {"A": "T", "C": "G", "G": "C", "T": "A"}
    reverseDna = [dnadict[c] for c in dna]
    return "".join(reverseDna[::-1])


def write_fastq(genome_segment, locs, sample_name, file_fp, rp_suffix=None):
    """
    rp_suffix: a pair of strings to append to R1 and R2 sequence IDs, for
    example ["/1", "/2"].  By default nothing is appended.
    """
    seqs1 = []  # read1: read out of forward strand
    seqs2 = []  # read2: read out of reverse strand
    for loc in locs:
        # random start point; fragment distribution: gaussian
        start_pos, frag_len = map(int, loc.split("\t"))

        frag_reverse = reverse_complement(
            genome_segment[start_pos : start_pos + frag_len]
        )

        r1 = genome_segment[start_pos : start_pos + 250]
        r2 = frag_reverse[0:250]

        if len(r1) != 250 or len(r2) != 250:
            print(len(r1), len(r2), frag_len)

        seqs1.append(r1)
        seqs2.append(r2)

    ## write dummy fastq files
    file1 = gzip.open(file_fp + "/PCMP_" + sample_name + "_R1.fastq.gz", "wt")
    for i, seq in enumerate(seqs1):
        header = "@D00728:28:C9W1KANXX:" + str(i)
        if rp_suffix:
            header += rp_suffix[0]
        file1.write("%s\n" % header)
        file1.write("%s\n" % seq)
        file1.write("%s\n" % "+")
        file1.write("%s\n" % "".join("G" * 250))
    file1.close()

    file2 = gzip.open(file_fp + "/PCMP_" + sample_name + "_R2.fastq.gz", "wt")
    for i, seq in enumerate(seqs2):
        header = "@D00728:28:C9W1KANXX:" + str(i)
        if rp_suffix:
            header += rp_suffix[1]
        file2.write("%s\n" % header)
        file2.write("%s\n" % seq)
        file2.write("%s\n" % "+")
        file2.write("%s\n" % "".join("G" * 250))
    file2.close()


# Using this to trigger a failure to assemble any contigs in the assembly
# rules.  Note the random.seed() call above so that we will always get the same
# sequences.
def write_random_fastq(sample_name, file_fp, numreads=4, readlen=100):
    """Write a pair of sample fastq files containing randomized sequences."""
    dna = "ATCG"
    randread = lambda: "".join(
        [dna[random.randint(0, len(dna) - 1)] for i in range(readlen)]
    )
    for r in (1, 2):
        with gzip.open(
            file_fp + "/PCMP_" + sample_name + "_R%s.fastq.gz" % r, "wt"
        ) as f:
            for s in range(numreads):
                f.write("@read%s\n" % s)
                f.write("%s\n" % randread())
                f.write("+\n")
                f.write("%s\n" % ("G" * readlen))


def generate_dummyfastq(rootpath, fastq_fn="data_files", rp_suffix=None):

    raw_fp = rootpath + "/raw"

    # read in the 10k fragments
    with open(raw_fp + "/GCF_Ecoli_10k_genomic.fna") as f:
        lines_e = f.read().splitlines()

    with open(raw_fp + "/GCF_Bfragilis_10k_genomic.fna") as f:
        lines_b = f.read().splitlines()

    ## read in 1000 random start position and fragment length
    with open(raw_fp + "/index_10k.txt") as f:
        locs = f.read().splitlines()

    genome_segment_e = lines_e[1]
    genome_segment_b = lines_b[1]

    fastq_fp = rootpath + "/" + fastq_fn
    if not os.path.exists(fastq_fp):
        os.makedirs(fastq_fp)

    write_fastq(genome_segment_e, locs, "dummyecoli", fastq_fp, rp_suffix)
    write_fastq(genome_segment_b, locs, "dummybfragilis", fastq_fp, rp_suffix)
    write_random_fastq("random", fastq_fp)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        rootpath = sys.argv[1]
        rp_suffix = None
        fastq_fn = "data_files"
        if len(sys.argv) == 5:
            fastq_fn = sys.argv[2]
            rp_suffix = sys.argv[3:5]
    else:
        raise SystemExit("Must specify an output directory")
    generate_dummyfastq(rootpath, fastq_fn, rp_suffix)
