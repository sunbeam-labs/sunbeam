# Chunyu Zhao 20161202
# https://www.ncbi.nlm.nih.gov/genome/?term=Bacteroides+fragilis
# https://www.ncbi.nlm.nih.gov/genome/167


import sys
import os

def reverse_complement(dna):
	dnadict = {'A':'T','C':'G','G':'C','T':'A'}
	reverseDna = [ dnadict[c] for c in dna ]
	return ''.join(reverseDna[::-1])

def write_fastq(genome_segment, locs, sample_name,file_fp):
	seqs1 = [] # read1: read out of forward strand
	seqs2 = [] # read2: read out of reverse strand 
	for loc in locs:
		# random start point; fragment distribution: gaussian
		start_pos, frag_len = map(int,loc.split("\t"))

		frag_reverse = reverse_complement(genome_segment[start_pos:start_pos+frag_len])

		r1 = genome_segment[start_pos:start_pos+250]
		r2 = frag_reverse[0:250]

		if len(r1) != 250 or len(r2) != 250:
			print(len(r1), len(r2), frag_len)

		seqs1.append(r1)
		seqs2.append(r2)

	## write dummy fastq files
	file1 = open(file_fp + '/PCMP_'+sample_name+'_R1.fastq', 'w')
	for i, seq in enumerate(seqs1):
		header = "@D00728:28:C9W1KANXX:" + str(i)
		file1.write("%s\n" % header)
		file1.write("%s\n" % seq)
		file1.write("%s\n" % "+")
		file1.write("%s\n" % ''.join("G"*250))
	file1.close()

	file2 = open(file_fp + '/PCMP_'+sample_name+'_R2.fastq', 'w')
	for i, seq in enumerate(seqs2):
		header = "@D00728:28:C9W1KANXX:" + str(i)
		file2.write("%s\n" % header)
		file2.write("%s\n" % seq)
		file2.write("%s\n" % "+")
		file2.write("%s\n" % ''.join("G"*250))
	file2.close()

def generate_dummyfastq(rootpath):

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

	fastq_fp = rootpath + "/data_files"
	if not os.path.exists(fastq_fp):
		os.makedirs(fastq_fp)

	write_fastq(genome_segment_e, locs, "dummyecoli", fastq_fp)
	write_fastq(genome_segment_b, locs, "dummybfragilis",fastq_fp)


if __name__ == '__main__':
	if len(sys.argv) == 2:
		rootpath = sys.argv[1]
	else:
		print("Use the same root path in test_config.yml; Update me when needed.")
		rootpath = "/home/chunyu/sunbeam/tests"

	generate_dummyfastq(rootpath)
