# Sunbeam test data generation

These reads are generated using `wgsim` to take 100 reads from each of four genomes: human, phiX 174, E. coli, and B. fragilis. These reference genomes are included under `tests/hosts` and `tests/raw`. Use the following commands to generate new reads:

```bash
cd tests/data/reads

wgsim ../raw/Bfragilis.fasta Bfragilis_R1.fq Bfragilis_R2.fq -N 100
wgsim ../raw/Ecoli.fasta Ecoli_R1.fq Ecoli_R2.fq -N 100
wgsim ../hosts/human.fasta human_R1.fq human_R2.fq -N 100
wgsim ../hosts/phix174.fasta phix174_R1.fq phix174_R2.fq -N 100

cat *_R1.fq > TEST_R1.fastq
cat *_R2.fq > TEST_R2.fastq

gzip *.fastq

rm *.fq
```

The `tests/data/single_end_reads` directory holds the `TEST_R1.fastq` file.