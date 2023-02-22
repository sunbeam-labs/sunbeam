# Sunbeam test data generation

These reads are generated using `wgsim` to take 100 reads from each of four genomes: human, phiX 174, E. coli, and B. fragilis. These reference genomes are included under `.tests/hosts` and `.tests/raw`. Use the following commands to generate new reads:

```bash
cd .tests/reads

wgsim ../raw/Bfragilis.fna Bfragilis_R1.fq Bfragilis_R2.fq -N 100
wgsim ../raw/Ecoli.fna Ecoli_R1.fq Ecoli_R2.fq -N 100
wgsim ../hosts/human.fasta human_R1.fq human_R2.fq -N 100
wgsim ../hosts/phix174.fasta phix174_R1.fq phix174_R2.fq -N 100

cat *_R1.fq > TEST_R1.fastq
cat *_R2.fq > TEST_R2.fastq

rm *.fq
```