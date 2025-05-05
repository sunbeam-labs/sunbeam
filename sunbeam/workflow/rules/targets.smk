# -*- mode: Snakemake -*-
#
# list the all input files for each step

# ---- Quality control
# FastQC reports
TARGET_FASTQC = [
    expand(
        QC_FP / "reports" / "{sample}_{rp}_fastqc" / "fastqc_data.txt",
        sample=Samples.keys(),
        rp=Pairs,
    ),
    str(QC_FP / "reports" / "fastqc_quality.tsv"),
]


# Quality-control reads
TARGET_CLEAN = expand(
    QC_FP / "cleaned" / "{sample}_{rp}.fastq.gz", sample=Samples.keys(), rp=Pairs
)


TARGET_QC = TARGET_CLEAN + TARGET_FASTQC


# Remove host reads
TARGET_DECONTAM = [
    expand(
        QC_FP / "decontam" / "{sample}_{rp}.fastq.gz", sample=Samples.keys(), rp=Pairs
    ),
    str(QC_FP / "reports" / "preprocess_summary.tsv"),
]


# ---- All targets
TARGET_ALL = TARGET_QC + TARGET_DECONTAM
