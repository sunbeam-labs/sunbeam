# -*- mode: Snakemake -*-
#
# ReportGeneration rules

rule all_reports:
    input:
        TARGET_REPORT


rule preprocess_report:
    """Combines the information from multiple preprocessing steps"""
    input:
        trim_files = expand(
            QC_FP/'log'/'trimmomatic'/'{sample}.out',
            sample=sorted(Samples.keys())),
        decontam_files = expand(
            QC_FP/'log'/'decontam'/'{sample}_1.txt',
            sample=sorted(Samples.keys())),
        komplexity_files = expand(
            QC_FP/'log'/'komplexity'/'{sample}.filtered_ids',
            sample=sorted(Samples.keys())),
    output:
        QC_FP/'reports'/'preprocess_summary.tsv'
    conda:
        "../../envs/reports.yml"
    script:
        "../../scripts/reports/preprocess_report.py"


rule fastqc_report:
    """ make fastqc reports """
    input:
        files = expand(
            QC_FP/'reports'/'{sample}_{rp}_fastqc/fastqc_data.txt',
            sample=Samples.keys(),rp=Pairs)
    output:
        QC_FP/'reports'/'fastqc_quality.tsv'
    conda:
        "../../envs/reports.yml"
    script:
        "../../scripts/reports/fastqc_report.py"

