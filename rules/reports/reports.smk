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
            str(QC_FP/'log'/'trimmomatic'/'{sample}.out'),
            sample=sorted(Samples.keys())),
        decontam_files = expand(
            str(QC_FP/'log'/'decontam'/'{sample}_1.txt'),
            sample=sorted(Samples.keys()))
    output:
        str(QC_FP/'reports'/'preprocess_summary.tsv')
    conda:
        "../../envs/pandas.yml"
    script:
        "../../scripts/reports/preprocess_report.py"


rule fastqc_report:
    """ make fastqc reports """
    input:
        files = expand(
            str(QC_FP/'reports'/'{sample}_{rp}_fastqc/fastqc_data.txt'),
            sample=Samples.keys(),rp=Pairs)
    output:
        str(QC_FP/'reports'/'fastqc_quality.tsv')
    conda:
        "../../envs/pandas.yml"
    script:
        "../../scripts/reports/fastqc_report.py"


rule multiqc_report:
    """Build multiqc report on qc step."""
    input:
        TARGET_FASTQC
    output:
        MULTIQC_REPORT
    params:
        title = Cfg['qc'].get('report_title', 'QC report'),
        outdir = str(QC_FP/'reports')
    conda:
        "../../envs/multiqc.yml"
    script:
        "../../scripts/reports/multiqc_report.py"
