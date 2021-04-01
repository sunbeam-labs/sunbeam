# -*- mode: Snakemake -*-
#
# Rules for running Kraken

rule all_classify:
    input:
        TARGET_CLASSIFY
        
rule kraken2_classify_report:
    input:
        expand(str(QC_FP/'decontam'/'{{sample}}_{rp}.fastq.gz'),rp = Pairs)
    output:
        raw = str(CLASSIFY_FP/'kraken'/'raw'/'{sample}-raw.tsv'),
        report = str(CLASSIFY_FP/'kraken'/'{sample}-taxa.tsv')
    params:
        db = Cfg['classify']['kraken_db_fp'],
        paired_end = "--paired" if Cfg['all']['paired_end'] else ""
    threads:
        Cfg['classify']['threads']
    shell:
        """
        kraken2 --gzip-compressed \
                --db {params.db} \
                --report {output.report} \
                {params.paired_end} {input} \
                > {output.raw}
        """

rule kraken2_biom:
    input:
        expand(str(CLASSIFY_FP/'kraken'/'{sample}-taxa.tsv'),
               sample=Samples.keys())
    output:
        str(CLASSIFY_FP/'kraken'/'all_samples.biom')
    shell:
        """
        kraken-biom --max D -o {output} {input}
        """

rule classic_k2_biom:
    input:
        str(CLASSIFY_FP/'kraken'/'all_samples.biom')
    output:
        str(CLASSIFY_FP/'kraken'/'all_samples.tsv')
    shell:
        """
        biom convert -i {input} -o {output} \
        --to-tsv --header-key=taxonomy --process-obs-metadata=taxonomy \
        --output-metadata-id="Consensus Lineage"
        """

