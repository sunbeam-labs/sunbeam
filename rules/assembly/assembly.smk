# -*- mode: Snakemake -*-
#
# Contig building and other assembly rules
#
# Requires Megahit.

rule all_assembly:
    """Build contigs for all samples."""
    input:
        TARGET_ASSEMBLY

ruleorder: megahit_paired > megahit_unpaired

rule megahit_paired:
    input:
        r1 = str(QC_FP/'decontam'/'{sample}_1.fastq.gz'),
        r2 = str(QC_FP/'decontam'/'{sample}_2.fastq.gz')
    output:
        str(ASSEMBLY_FP/'megahit'/'{sample}_asm'/'final.contigs.fa')
    params:
        out_fp = str(ASSEMBLY_FP/'megahit'/'{sample}_asm')
    threads:
        Cfg['assembly']['threads']
    conda:
        "../../envs/megahit.yml"
    shell:
        """
        ## turn off bash strict mode
        set +o pipefail

        ## sometimes the error is due to lack of memory
        exitcode=0
        ls -al {params.out_fp} > /mnt/d/Penn/sunbeam/log.txt
        if [ -d {params.out_fp} ]
        then
            rm -rf {params.out_fp}
        fi
        megahit -t {threads} -1 {input.r1} -2 {input.r2} -o {params.out_fp} --continue || exitcode=$?

        echo $exitcode
        if [ $exitcode -eq 255 ]
        then
            echo "Empty contigs"
            touch {output}
        elif [ $exitcode -gt 1 ]
        then
            echo "Check your memory"
        fi
        """

rule megahit_unpaired:
    input:
        str(QC_FP/'decontam'/'{sample}_1.fastq.gz')
    output:
        str(ASSEMBLY_FP/'megahit'/'{sample}_asm'/'final.contigs.fa')
    params:
        out_fp = str(ASSEMBLY_FP/'megahit'/'{sample}_asm')
    threads:
        Cfg['assembly']['threads']
    conda:
        "../../envs/megahit.yml"
    shell:
        """
        ## turn off bash strict mode
        set +o pipefail

        ## sometimes the error is due to lack of memory
        exitcode=0
        megahit -t {threads} -r {input} -o {params.out_fp} --continue || exitcode=$?

        if [ $exitcode -eq 255 ]
        then
            echo "Empty contigs"
            touch {output}
        elif [ $exitcode -gt 1 ]
        then
            echo "Check your memory"
        fi
        """

rule final_filter:
    input:
        str(ASSEMBLY_FP/'megahit'/'{sample}_asm'/'final.contigs.fa')
    output:
        str(ASSEMBLY_FP/'contigs'/'{sample}-contigs.fa')
    params:
        len = Cfg['assembly']['min_length']
    log:
        str(ASSEMBLY_FP/'log'/'vsearch'/'{sample}.log')
    conda:
        "../../envs/vsearch.yml"
    script:
        "../../scripts/assembly/final_filter.py"

rule clean_assembly:
    input:
        M = str(ASSEMBLY_FP/'megahit'),
    shell:
        """
        rm -rf {input.M} && echo "Cleanup assembly finished."
        """
