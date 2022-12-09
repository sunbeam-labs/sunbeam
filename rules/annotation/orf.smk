# -*- mode: Snakemake -*-
#
# Contig annotation:
# 	Rules for finding and extracting ORFs.
#
# See Readme.md

rule all_annotation:
    input: TARGET_ANNOTATION

rule prodigal:
    """Use Progial for coding genes predictions in contigs."""
    input:
        ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa",
    output:
        gff=ANNOTATION_FP / "genes" / "prodigal" / "{sample}_genes.gff",
        faa=ANNOTATION_FP / "genes" / "prodigal" / "{sample}_genes_prot.fa",
        fna=ANNOTATION_FP / "genes" / "prodigal" / "{sample}_genes_nucl.fa",
        log=ANNOTATION_FP / "genes" / "prodigal" / "log" / "{sample}.out",
    benchmark:
        BENCHMARK_FP / "prodigal_{sample}.tsv"
    conda:
        "../../envs/annotation.yml"
    shell:
        """
        if [[ -s {input} ]]; then
          prodigal -i {input} -o {output.gff} \
          -a {output.faa} -d {output.fna} -p meta  &> {output.log}
        else
          touch {output.faa}
          touch {output.gff}
          touch {output.fna}
          touch {output.log}
        fi
        """


rule _test_prodigal:
    input:
        expand(
            ANNOTATION_FP / "genes" / "prodigal" / "{sample}_genes_{suffix}.fa",
            sample=Samples.keys(),
            suffix=["prot", "nucl"],
        ),
