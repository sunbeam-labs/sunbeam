# -*- mode: Snakemake -*-
#
# Contig annotation.
#
# See Readme.md


rule all_annotate:
    input:
        TARGET_ANNOTATE,


rule aggregate_results:
    input:
        contigs=ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa",
        contig_results=expand(
            ANNOTATION_FP / "blastn" / "{db}" / "contig" / "{{sample}}.btf",
            db=Blastdbs["nucl"],
        ),
        gene_results=expand(
            ANNOTATION_FP / "{blastpx}" / "{db}" / "{orf_finder}" / "{{sample}}.btf",
            blastpx=["blastp", "blastx"],
            db=Blastdbs["prot"],
            orf_finder=["prodigal"],
        ),
    output:
        ANNOTATION_FP / "summary" / "{sample}.tsv",
    benchmark:
        BENCHMARK_FP / "aggregate_results_{sample}.tsv"
    params:
        dbs=list(Blastdbs["nucl"].keys()) + list(Blastdbs["prot"].keys()),
        nucl=Blastdbs["nucl"],
        prot=Blastdbs["prot"],
    conda:
        "../../envs/annotation.yml"
    script:
        "../../scripts/annotation/aggregate_results.py"


rule aggregate_all:
    input:
        expand(ANNOTATION_FP / "summary" / "{sample}.tsv", sample=Samples.keys()),
    output:
        ANNOTATION_FP / "all_samples.tsv",
    run:
        with open(output[0], "w") as out:
            out.writelines(open(input[0]).readlines())
            for infile in input[1:]:
                out.writelines(l for i, l in enumerate(open(infile)) if i > 0)
