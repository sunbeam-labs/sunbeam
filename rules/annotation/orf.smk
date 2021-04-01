# -*- mode: Snakemake -*-
#
# Contig annotation:
# 	Rules for finding and extracting ORFs.
#
# See Readme.md


rule prodigal:
    """Use Progial for coding genes predictions in contigs."""
    input:
        str(ASSEMBLY_FP/'contigs'/'{sample}-contigs.fa')
    output:
        gff = str(ANNOTATION_FP/'genes'/'prodigal'/'{sample}_genes.gff'),
        faa = str(ANNOTATION_FP/'genes'/'prodigal'/'{sample}_genes_prot.fa'),
        fna = str(ANNOTATION_FP/'genes'/'prodigal'/'{sample}_genes_nucl.fa'),
        log = str(ANNOTATION_FP/'genes'/'prodigal'/'log'/'{sample}.out')
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
        expand(str(ANNOTATION_FP/'genes'/'prodigal'/'{sample}_genes_{suffix}.fa'),
        sample=Samples.keys(), suffix=['prot','nucl'])


        
                
        
