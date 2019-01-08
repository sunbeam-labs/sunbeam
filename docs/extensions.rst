.. _extensions:

==============
Sunbeam Extensions
==============

Sunbeam extensions allow you to add new features to the pipeline,
which can be re-distributed to other researchers to facilitate
reproducible analyses.

Input files for extensions:

+-----------------------+----------------------------------------------------------------+
| Sequence data files   | Target                                                         |
+=======================+================================================================+
| Quality-controlled,   | str(QC_FP/'cleaned'/'{sample}_{rp}.fastq.gz')                  |
| non-decontaminated    |                                                                |
| sequences             |                                                                |
+-----------------------+----------------------------------------------------------------+
| Quality-controlled,   | str(QC_FP/'decontam'/'{sample}_{rp}.fastq.gz')                 |
| decontaminated        |                                                                |
| sequences             |                                                                |
+-----------------------+----------------------------------------------------------------+
| Contig sequences      | str(ASSEMBLY_FP/'contigs'/'{sample}-contigs.fa')               |
+-----------------------+----------------------------------------------------------------+
| Open reading frame    | str(ANNOTATION_FP/'genes'/'prodigal'/'{sample}_genes_nucl.fa') |
| nucleotide sequences  |                                                                |
+-----------------------+----------------------------------------------------------------+
| Open reading frame    | str(ANNOTATION_FP/'genes'/'prodigal'/'{sample}_genes_prot.fa') |
| protein sequences     |                                                                |
+-----------------------+----------------------------------------------------------------+


+-----------------------+-----------------------------------------------+
| Summary tables        | Target                                        |
+=======================+===============================================+
| Attrition from        | str(QC_FP/'reports'/'preprocess_summary.tsv') |
| decontamination and   |                                               |
| quality control       |                                               |
+-----------------------+-----------------------------------------------+
| Sequence              | str(QC_FP/'reports'/'fastqc_quality.tsv')     |
| quality scores        |                                               |
+-----------------------+-----------------------------------------------+
| Taxonomic assignments | str(CLASSIFY_FP/'kraken'/'all_samples.tsv')   |
| from Kraken           |                                               |
+-----------------------+-----------------------------------------------+
