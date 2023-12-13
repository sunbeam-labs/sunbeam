#!/usr/bin/env bash

set +e

read1="${snakemake_input[read1]}"
read2="${snakemake_input[read2]}"
ids="${snakemake_input[ids]}"
out1="${snakemake_output[out1]}"
out2="${snakemake_output[out2]}"
log_fp="$(dirname "${ids}")"
base_name="$(basename "${ids}")"
SAMPLEID=${base_name%.filtered_ids}

echo "make list of trimmomatic output IDs"
zgrep "^@" $read1 > ${log_fp}/${SAMPLEID}.trimm_verbose_ids
zgrep "^@" $read2 >> ${log_fp}/${SAMPLEID}.trimm_verbose_ids
sed 's/ .*$//g' ${log_fp}/${SAMPLEID}.trimm_verbose_ids | sort -u > ${log_fp}/${SAMPLEID}.trimm_ids
echo "grep -v the komplexity ids to get subsample to keep"
grep -v -f ${ids} ${log_fp}/${SAMPLEID}.trimm_ids > ${log_fp}/${SAMPLEID}.komplexity_keep_ids
echo "filter reads with zgrep"
komp_fp="$(dirname "${out1}")"
mkdir -p $komp_fp &>/dev/null # be silent
zgrep -A 3 -f ${log_fp}/${SAMPLEID}.komplexity_keep_ids $read1 | gzip > $out1
zgrep -A 3 -f ${log_fp}/${SAMPLEID}.komplexity_keep_ids $read2 | gzip > $out2

exitcode=$?

if [ $exitcode -eq 1 ]
then
    exit 1
else
    exit 0
fi
