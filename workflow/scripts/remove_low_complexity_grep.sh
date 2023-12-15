#!/usr/bin/env bash

set +e

#read1="${snakemake_input[read1]}"
read1="${snakemake_input[reads]}"
#read2="${snakemake_input[read2]}"
echo $reads

ids="${snakemake_input[ids]}"
out1="${snakemake_output}"
log_fp="$(dirname "${ids}")"
base_name="$(basename "${read1}")"
SAMPLEID=${base_name%.fastq.gz}

echo "make list of trimmomatic output IDs"
zgrep "^@" $read1 > ${log_fp}/${SAMPLEID}.trimm_verbose_ids
sed 's/ .*$//g' ${log_fp}/${SAMPLEID}.trimm_verbose_ids | sed 's/\/[1-2]$//g' | sort -u > ${log_fp}/${SAMPLEID}.trimm_ids
sed 's/ .*$//g' ${ids} | sed 's/\/[1-2]$//g' | sort -u > ${ids}_unique
echo "grep -v the komplexity ids to get subsample to keep"
grep -v -f ${ids}_unique ${log_fp}/${SAMPLEID}.trimm_ids > ${log_fp}/${SAMPLEID}.komplexity_keep_ids
echo "filter reads with zgrep"
komp_fp="$(dirname "${out1}")"
mkdir -p $komp_fp &>/dev/null # be silent
zgrep -A 3 -f ${log_fp}/${SAMPLEID}.komplexity_keep_ids $read1 | sed '/^--$/d' | gzip > $out1

exitcode=$?

newheaders=$( zgrep -c "^@" $out1 )
newlines=$( zcat $out1 | wc -l )
numids=$(< ${log_fp}/${SAMPLEID}.komplexity_keep_ids wc -l )
explines=$(( "$numids" + "$numids" + "$numids" + "$numids" ))
echo $newheaders
echo $newlines
echo $explines
echo $numids
if [ "$newheaders" -eq "$numids" ]; then
	if [ "$newlines" -eq "$explines" ]; then
		continue
	fi
else	
	exitcode=$(( "$exitcode" + 1 ))	
	echo "Your filtered read 1 file does not equal the expected length"
fi

if [ $exitcode -eq 1 ]
then
    exit 1
else
    exit 0
fi
