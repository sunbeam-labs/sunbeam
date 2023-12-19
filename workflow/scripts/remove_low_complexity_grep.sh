#!/usr/bin/env bash

set +e

read1="${snakemake_input[reads]}"
log="${snakemake_log}"
ids="${snakemake_input[ids]}"
out1="${snakemake_output}"
log_fp="$(dirname "${ids}")"
base_name="$(basename "${read1}")"
SAMPLEID=${base_name%.fastq.gz}
#NAME=${base_name%_[1-2].fastq.gz}

#echo $NAME
echo $read1 &>> $log
echo $SAMPLEID &>> $log
zgrep -F "^@" $read1 > ${log_fp}/${SAMPLEID}.trimm_verbose_ids
sed 's/ .*$//g' ${log_fp}/${SAMPLEID}.trimm_verbose_ids | sed 's/\/[1-2]$//g' | sort -u > ${log_fp}/${SAMPLEID}.trimm_ids
sed 's/ .*$//g' ${ids} | sed 's/\/[1-2]$//g' | sort -u > ${ids}_unique
grep -vF -f ${ids}_unique ${log_fp}/${SAMPLEID}.trimm_ids > ${log_fp}/${SAMPLEID}.komplexity_keep_ids
komp_fp="$(dirname "${out1}")"
mkdir -p $komp_fp &>/dev/null # be silent
zgrep -F -A 3 -f ${log_fp}/${SAMPLEID}.komplexity_keep_ids $read1 | sed '/^--$/d' | gzip > $out1
mistakes=$( zgrep -F -c -f -m1 $ids $out1 ) 2>/dev/null
newheaders=$( zgrep -cF "^@" $out1 )
newlines=$( zcat $out1 | wc -l )
numids=$(< ${log_fp}/${SAMPLEID}.komplexity_keep_ids wc -l )
explines=$(( "$numids" + "$numids" + "$numids" + "$numids" ))
echo $newheaders &>> $log
echo $numids &>> $log
echo $newlines &>> $log
echo $explines &>> $log

if [ "$mistakes" -gt 0 ]; then
	echo "ERROR: Your filtered fastq file contains illegal reads. Try increasing memory and threads" &>> $log
	exit 1
fi

if [ "$newheaders" -ne "$numids" ]; then
	echo "ERROR: Your filtered list of IDs does not have the expected length" &>> $log
	exit 1
fi
if [ "$newlines" -ne "$explines" ]; then
	echo "ERROR: Your filtered fastq file does not have the expected length" &>> $log
	exit 1
else
	echo "Your filtered list of IDs and output fastq file have the expected length" &>> $log
fi

exitcode=$?
if [ $exitcode -eq 1 ]
then
    exit 1
else
    exit 0
fi
