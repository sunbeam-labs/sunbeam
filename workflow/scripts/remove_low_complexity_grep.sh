#!/usr/bin/env bash
bash --version
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
zgrep "^@" $read1 > ${log_fp}/${SAMPLEID}.trimm_verbose_ids
sed 's/ .*$//g' ${log_fp}/${SAMPLEID}.trimm_verbose_ids | sed 's/\/[1-2]$//g' | sort -u > ${log_fp}/${SAMPLEID}.trimm_ids
sed 's/ .*$//g' ${ids} | sed 's/\/[1-2]$//g' | sort -u > ${ids}_unique
grep -v -f ${ids}_unique ${log_fp}/${SAMPLEID}.trimm_ids > ${log_fp}/${SAMPLEID}.komplexity_keep_ids
komp_fp="$(dirname "${out1}")"
mkdir -p $komp_fp &>/dev/null # be silent
zgrep -A 3 -f ${log_fp}/${SAMPLEID}.komplexity_keep_ids $read1 | sed '/^--$/d' | gzip > $out1

exitcode=$?

newheaders=$( zgrep -c "^@" $out1 )
newlines=$( zcat $out1 | wc -l )
numids=$(< ${log_fp}/${SAMPLEID}.komplexity_keep_ids wc -l )
explines=$(( "$numids" + "$numids" + "$numids" + "$numids" ))
echo $newheaders &>> $log
echo $numids &>> $log
echo $newlines &>> $log
echo $explines &>> $log
if [ "$newheaders" -eq "$numids" ]; then
	if [ "$newlines" -ne "$explines" ]; then
		exitcode=1	
		echo "Your filtered file does not have the expected length" >> $log
	else
		echo "Your filtered file has the expected length" &>> $log
	fi
fi

if [ $exitcode -eq 1 ]
then
    exit 1
else
    exit 0
fi
