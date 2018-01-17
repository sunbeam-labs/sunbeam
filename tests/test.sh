#!/bin/bash
# Bash flags: Do not commit to repo with these commented out
set -e # Stop on errors
set -x # Display all commands

# Ensure we can activate the environment
export PATH=$PATH:$HOME/miniconda3/bin

# Set up paths
ROOT=`pwd`

if [ $# -ne 1 ]; then
    TEMPDIR=`mktemp -d`
else
    echo "Write sunbeam test result to provided path"
    TEMPDIR=`readlink -f $1`
fi

# Activate the sunbeam environment
source activate sunbeam
command -v snakemake

mkdir -p $TEMPDIR/data_files

function cleanup {
    # Remove temporary directory if it exists
    # (must be careful with rm -rf and variables)
    # Keep retcode from any previous command
    RETCODE=$?
    echo "Exiting with return code $RETCODE"
    [ -z ${TEMPDIR+x} ] || rm -rf "$TEMPDIR"; exit $RETCODE
}

# Calls cleanup when the script exits
if [ $# -ne 1 ]; then
    trap cleanup EXIT
fi

pushd tests
# Copy data into the temporary directory
cp -r indexes $TEMPDIR
cp -r raw $TEMPDIR
cp -r truncated_taxonomy $TEMPDIR
cp seqid2taxid.map $TEMPDIR

python generate_dummy_data.py $TEMPDIR
# Create a version of the config file customized for this tempdir
sunbeam_init $TEMPDIR --defaults testing > $TEMPDIR/tmp_config.yml
popd

pushd $TEMPDIR

# Build fake kraken data
kraken-build --db mindb --add-to-library raw/GCF_Bfragilis_10k_genomic.fna
kraken-build --db mindb --add-to-library raw/GCF_Ecoli_10k_genomic.fna
mv truncated_taxonomy mindb/taxonomy
cp seqid2taxid.map mindb
kraken-build --db mindb --build --kmer-len 16 --minimizer-len 1
kraken-build --db mindb --clean

# Build fake blast database
mkdir -p local/blast
cat raw/*.fna > local/blast/bacteria.fa
makeblastdb -dbtype nucl -in local/blast/bacteria.fa
cp indexes/card.fa local/blast
makeblastdb -dbtype prot -in local/blast/card.fa
popd

# Running snakemake

echo " ===== CONFIG FILE ====="
cat $TEMPDIR/tmp_config.yml
echo " ===== END CONFIG FILE ===== "

# Integration tests: add as needed
# Make each test a function, then call it so that the error can be isolated
# ===================================

# Testing correct expected behavior
function test_all {
echo "Now testing snakemake: "
snakemake --configfile=$TEMPDIR/tmp_config.yml -p
snakemake --configfile=$TEMPDIR/tmp_config.yml clean_assembly -p

# Check contents
echo "Now checking whether we hit the expected genome:"
awk '/NC_000913.3|\t2/  {rc = 1; print}; END { exit !rc }' $TEMPDIR/sunbeam_output/annotation/summary/dummybfragilis.tsv

# Check targets
python tests/find_targets.py --prefix $TEMPDIR/sunbeam_output tests/targets.txt 
}

test_all

# Fix for #38: Make Cutadapt optional
# -- Remove adapter sequences and check to make sure qc proceeds correctly
function test_optional_cutadapt {
sed 's/adapters: \[.*\]/adapters: \[\]/g' $TEMPDIR/tmp_config.yml > $TEMPDIR/tmp_config_nocutadapt.yml
rm -rf $TEMPDIR/sunbeam_output/qc
snakemake --configfile=$TEMPDIR/tmp_config_nocutadapt.yml all_decontam
[ -f $TEMPDIR/sunbeam_output/qc/decontam/dummyecoli_R1.fastq.gz ]
[ -f $TEMPDIR/sunbeam_output/qc/decontam/dummyecoli_R2.fastq.gz ]
}

test_optional_cutadapt

# Test for template option for sunbeamlib: #54
function test_template_option {
pushd tests
# Create a version of the config file customized for this tempdir
# Provide the sunbeamlib package config file manually
CONFIG_FP=$HOME/miniconda3/envs/sunbeam/lib/python3.5/site-packages/sunbeamlib/data/default_config.yml
sunbeam_init $TEMPDIR --template $CONFIG_FP --defaults testing > $TEMPDIR/tmp_config_2.yml
popd
rm -r $TEMPDIR/sunbeam_output
echo "Now re-run snakemake with custom config file: "
snakemake --configfile=$TEMPDIR/tmp_config_2.yml 
snakemake --configfile=$TEMPDIR/tmp_config_2.yml clean_assembly
python tests/find_targets.py --prefix $TEMPDIR/sunbeam_output tests/targets.txt
}

test_template_option


# Test for barcodes file
function test_barcode_file {
sunbeam_mod_config --config $TEMPDIR/tmp_config.yml --mod_str 'all: {samplelist_fp: barcodes.txt}' > $TEMPDIR/tmp_config_barcode.yml
echo -e "dummybfragilis\tTTTTTTTT\ndummyecoli\tTTTTTTTT" > $TEMPDIR/barcodes.txt
rm -rf $TEMPDIR/sunbeam_output/qc/decontam*
echo "CONFIG START"
cat $TEMPDIR/tmp_config_barcode.yml
echo "CONFIG END"
snakemake --configfile=$TEMPDIR/tmp_config_barcode.yml all_decontam
[ -f $TEMPDIR/sunbeam_output/qc/decontam/dummyecoli_R1.fastq.gz ]
[ -f $TEMPDIR/sunbeam_output/qc/decontam/dummyecoli_R2.fastq.gz ]
}

test_barcode_file
