# Test normal behavior
function test_all {
    snakemake --configfile=$TEMPDIR/tmp_config.yml -p

    # Check contents
    awk '/NC_000913.3|\t2/  {rc = 1; print}; END { exit !rc }' $TEMPDIR/sunbeam_output/annotation/summary/dummyecoli.tsv

    # Check targets
    python tests/find_targets.py --prefix $TEMPDIR/sunbeam_output tests/targets.txt 
}


# Fix for #38: Make Cutadapt optional
# Remove adapter sequences and check to make sure qc proceeds correctly
function test_optional_cutadapt {
    sed 's/adapters: \[.*\]/adapters: \[\]/g' $TEMPDIR/tmp_config.yml > $TEMPDIR/tmp_config_nocutadapt.yml
    rm -rf $TEMPDIR/sunbeam_output/qc
    snakemake --configfile=$TEMPDIR/tmp_config_nocutadapt.yml all_decontam
    [ -f $TEMPDIR/sunbeam_output/qc/decontam/dummyecoli_R1.fastq.gz ]
    [ -f $TEMPDIR/sunbeam_output/qc/decontam/dummyecoli_R2.fastq.gz ]
}

# Fix for #54
# Test for template option for sunbeamlib
function test_template_option {
    pushd tests
    # Create a version of the config file customized for this tempdir
    # Provide the sunbeamlib package config file manually
    CONFIG_FP=cfg_template.yml
    sunbeam init $TEMPDIR --template $CONFIG_FP --defaults testing | grep 'from_template:' || exit 1
    popd
}


# Test for barcodes file
function test_barcode_file {
    sunbeam_mod_config --config $TEMPDIR/tmp_config.yml --mod_str 'all: {samplelist_fp: barcodes.txt}' > $TEMPDIR/tmp_config_barcode.yml
    echo -e "dummybfragilis\tTTTTTTTT\ndummyecoli\tTTTTTTTT" > $TEMPDIR/barcodes.txt
    snakemake --configfile=$TEMPDIR/tmp_config_barcode.yml all_decontam
    [ -f $TEMPDIR/sunbeam_output/qc/decontam/dummyecoli_R1.fastq.gz ]
    [ -f $TEMPDIR/sunbeam_output/qc/decontam/dummyecoli_R2.fastq.gz ]
}


# Test for version check
function test_version_check {
    sunbeam_mod_config --config $TEMPDIR/tmp_config.yml --mod_str 'all: {version: 9999.9.9}' > $TEMPDIR/too_high_config.yml
    # this should produce a nonzero exit code and fail if it does not
    if snakemake --configfile $TEMPDIR/too_high_config.yml; then
	exit 1
    fi
}

# Test that we can 
function test_extensions {
    snakemake --configfile $TEMPDIR/tmp_config.yml sbx_test | grep "SBX_TEST"
}
