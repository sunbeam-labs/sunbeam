# Test normal behavior
function test_all {
    sunbeam run -- --configfile=$TEMPDIR/tmp_config.yml -p

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
    sunbeam run --configfile=$TEMPDIR/tmp_config_nocutadapt.yml all_decontam
    [ -f $TEMPDIR/sunbeam_output/qc/decontam/dummyecoli_1.fastq.gz ]
    [ -f $TEMPDIR/sunbeam_output/qc/decontam/dummyecoli_2.fastq.gz ]
}

# Fix for #54
# Test for template option for sunbeamlib
function test_template_option {
    pushd tests
    CONFIG_FP=cfg_template.yml
    # Create a version of the config file customized for this tempdir
    # Provide the sunbeamlib package config file manually
    sunbeam init --force --output 54.yml --template $CONFIG_FP $TEMPDIR
    grep 'from_template:' $TEMPDIR/54.yml || exit 1
    popd
}


# Test for version check
function test_version_check {
    sunbeam config modify --str 'all: {version: 9999.9.9}' \
	    $TEMPDIR/tmp_config.yml > $TEMPDIR/too_high_config.yml
    # this should produce a nonzero exit code and fail if it does not
    if sunbeam run --configfile $TEMPDIR/too_high_config.yml; then
	exit 1
    fi
}

# Test that we detect and run extensions
function test_extensions {
    sunbeam run --configfile $TEMPDIR/tmp_config.yml sbx_test | grep "SBX_TEST"
}

# Test that single-end sequencing configurations work
function test_single_end {
    rm -rf $TEMPDIR/sunbeam_output/qc
    sunbeam config modify --str 'all: {paired_end: false}' \
	    $TEMPDIR/tmp_config.yml > $TEMPDIR/single_end_config.yml
    sunbeam run --configfile $TEMPDIR/single_end_config.yml
    python tests/find_targets.py --prefix $TEMPDIR/sunbeam_output tests/targets_singleend.txt
}

# Fix for #131
# Test that paired-end qc rules produce files with the same number of reads
function test_pair_concordance {
    rm -rf $TEMPDIR/sunbeam_output/qc
    sunbeam run --configfile $TEMPDIR/tmp_config.yml all_decontam
    for r1 in $TEMPDIR/sunbeam_output/qc/cleaned/*_1.fastq.gz; do
	r1_lines=$(zcat $r1 | wc -l)
	r2=${r1%_1.fastq.gz}_2.fastq.gz
	r2_lines=$(zcat $r2 | wc -l)
	if [ $r1_lines -ne $r2_lines ]; then
	    exit 1
	fi
    done
}
