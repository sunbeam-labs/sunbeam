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

# Test that we can guess a variety of sample names correctly
# Correct behavior for only two samples
function test_guess_with_two_samples {
    mkdir -p $TEMPDIR/only_two_samples
    touch $TEMPDIR/only_two_samples/sample_1.fastq.gz
    touch $TEMPDIR/only_two_samples/sample_2.fastq.gz
    sunbeam list_samples $TEMPDIR/only_two_samples 2> >(tee out.txt >&2)
    grep '{sample}_{rp}.fastq.gz' out.txt
}

# Correct behavior for samples with inconsistent _ or .
function test_guess_with_inconsistent_samples {
    mkdir -p $TEMPDIR/inconsistent_samples
    touch $TEMPDIR/inconsistent_samples/asdf_123_R1.fastq.gz
    touch $TEMPDIR/inconsistent_samples/asdf_123_R2.fastq.gz
    touch $TEMPDIR/inconsistent_samples/asddf_R1.fastq.gz
    touch $TEMPDIR/inconsistent_samples/asddf_R2.fastq.gz
    sunbeam list_samples $TEMPDIR/inconsistent_samples 2> >(tee out.txt >&2)
    grep '{sample}_R{rp}.fastq.gz' out.txt
    rm -r $TEMPDIR/inconsistent_samples
}

# Correct behavior for folders that still have index files
function test_guess_with_index_files_present {
    mkdir -p $TEMPDIR/idx_files_present
    touch $TEMPDIR/idx_files_present/asdf_123_R1.fastq.gz
    touch $TEMPDIR/idx_files_present/asdf_123_R2.fastq.gz
    touch $TEMPDIR/idx_files_present/asddf_R1.fastq.gz
    touch $TEMPDIR/idx_files_present/asddf_R2.fastq.gz
    touch $TEMPDIR/idx_files_present/asdf_123_I1.fastq.gz
    touch $TEMPDIR/idx_files_present/asdf_123_I2.fastq.gz
    touch $TEMPDIR/idx_files_present/asddf_I1.fastq.gz
    touch $TEMPDIR/idx_files_present/asddf_I2.fastq.gz
    sunbeam list_samples $TEMPDIR/idx_files_present 2> >(tee out.txt >&2)
    grep '{sample}_R{rp}.fastq.gz' out.txt
#    rm -r $TEMPDIR/idx_files_present
}

# Fix for #153
# Empty strings for _fp entries in the config should result in the relevant
# files being ignored, even if they would normally match the pattern used.
function test_blank_fp_behavior {
    # Move indexes to top-level
    for file in $TEMPDIR/indexes/*; do
        mv $file $TEMPDIR/indexes_${file##*/}
    done
    # Run sunbeam with a blank string for the mapping step's genomes path, and
    # with fasta files lying around in the top-level directory.
    sunbeam config modify --str 'mapping: {genomes_fp: ""}' \
        $TEMPDIR/tmp_config.yml > $TEMPDIR/blank_fp_config.yml
    sunbeam run --configfile $TEMPDIR/blank_fp_config.yml
    # Move indexes back to original location
    for file in $TEMPDIR/indexes_*; do
        mv $file ${file/indexes_/indexes\//}
    done
    # These files should *not* exist, because we don't accept the top-level
    # directory as a genomes_fp directory.
    [ ! -f $TEMPDIR/sunbeam_output/mapping/indexes_human/coverage.csv ]
    [ ! -f $TEMPDIR/sunbeam_output/mapping/indexes_phix174/coverage.csv ]
}

# Test that mapping reports aligning read pairs for a contrived example.
# This could maybe be rolled into test_all (as the annotation check already is)
# if the data setup were reorganized.
function test_mapping {
    # Create two read pairs using lines from the human genome fasta.
    r1_1=$TEMPDIR/data_files/PCMP_stub_human_R1.fastq.gz
    r2_1=$TEMPDIR/data_files/PCMP_stub_human_R2.fastq.gz
    r1_2=$TEMPDIR/data_files/PCMP_stub2_human_R1.fastq.gz
    r2_2=$TEMPDIR/data_files/PCMP_stub2_human_R2.fastq.gz
    human=$TEMPDIR/indexes/human.fasta
    (
        echo "@read0"
        sed -n 2p $human
        echo "+"
        sed -n 's:.:G:g;2p' $human
    ) | gzip > $r1_1
    (
        echo "@read0"
        sed -n 2p $human | rev | tr '[ACTG]' '[TGAC]'
        echo "+"
        sed -n 's:.:G:g;2p' $human
    ) | gzip > $r2_1
    (
        echo "@read0"
        sed -n 15p $human | cut -c 1-40
        echo "+"
        sed -n 's:.:G:g;2p' $human | cut -c 1-40
    ) | gzip > $r1_2
    (
        echo "@read0"
        sed -n 15p $human | cut -c 1-40 | rev | tr '[ACTG]' '[TGAC]'
        echo "+"
        sed -n 's:.:G:g;2p' $human | cut -c 1-40
    ) | gzip > $r2_2
    # Run sunbeam mapping rules with these two samples defined.
    (
	    echo "stub_human,$r1_1,$r2_1"
	    echo "stub2_human,$r1_2,$r2_2"
    ) > $TEMPDIR/samples_test_mapping.csv
    sunbeam config modify --str 'all: {samplelist_fp: "samples_test_mapping.csv"}' \
        $TEMPDIR/tmp_config.yml > $TEMPDIR/test_mapping_config.yml
    # There should be two lines in the human coverage summary and none at all
    # in the phix174 summary.  The human.csv lines should be sorted in standard
    # alphanumeric order; stub2_human will come before stub_human.
    sunbeam run --configfile $TEMPDIR/test_mapping_config.yml all_mapping
    md5sum --check --status <(
    echo "c624406eb2582cac5e0cfb160c79a900  $TEMPDIR/sunbeam_output/mapping/human/coverage.csv"
    echo "1aee435ade0310a6b3c63d44cbdc2029  $TEMPDIR/sunbeam_output/mapping/phix174/coverage.csv"
    )
}

# Test for sunbeam init for SRA
# Make sure samples.csv contains the correct number of samples.

function test_sunbeam_get {
    mkdir -p $TEMPDIR/test_sunbeam_get
    sunbeam init --force --output sunbeam_config_SRA.yml $TEMPDIR/test_sunbeam_get --data_acc SRP021545
    test `wc -l < $TEMPDIR/test_sunbeam_get/samples.csv` -eq 89
}

# Test for sunbeam init for SRA -- study with paired and unpaired samples
# Make sure Sunbeam exits with nonzero exit code if a study contains paired and unpaired reads
# Both sets should be written separately to a config/samples.csv pair of files

function test_get_paired_unpaired {
    dp=$TEMPDIR/test_get_paired_unpaired
    mkdir -p $dp
    # "!" because we *expect* this to exit nonzero.
    ! sunbeam init --force --output sunbeam_config_SRA.yml $dp --data_acc ERP020555
    # Check contents of the two config files
    grep '^  samplelist_fp: samples_unpaired.csv$' $dp/unpaired_sunbeam_config_SRA.yml
    grep '^  paired_end: false$'                   $dp/unpaired_sunbeam_config_SRA.yml
    grep '^  samplelist_fp: samples_paired.csv$'   $dp/paired_sunbeam_config_SRA.yml
    grep '^  paired_end: true$'                    $dp/paired_sunbeam_config_SRA.yml
    # Check contents of the two samples csv files
    test `wc -l < $dp/samples_unpaired.csv` -eq  1
    test `wc -l < $dp/samples_paired.csv`   -eq  13
}
