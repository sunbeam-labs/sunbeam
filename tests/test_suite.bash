# Test normal behavior
function test_all {
    sunbeam run -- --configfile=$TEMPDIR/tmp_config.yml -p

    # Check contents
    annot_summary=sunbeam_output/annotation/summary/dummyecoli.tsv
    awk '/NC_000913.3|\t2/  {rc = 1; print}; END { exit !rc }' $TEMPDIR/$annot_summary || (
        # stderr will show up on the summary output for the test suite as well
        # as in the .err file.  false will cause the shell to exit assuming -e
        # is in effect here.
        echo "Check failed on $annot_summary" > /dev/stderr
        false
    )

    # Check targets
    python tests/find_targets.py --prefix $TEMPDIR/sunbeam_output tests/targets.txt 
}

# For #221: full test using old-style Illumina paired sequence IDs (/1 and /2)
function test_all_old_illumina {
    # init will always write to samples.csv so we'll stash the old one and then
    # restore it.
    mv $TEMPDIR/samples.csv $TEMPDIR/samples_orig.csv
    sunbeam init \
            --force \
            --output tmp_config_old_illumina.yml \
            --defaults <(sed s/sunbeam_output/sunbeam_output_old_illumina/ $TEMPDIR/tmp_config.yml) \
            --data_fp $TEMPDIR/data_files_old_illumina \
            $TEMPDIR
    # Add config entry for suffix to remove from sequence IDs, and run just
    # like before.
    sed -i -r 's/^( +seq_id_ending: *).*$/\1"\/[12]"/' $TEMPDIR/tmp_config_old_illumina.yml
    sunbeam run -- --configfile=$TEMPDIR/tmp_config_old_illumina.yml -p
    mv $TEMPDIR/samples_orig.csv $TEMPDIR/samples.csv

    # Check contents
    annot_summary=sunbeam_output_old_illumina/annotation/summary/dummyecoli.tsv
    awk '/NC_000913.3|\t2/  {rc = 1; print}; END { exit !rc }' $TEMPDIR/$annot_summary || (
        # stderr will show up on the summary output for the test suite as well
        # as in the .err file.  false will cause the shell to exit assuming -e
        # is in effect here.
        echo "Check failed on $annot_summary" > /dev/stderr
        false
    )

    # Check targets
    python tests/find_targets.py --prefix $TEMPDIR/sunbeam_output_old_illumina tests/targets.txt
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

# Test that we have the updated snakemake that uses "conda activate"
function test_use_conda {
    sunbeam run --configfile $TEMPDIR/tmp_config.yml --use-conda sbx_test | grep "SBX_TEST"
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
    # Move human host files to top-level, since we're using that genome for
    # mapping in this test and shouldn't decontaminate using it as well.
    for file in $TEMPDIR/hosts/human*; do
        mv $file $TEMPDIR/hosts_${file##*/}
    done
    sunbeam run --configfile $TEMPDIR/test_mapping_config.yml all_mapping
    # Move human host files back to original location
    for file in $TEMPDIR/hosts_*; do
        mv $file ${file/hosts_/hosts\//}
    done
    # After the header line, there should be two lines in the human and phix
    # coverage summaries, with two reads mapping for human and none for phix.
    # The lines should be sorted in standard alphanumeric order; stub2_human
    # will come before stub_human.
    (
	    csv_human=$TEMPDIR/sunbeam_output/mapping/human/coverage.csv
	    csv_phix=$TEMPDIR/sunbeam_output/mapping/phix174/coverage.csv
	    function col3 { cut -f3 -d, | tr '\n' : ; }
	    function col5 { cut -f5 -d, | tr '\n' : ; }
	    test "Sample:stub2_human:stub_human:" == $(col3 < "$csv_human")
	    test "Max:2:2:" == $(col5 < "$csv_human")
	    test "Sample:stub2_human:stub_human:" == $(col3 < "$csv_phix")
	    test "Max:0:0:" == $(col5 < "$csv_phix")
    ) || (
	    echo "Unexpected coverage.csv content from mapping rules" > /dev/stderr
	    false
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

# Fix for #185:
# While the core Sunbeam rules keep a simple directory structure, using more
# complicated nested subdirectories can complicate the wildcard/graph
# resolution in Snakemake and result in unexpected wildcard patterns being
# evaluated.  Enforcing a pattern on our sample names (no slashes allowed)
# avoids this.
function test_subdir_patterns {
    # All we need to check is that the graph resolution works.
    sunbeam run --configfile $TEMPDIR/tmp_config.yml sbx_test_subdir -n
}

# Fix for #167:
# Check that if megahit gives a nonzero exit code it is handled appropriately.
# The two main cases are 255 (empty contigs) and anything else nonzero
# (presumed to be memory-related in the assembly rules).
# Checking for successful behavior is already handled in test_all.
function test_assembly_failures {
    # Up to just before the assembly rules, things should work fine.
    sunbeam run -- --configfile=$TEMPDIR/tmp_config.yml -p all_decontam
    # Remove previous assembly files, if they exist.
    rm -rf $TEMPDIR/sunbeam_output/assembly
    # If megahit gives an exit code != 0 and != 255 it is an error.
    mkdir -p "$TEMPDIR/megahit_137"
    echo -e '#!/usr/bin/env bash\nexit 137' > $TEMPDIR/megahit_137/megahit
    chmod +x $TEMPDIR/megahit_137/megahit
    (
    export PATH="$TEMPDIR/megahit_137:$PATH"
    # (This command should *not* exit successfully.)
    ! txt=$(sunbeam run -- --configfile=$TEMPDIR/tmp_config.yml -p all_assembly)
    echo "$txt" | grep "Check your memory"
    )
    # If megahit exits with 255, it implies no contigs were built.
    mkdir -p "$TEMPDIR/megahit_255"
    echo -e '#!/usr/bin/env bash\nexit 255' > $TEMPDIR/megahit_255/megahit
    chmod +x $TEMPDIR/megahit_255/megahit
    (
    export PATH="$TEMPDIR/megahit_255:$PATH"
    txt=$(sunbeam run -- --configfile=$TEMPDIR/tmp_config.yml -p all_assembly)
    echo "$txt" | grep "Empty contigs"
    )
}

# For #150 and #152: make sure sunbeam config update works
function test_sunbeam_config_update {
    # Make a new config file
    cp $TEMPDIR/tmp_config.yml $TEMPDIR/tmp_config_test_config_update.yml

    # Break config
    sed -i 's/download_reads: false:/download_reads: BROKEN/' $TEMPDIR/tmp_config_test_config_update.yml
    sunbeam config update $TEMPDIR/tmp_config_test_config_update.yml
    test `cat $TEMPDIR/tmp_config_test_config_update.yml | grep "BROKEN" | wc -l` -eq 0

}

# For #247: test to see whether extension config is included in the main configfile on initialization
function test_extension_config_init {

    if [ `echo $HOME | grep "/home/circleci" | wc -l` -eq 1 ]; then
        echo "Tests running on CircleCI, adding config for sbx_test"
        echo "sbx_test:" > $SUNBEAM_DIR/extensions/sbx_test/config.yml
        echo "  test_param: ''" >> $SUNBEAM_DIR/extensions/sbx_test/config.yml
    fi

    sunbeam init \
            --force \
            --output tmp_config_inclsbx.yml \
            --data_fp $TEMPDIR/data_files \
            $TEMPDIR

    echo "sbx_test found in config" `cat $TEMPDIR/tmp_config_inclsbx.yml | grep "sbx_test:" | wc -l` "time(s)" 
    test `cat $TEMPDIR/tmp_config_inclsbx.yml | grep "sbx_test:" | wc -l` -eq 1

}

# For #247: make sure `sunbeam config update` includes extension info
function test_extension_config_update {

    # Make a config file copy
    cp $TEMPDIR/tmp_config_inclsbx.yml $TEMPDIR/tmp_config_extension_config_update.yml
    cat $TEMPDIR/tmp_config_extension_config_update.yml | grep "test"

    # Remove extension config entry
    sed -i 's/sbx_test://' $TEMPDIR/tmp_config_extension_config_update.yml
    sed -i "s/  test_param: ''//" $TEMPDIR/tmp_config_extension_config_update.yml
    echo "edited config instances of 'test':" `cat $TEMPDIR/tmp_config_extension_config_update.yml | grep "test" | wc -l`
    # Update it again
    sunbeam config update -i $TEMPDIR/tmp_config_extension_config_update.yml
    echo "sbx_test found in config" `cat $TEMPDIR/tmp_config_extension_config_update.yml | grep "sbx_test:"` "time(s)"
    test `cat $TEMPDIR/tmp_config_extension_config_update.yml | grep "sbx_test:" | wc -l` -eq 1
}

# For #251: test sunbeam extend

function test_all_sunbeam_extend {
    sunbeam extend https://github.com/sunbeam-labs/sbx_coassembly
    sunbeam config update -i $TEMPDIR/tmp_config.yml
    sunbeam run --use-conda --configfile=$TEMPDIR/tmp_config.yml -p all_coassemble
    test `ls $TEMPDIR/sunbeam_output/assembly | grep "coassembly" | wc -l` -eq 1
}

# For #261: handle URLs with a trailing slash

function test_extend_trailing_slash {

    sunbeam extend https://github.com/sunbeam-labs/sbx_metaphlan/

    rm -rf $SUNBEAM_DIR/extensions/sbx_metaphlan/

}
