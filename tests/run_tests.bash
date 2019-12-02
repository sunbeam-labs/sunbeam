#!/bin/bash

# setup

set -e

STARTING_DIR=$(pwd)

# Ensure we're running in the correct directory
case $BASH_SOURCE in
    tests/*)
	;;
    *run_tests.bash)
	echo "Error: must be run from sunbeam directory using tests/run_test.sh"
	exit 1
	;;
    *)
	echo "Unrecognized script path, may cause errors..."
	;;
esac

# Colors
GREEN="\x1B[32m"
RED="\x1B[31m"
RESET="\x1B[0m"

PASS="${GREEN}\u2714${RESET}"
FAIL="${RED}x${RESET}"

# All testing info will be output here
TTY=${TTY:-/dev/tty}

# Testing options
USE_TMPDIR=true
USE_TMPENV=true
SBX_FP=extensions
VERBOSE=false

while getopts "d:e:t:vh" opt; do
    case $opt in
	d)
	    USE_TMPDIR=false
	    TEMPDIR=`readlink -f $OPTARG`
	    mkdir -p $TEMPDIR
	    ;;
	e)
	    USE_TMPENV=false
	    SUNBEAM_ENV=$OPTARG
	    ;;
	t)
	    RUN_TEST=$OPTARG
	    ;;
	v)
	    VERBOSE=true
	    ;;
	h)
	    echo "Run the Sunbeam test suite."
	    echo "  -d DIR       Use DIR rather than a temporary directory (remains after tests finish)"
	    echo "  -e ENV_NAME  Use a pre-existing Conda environment rather than creating one (remains after tests finish)"
	    echo "  -t TEST      Run a specific test from tests/test_suite.bash only"
	    echo "  -v           Show command output while running"
	    echo "  -h           Display this message and exit"
	    exit 1
	    ;;
	\?)
	    echo "Unrecognized option -'$OPTARG'"
	    exit 1
	    ;;
    esac
done

function verbose {
    if [ "$VERBOSE" = true ]; then
	echo -ne "${1}"
    fi
} &>$TTY

function msg {
    echo -ne "${1}"
} &>$TTY

# debug
verbose "use_tmpdir: $USE_TMPDIR\n"
verbose "use_tmpenv: $USE_TMPENV\n"

function broke {
    local RETCODE=$?
    msg "\nFailed command error output:\n`cat ${2}.err`\n"
    msg "${FAIL} (log: ${LOGFILE}.[out/err])\n"
    cleanup 1
} 

function capture_output {
    msg "Running ${1}... "
    if [ -z ${TEMPDIR+x} ]; then
	LOGFILE="${1}"
    else
	LOGFILE="${TEMPDIR}/${1}"	
    fi

    set -o pipefail
    if [ "$VERBOSE" = true ]; then
	OUTPUT_STRING="> >(tee ${LOGFILE}.out) 2> >(tee ${LOGFILE}.err >&2)"
    else
	OUTPUT_STRING="> ${LOGFILE}.out 2> ${LOGFILE}.err"
    fi
    trap "broke ${1} ${LOGFILE} $?" exit
    eval "${1} ${OUTPUT_STRING}"
    set +o pipefail
    trap "cleanup $?" exit
    msg "${PASS}\n"
}

function setup {
    # Create temporary directory
    if [ "$USE_TMPDIR" = true ]; then
	TEMPDIR=`mktemp -d`
    fi
    
    verbose "\n\t${GREEN}Test directory${RESET}: ${TEMPDIR}"

    export PATH=$PATH:$HOME/miniconda3/bin
    # Allow conda [de]activate
    CONDA_BASE=$(conda info --base) # see https://github.com/conda/conda/issues/7980
    source $CONDA_BASE/etc/profile.d/conda.sh

    # Install Sunbeam (maybe)
    if [ "$USE_TMPENV" = true ]; then
	SUNBEAM_ENV="sunbeam-`basename $TEMPDIR`"
	bash install.sh -e $SUNBEAM_ENV 
    fi
    
    verbose "\n\t${GREEN}Conda environment${RESET}: ${SUNBEAM_ENV}\n"

    # Activate Sunbeam
    conda activate $SUNBEAM_ENV

    # Move extensions out of the way temporarily
    if [ -d $SBX_FP ]; then
	OLD_SBX_FP="${SBX_FP}_moved_for_testing"
	mv $SBX_FP $OLD_SBX_FP
    fi
    mkdir $SBX_FP
}
    
function cleanup {
    local TMPRC=$?
    local RETCODE=$TMPRC
    if [ ${1} -gt ${TMPRC} ]; then
	RETCODE=${1}
    else
	RETCODE=${TMPRC}
    fi
    cd $STARTING_DIR
    if [ $RETCODE -ne 0 ]; then
	msg "${RED}-- TESTS FAILED --${RESET}\n"
    else
	msg "${GREEN}-- TESTS SUCCEEDED --${RESET}\n"
	if [ ! -z ${LOGFILE+x} ]; then
	    [ -f "${LOGFILE}.err" ] && rm "${LOGFILE}.err"
	    [ -f "${LOGFILE}.out" ] && rm "${LOGFILE}.out"
	fi
    fi
    conda deactivate
    # Remove Sunbeam environment if created
    if [ "$INSTALL_SUNBEAM" = true ]; then
	verbose "Deleting temporary Sunbeam environment ${SUNBEAM_ENV} \n"
	conda env remove -yn $SUNBEAM_ENV
    fi
    # Remove temp directory if created
    if [ "$USE_TMPDIR" = true ]; then
	if [ ! -z ${TEMPDIR+x} ]; then
	    verbose "Deleting temporary directory ${TEMPDIR}\n"
	    rm -rf $TEMPDIR
	fi
    fi
    # Restore old extension folder
    if [ ! -z ${OLD_SBX_FP+x} ]; then
	verbose "Restoring extensions folder\n"
	rm -rf "$SBX_FP" && mv "$OLD_SBX_FP" "$SBX_FP"
    fi
    # Exit, maintaining previous return code
    exit $RETCODE
}

function build_test_data {
    # Copy data into the temporary directory
    pushd tests
    cp -r indexes $TEMPDIR
    cp -r raw $TEMPDIR
    cp -r truncated_taxonomy $TEMPDIR
    cp -r sbx_test $STARTING_DIR/$SBX_FP/sbx_test
    cp -r sbx_test_subdir $STARTING_DIR/$SBX_FP/sbx_test_subdir
    cp seqid2taxid.map $TEMPDIR

    mkdir -p $TEMPDIR/hosts
    cp indexes/*.fasta $TEMPDIR/hosts
    python generate_dummy_data.py $TEMPDIR
    python generate_dummy_data.py $TEMPDIR data_files_old_illumina /1 /2
    # Create a version of the config file customized for this tempdir
    sunbeam init \
	    --force \
	    --output tmp_config.yml \
	    --defaults test_config_defaults.yml \
	    --data_fp $TEMPDIR/data_files \
	    $TEMPDIR
    popd
    
    # Build fake kraken data
    pushd $TEMPDIR
    kraken2-build --db mindb --add-to-library \
		 raw/GCF_Bfragilis_10k_genomic.fna
    kraken2-build --db mindb --add-to-library \
		 raw/GCF_Ecoli_10k_genomic.fna
    mv truncated_taxonomy mindb/taxonomy
    cp seqid2taxid.map mindb
    kraken2-build --db mindb --build --kmer-len 16 --minimizer-len 1 --minimizer-spaces 0
    kraken2-build --db mindb --clean

    # Build fake blast database
    mkdir -p local/blast
    cat raw/*.fna > local/blast/bacteria.fa
    makeblastdb -dbtype nucl -in local/blast/bacteria.fa
    cp indexes/card.fa local/blast 
    makeblastdb -dbtype prot -in local/blast/card.fa
    popd
}

function run_test_suite {
    for testcase in $(declare -f | grep -o "^test[a-zA-Z_]*") ; do
	capture_output ${testcase}
    done
}

trap cleanup exit

capture_output setup
capture_output build_test_data

source tests/test_suite.bash

# Run single test, if specified, or all detected tests otherwise
if [ ! -z ${RUN_TEST+x} ]; then
    capture_output ${RUN_TEST}
else
    run_test_suite
fi



