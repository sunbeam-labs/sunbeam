#!/bin/bash
set -e

# Colors
GREEN="\x1B[32m"
RED="\x1B[31m"
RESET="\x1B[0m"

function truefalse {
    if [ $1 = true ]; then
	echo -ne "\u2713"
    else
	echo -ne "\u2717"
    fi
}
	
SCRIPT_DIR=$(dirname $(readlink -f $BASH_SOURCE))
UPDATE=false

function msg {
    echo -ne "${1}"
} &>$TTY


while getopts "e:s:c:u:qh" opt; do
    case $opt in
	e)
	    SUNBEAM_ENV=$OPTARG
	    ;;
	s)
	    SUNBEAM_DIR=$OPTARG
	    ;;
	c)
	    PREFIX=$OPTARG
	    ;;
	u)
	    UPDATE=$OPTARG
	    ;;
	q)
	    OUTPUT=/dev/null
	    ;;
	h)
	    echo "Installs the Sunbeam metagenomics pipeline."
	    echo "  -e ENV_NAME      Conda environment name to create/update (default: sunbeam)"
	    echo "  -s SUNBEAM_DIR    Path containing Sunbeam rules (default: this directory)"
	    echo "  -c CONDA_DIR      Path to Anaconda/miniconda install (default: $HOME/miniconda3)"
	    echo "  -u [lib/env/all]  Update sunbeam[lib], conda [env], or both [all]"
	    echo "  -q                Suppress most Conda output"
	    echo "  -h                Display this message and exit"
	    exit 1
	    ;;
    esac
done

# Set defaults if not set by user
PREFIX=${PREFIX:-"${HOME}/miniconda3"}
SUNBEAM_ENV=${SUNBEAM_ENV:-sunbeam}
SUNBEAM_DIR=${SUNBEAM_DIR:-$SCRIPT_DIR}
UPDATE_ENV=false
UPDATE_LIB=false
case $UPDATE in
    lib)
	UPDATE_LIB=true	;;
    env)
	UPDATE_ENV=true	;;
    both)
	UPDATE_LIB=true
	UPDATE_ENV=true	;;
    *)
	;;
esac
OUTPUT=${OUTPUT:-/dev/stdout}

echo "Installation parameters:"
echo "  Conda installation:  ${PREFIX}"
echo "  Sunbeam environment: ${SUNBEAM_ENV}"
echo "  Sunbeam directory:   ${SUNBEAM_DIR}"
echo "  Update environment:  $(truefalse $UPDATE_ENV)"
echo "  Update sunbeamlib:   $(truefalse $UPDATE_LIB)"

export PATH=$PATH:$PREFIX/bin

install_conda () {
    wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $PREFIX >> $OUTPUT
    export PATH=$PATH:$PREFIX/bin
    command -v conda >/dev/null 2>&1 || {
	echo "Conda still isn't on the path, try installing manually"; exit 1;
    }
}

function install_env_vars {
    source activate $SUNBEAM_ENV
    echo -ne "#/bin/sh\nexport SUNBEAM_DIR=${SUNBEAM_DIR}" > \
	 ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh
    echo -ne "#/bin/sh\nunset SUNBEAM_DIR" > \
	 ${CONDA_PREFIX}/etc/conda/deactivate.d/env_vars.sh
}

function install_sunbeamlib {
    source activate $SUNBEAM_ENV
    pip install --upgrade --editable . >> $OUTPUT
    command -v sunbeam_init >/dev/null 2>&1 || {
	echo "Couldn't install Sunbeam; please report this as a bug."; exit 1;
    }
    echo "Sunbeam successfully installed.";
}

# Install conda if it doesn't show up on the path
command -v conda >/dev/null 2>&1 || {
    echo "Conda not installed, installing now"
    install_conda
    echo "Finished installing Conda."
}

$(conda env list | cut -f1 -d' ' | grep -Fxq $SUNBEAM_ENV >> $OUTPUT) || false && true
ENV_EXISTS=$?
ENV_CHANGED=false

if [ $UPDATE_ENV = true ]; then
    echo "Updating Sunbeam environment '${SUNBEAM_ENV}'"
    conda env update --name=$SUNBEAM_ENV --quiet -f environment.yml >> $OUTPUT
    ENV_CHANGED=true
elif [ $ENV_EXISTS -ne 0 ]; then
    echo "Creating new Sunbeam environment '${SUNBEAM_ENV}'"
    conda env create --name=$SUNBEAM_ENV --quiet -f environment.yml >> $OUTPUT
    ENV_CHANGED=true
else
    echo -ne "Skipping conda environment creation (${SUNBEAM_ENV} already exists). "
    echo "Re-run with '-u env' option to force upgrade."
fi

if [[ $UPDATE_ENV = true ]] || [[ $ENV_CHANGED ]]; then
    # Skip if we created a new environment, as it was just installed
    install_sunbeamlib >> $OUTPUT
fi

install_env_vars

echo "To get started, ensure ${PREFIX}/bin is in your PATH and run 'source activate $SUNBEAM_ENV'"

