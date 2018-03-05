#!/bin/bash
set -e

PREFIX=$HOME/miniconda3

while getopts "e:s:uqh" opt; do
    case $opt in
	e)
	    SUNBEAM_ENV=$OPTARG
	    ;;
	s)  SUNBEAM_DIR=$OPTARG
	    ;;
	u)
	    UPDATE=true
	    ;;
	q)
	    OUTPUT=/dev/null
	    ;;
	h)
	    echo -ne "Installs the Sunbeam metagenomics pipeline.\n"
	    echo -ne "  -e ENV_NAME     Conda environment name to create/update (default: sunbeam).\n"
	    echo -ne "  -s SUNBEAM_DIR  Path containing Sunbeam rules (default: this directory).\n"
	    echo -ne "  -u              Update an existing Sunbeam install.\n"
	    echo -ne "  -q              Suppress most Conda output.\n"
	    echo -ne "  -h              Display this message and exit.\n\n"
	    exit 1
	    ;;
    esac
done

# Set defaults if not set by user
SUNBEAM_ENV=${SUNBEAM_ENV:-sunbeam}
SUNBEAM_DIR=${SUNBEAM_DIR:-$BASH_SOURCE}
UPDATE=${UPDATE:-false}
OUTPUT=${OUTPUT:-/dev/tty}

echo "Installation parameters:"
echo "  Environment (-e):       ${SUNBEAM_ENV}"
echo "  Sunbeam directory (-s): ${SUNBEAM_DIR}"
echo "  Perform update (-u):    ${UPDATE}"

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
}

conda env list | cut -f1 -d' ' | grep -Fxq $SUNBEAM_ENV > /dev/null
ENV_EXISTS=$?

if [ $UPDATE = true ]; then
    echo "Updating Sunbeam environment '${SUNBEAM_ENV}'"
    conda env update --name=$SUNBEAM_ENV --quiet -f environment.yml >> $OUTPUT
    install_env_vars
    install_sunbeamlib
elif [ $ENV_EXISTS -ne 0 ]; then
    echo "Creating new Sunbeam environment '${SUNBEAM_ENV}'"
    conda env create --name=$SUNBEAM_ENV --quiet -f environment.yml >> $OUTPUT
    install_env_vars
    install_sunbeamlib
else
    install_env_vars
    echo "Skipping further installation (${SUNBEAM_ENV} already exists)".
    echo "Re-run with -u option to force upgrade."
fi

echo "To get started, ensure ${PREFIX}/bin is in your path and run 'source activate $SUNBEAM_ENV'"

