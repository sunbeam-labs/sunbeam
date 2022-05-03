#!/usr/bin/env bash

read -r -d '' __usage <<-'EOF'
  -e --environment  [arg] Sunbeam environment name. Default: "sunbeam"
  -s --sunbeam_dir  [arg] Location of Sunbeam source code. Default: this directory
  -c --conda  [arg]       Location of Conda installation. Default: ${PREFIX}
  -k --keep [arg]         Don't uninstall sunbeam [lib]rary, conda [env], or [all].
  -v --verbose            Show subcommand output
  -d --debug              Run in debug mode.
  -h --help               Display this message and exit.
EOF

# Load BASH3Boilerplate for command-line parsing and logging
source etc/b3bp.sh

function __err_report() {
    local error_code
    error_code=${?}
    error "Error in ${__file} in function ${1} on line ${2}"
    exit ${error_code}
}
trap '__err_report "${FUNCNAME:-.}" ${LINENO}' ERR  

# help mode
if [[ "${arg_h:?}" = "1" ]]; then
    # Help exists with code 1
    help "Help using ${0}:"
fi

# verbose mode
if [[ "${arg_v:?}" = "1" ]]; then
    LOG_LEVEL="7"
fi

# debug mode
if [[ "${arg_d:?}" = "1" ]]; then
    set -o xtrace
    LOG_LEVEL="7"
fi

function debug_capture () {
    debug "$(echo -e "$(${@})")"
}

function installation_error () {
    error "${1} failed!"
    if [[ "${arg_v:?}" != 1 && "${arg_d:?}" != 1 ]]; then
	error "Try re-running with -v or -d, or file an issue on Github."
    fi
    exit 1
}

# Set variables
__conda_path="${arg_c:-${HOME}/miniconda3}"
__sunbeam_dir="${arg_s:-$(readlink -f ${__dir})}"
__sunbeam_env="${arg_e:-sunbeam}"
__keep_lib=false
__keep_env=false
if [[ "${arg_k}" = "all" ]]; then
    __keep_lib=true
    __keep_env=true
elif [[ "${arg_k}" = "env" ]]; then
    __keep_env=true
elif [[ "${arg_k}" = "lib" ]]; then
    __keep_lib=true
fi

__old_path=$PATH
PATH=$PATH:${__conda_path}/bin

function __test_env() {
    if [[ $(__test_conda) = true ]]; then
	$(conda env list \
		 | cut -f1 -d' ' \
		 | grep -Fxq $__sunbeam_env > /dev/null) && \
	   echo true || echo false
    else
	echo false
    fi
}

function __test_sunbeam() {
    if [[ $(__test_env) = true ]]; then
	activate_sunbeam
	command -v sunbeam &> /dev/null && echo true || echo false
	deactivate_sunbeam
    else
	echo false
    fi
}

function enable_conda_activate () {
    # Allow conda [de]activate in this script
    CONDA_BASE=$(conda info --base) # see https://github.com/conda/conda/issues/7980
    source $CONDA_BASE/etc/profile.d/conda.sh
}

function deactivate_sunbeam () {
    enable_conda_activate
    set +o nounset
    conda deactivate
    set -o nounset
}

function uninstall_environment () {
    debug_capture conda env remove --name=$__sunbeam_env --quiet

    if [[ $(__test_env) == true ]]; then
	installation_error "Environment removal"
    fi
}

function uninstall_library () {
    debug_capture pip uninstall $__sunbeam_dir 2>&1
    if [[ $(__test_sunbeam) == true ]]; then
	installation_error "Library installation"
    fi
}

info "Starting Sunbeam removal..."
info "    Conda path:   ${__conda_path}"
info "    Sunbeam src:  ${__sunbeam_dir}"
info "    Sunbeam env:  '${__sunbeam_env}'"

# Remove Conda environment for Sunbeam
if [[ $__keep_env = false ]]; then
    info "Not removing conda environment"
else
    info "Removing Sunbeam environment..."
    uninstall_environment
fi

# Uninstalling Sunbeam lib
if [[ $__keep_lib = false ]]; then
    info "Not removing Sunbeam lib"
else
    info "Removing Sunbeam lib..."
    uninstall_library
fi

# Finalize removal
if [[ $__keep_lib = true ]]; then
    mv sunbeamlib/ ../
fi
cd ../
rm -r $__sunbeam_dir
    