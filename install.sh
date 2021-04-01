#!/usr/bin/env bash

__conda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

read -r -d '' __usage <<-'EOF'
  -e --environment  [arg] Environment to install to. Default: "sunbeam"
  -s --sunbeam_dir  [arg] Location of Sunbeam source code. Default: this directory
  -c --conda  [arg]       Location of Conda installation. Default: ${PREFIX}
  -u --update [arg]       Update sunbeam [lib]rary, conda [env], or [all].
  -m --mamba              Install and use mamba in base environment as alternative dependency solver
  -v --verbose            Show subcommand output
  -d --debug              Run in debug mode.
  -h --help               Display this message and exit.
EOF

read -r -d '' __helptext <<-'EOF'
 This script installs or upgrades Sunbeam, including Conda (if not installed).
 To upgrade, pass the '--upgrade all' option, then be sure to update your config
 files using 'sunbeam config update'.
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
__update_lib=false
__update_env=false
__install_mamba=false
if [[ "${arg_u}" = "all" || "${arg_u}" = "env" ]]; then
    __update_lib=true
    __update_env=true
elif [[ "${arg_u}" = "lib" ]]; then
    __update_lib=true
fi

if [[ "${arg_m:?}" = "1" ]]; then
    __install_mamba=true
fi

__old_path=$PATH
PATH=$PATH:${__conda_path}/bin

function __git_dir_exists() {
    if [ -d ".git" ]; then
      echo true
    else
      echo false
    fi
}

function __test_conda() {
    command -v conda &> /dev/null && echo true || echo false
}

function __test_mamba() {
    command -v mamba &> /dev/null && echo true || echo false
}

function __detect_conda_install() {
    local discovered=$(__test_conda)
    if [[ $discovered = true ]]; then
	local conda_path="$(which conda)"
	echo ${conda_path%'/bin/conda'}
    fi
}    

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

function activate_sunbeam () {
    enable_conda_activate
    set +o nounset
    conda activate $__sunbeam_env
    set -o nounset
}

function deactivate_sunbeam () {
    enable_conda_activate
    set +o nounset
    conda deactivate
    set -o nounset
}

function install_conda () {
    local tmpdir=$(mktemp -d)
    debug "Downloading miniconda..."
    debug_capture wget -nv ${__conda_url} -O ${tmpdir}/miniconda.sh 2>&1
    debug "Installing miniconda..."
    debug_capture bash ${tmpdir}/miniconda.sh -b -p ${__conda_path} 2>&1
    if [[ $(__test_conda) != true ]]; then
	installation_error "Environment creation"
    fi
    rm ${tmpdir}/miniconda.sh
}

function install_environment () {
    if [[ $(__test_mamba) = true ]]; then
        cmd=mamba
    else
        cmd=conda
    fi
    debug_capture $cmd env update --name=$__sunbeam_env \
			  --quiet --file environment.yml
    if [[ $(__test_env) != true ]]; then
	installation_error "Environment creation"
    fi
}

function install_env_vars () {
    activate_sunbeam
    echo -ne "#/bin/sh\nexport SUNBEAM_DIR=${__sunbeam_dir}" > \
	 ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh
    echo -ne "#/bin/sh\nunset SUNBEAM_DIR" > \
	 ${CONDA_PREFIX}/etc/conda/deactivate.d/env_vars.sh
}

function install_sunbeamlib () {
    activate_sunbeam
    if [[ $(__git_dir_exists) != true ]]; then
      installation_error "Sunbeam requires a git clone \
(instead of a compressed archive) to detect its version"
    fi
    debug_capture pip install --upgrade $__sunbeam_dir 2>&1
    if [[ $(__test_sunbeam) != true ]]; then
	installation_error "Library installation"
    fi
}

info "Starting Sunbeam installation..."
info "    Conda path:   ${__conda_path}"
info "    Sunbeam src:  ${__sunbeam_dir}"
info "    Sunbeam env:  '${__sunbeam_env}'"

debug "Components detected:"
__conda_installed=$(__test_conda)
debug "    Conda:       ${__conda_installed}"
__mamba_installed=$(__test_mamba)
debug "    Mamba:       ${__mamba_installed}"
__env_exists=$(__test_env)
debug "    Environment: ${__env_exists}"
__sunbeam_installed=$(__test_sunbeam)
debug "    Library:     ${__sunbeam_installed}"

__env_changed=false

# Install Conda if necessary
if [[ $__conda_installed = true ]]; then
    if [[ $(__detect_conda_install) != $__conda_path ]]; then
	warning "Found pre-existing Conda installation in $(__detect_conda_install)".
	warning "Ignoring specified Conda path in favor of existing Conda install."
	__conda_path=$(__detect_conda_install)
    fi
    info "Conda already installed."
else
    info "Installing Conda..."
    install_conda
    __env_changed=true
fi

# Install mamba if necessary
if [[ $__install_mamba = true ]]; then
    info "Installing mamba..."
    conda install --yes --quiet -n base -c conda-forge mamba
fi

# Create Conda environment for Sunbeam
if [[ $__env_exists = true && $__update_env = false ]]; then
    info "Specified environment already exists (use '--update env' to update)"
else
    info "Creating Sunbeam environment..."
    install_environment
    __env_changed=true
fi

# Check again to ensure success

# Install sunbeamlib into environment if changed or requested
if [[ $__env_changed = true ]]; then
    info "Environment installed/updated; (re)installing Sunbeam library..."
    install_sunbeamlib
elif [[ $__sunbeam_installed = false ]]; then
    info "Installing Sunbeam library..."
    install_sunbeamlib
elif [[ $__update_lib = true ]]; then
    info "Updating Sunbeam library..."
    install_sunbeamlib
else
    info "Sunbeam library already installed (use '--update lib' to update)"
fi

# Always update the env_vars.sh in the sunbeam environment
debug "Updating \$SUNBEAM_DIR variable to point to ${__sunbeam_dir}"
install_env_vars

# Check if on pre-existing path
if [[ $__old_path != *"${__conda_path}/bin"* ]]; then
    warning "** Conda was not detected on your PATH. **"
    warning "This is normal if you haven't installed Conda before."
    warning "To add it to your path, run "
    warning "   'echo \"export PATH=\$PATH:${__conda_path}/bin\" >> ~/.bashrc'"
    warning "and close and re-open your terminal session to apply."
    warning "When finished, run 'conda activate ${__sunbeam_env}' to begin."
else
    info "Done. Run 'conda activate ${__sunbeam_env}' to begin."
fi

   


