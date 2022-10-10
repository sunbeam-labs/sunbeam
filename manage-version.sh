#!/usr/bin/env bash

read -r -d '' __usage <<-'EOF'
  -l --list   [arg]       List all [installed] or all [available] versions of sunbeam.
  -c --current            List environment for the code currently installed (current branch tag).
  -s --switch [arg]       Switch to a new version of sunbeam (install if not installed).
  -r --remove [arg]       Uninstall the specified version of sunbeam.
  -v --verbose            Show subcommand output.
  -d --debug              Run in debug mode.
  -h --help               Display this message and exit.
EOF

read -r -d '' __helptext <<-'EOF'
 This script manages your sunbeam version, switching between versions, installing new 
 ones, or listing what's currently installed.
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

function __test_conda() {
    command -v conda &> /dev/null && echo true || echo false
}

function __test_env() {
    if [[ $(__test_conda) = true ]]; then
	$(conda env list \
		 | cut -f1 -d' ' \
		 | grep -Fxq $1 > /dev/null) && \
	   echo true || echo false
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
    conda activate $1
    set -o nounset
}

function deactivate_sunbeam () {
    enable_conda_activate
    set +o nounset
    conda deactivate
    set -o nounset
}

debug_capture git pull

# list current
if [[ "${arg_c:?}" = "1" ]]; then
    TAG=$(git describe --tag)
    ENV_NAME="sunbeam${TAG:1}"
    info "Environment for installed code is ${ENV_NAME}."
    if [[ $(__test_env $ENV_NAME) = true ]]; then
        info "Activate with 'conda activate ${ENV_NAME}'."
    else
        info "Environment doesn't exist, install by running './install.sh'."
    fi
    exit 0
fi

if [[ "${arg_l}" = "installed" ]]; then
    info "Installed sunbeam envs:"

    conda env list |
    while read line
    do
        if [[ "${line:0:7}" = "sunbeam" ]]; then
            array_form=(${line})
            echo "${array_form[0]}"
        fi
    done
elif [[ "${arg_l}" = "available" ]]; then
    git tag --list
    git branch | while read line; do echo ${line}; done
fi

if [[ ! -z "${arg_s}" ]]; then
    # Switch to new branch
    if [[ "${arg_s}" = "dev" ]]; then
        info "Switching to branch dev ..."
        git checkout dev
    elif [[ "${arg_s}" = "stable" ]]; then
        info "Switching to branch stable ..."
        git checkout stable
    else
        __is_branch=false
        __is_tag=false
        __cleaned_name="${arg_s}"
        if [[ "${arg_s:0:7}" = "sunbeam" ]]; then
            __cleaned_name="${arg_s:7}"
        fi
        if [[ "${__cleaned_name:0:1}" = "v" ]]; then
            __cleaned_name="${__cleaned_name:1}"
        fi

        git branch |
        while read line
        do
            if [[ "${__cleaned_name}" = "${line}" ]]; then
                info "Switching to branch ${__cleaned_name} ..."
                git checkout $__cleaned_name
                __is_branch=true
                break
            fi
        done

        git tag --list |
        while read line
        do
            if [[ "$__cleaned_name" = "${line}" ]]; then
                info "Switching to release v${__cleaned_name} ..."
                git checkout tags/v${__cleaned_name} -b ${__cleaned_name}
                __is_tag=true
                break
            fi
        done

        # Check if no cases were hit
        if [[ $__is_branch = false ]]; then
            if [[ $__is_tag = false ]]; then
                error "Couldn't find ${arg_s}. Make sure to use a valid branch name or version tag."
                exit 1
            fi
        fi
    fi

    # Switch conda env
    __env_tag=$(git describe --tag)
    __env_name="sunbeam${__env_tag:1}"
    deactivate_sunbeam

    if [[ $(__test_env ${__env_name}) == true ]]; then
        info "Found existing environment, activate with 'conda activate ${__env_name}'"
    else
        info "Couldn't find environment for ${__env_name}, installing..."
        ./install.sh -e ${__env_name}
    fi

    #activate_sunbeam $__env_name # Only works within script, won't change user's context
fi

if [[ ! -z "${arg_r}" ]]; then
    __env_name="${arg_r}"
    if [[ "${arg_r:0:7}" != "sunbeam" ]]; then
        __env_name="sunbeam${arg_r}"
    fi

    if [[ "$CONDA_DEFAULT_ENV" = "${__env_name}" ]]; then
        enable_conda_activate
        conda deactivate
    fi

    if [[ $(__test_env ${__env_name}) == true ]]; then
        info "Removing ${__env_name}"
        conda env remove -n ${__env_name}
        if [[ $(__test_env ${__env_name}) != true ]]; then
            info "Removed ${__env_name}"
        else
            error "Failed to remove ${__env_name}"
        fi
    else
        error "Can't find env ${__env_name}"
    fi
fi

   


