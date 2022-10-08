#!/usr/bin/env bash

read -r -d '' __usage <<-'EOF'
  -l --list               List all installed versions of sunbeam.
  -s --switch [arg]       Switch to a new version of sunbeam (install if not installed).
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

