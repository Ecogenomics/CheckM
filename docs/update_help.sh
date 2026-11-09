#!/usr/bin/env bash
#
# Regenerate the command help text included in the documentation.
#
# Writes the output of `checkm -h` to src/commands/help/checkm.txt and the
# output of `checkm <command> -h` to src/commands/help/<command>.txt for every
# subcommand defined in the checkm script.
#
# Usage:
#   export CHECKM_DATA_PATH=/path/to/checkm_data   # must be set (see below)
#   ./update_help.sh [path/to/checkm]              # default: checkm on PATH
#
# CHECKM_DATA_PATH must point to a CheckM reference data directory. Without
# it, CheckM reads the data path from its DATA_CONFIG file and, if no data is
# found there, starts an interactive prompt that rewrites DATA_CONFIG.

set -euo pipefail

CHECKM=${1:-$(command -v checkm || true)}
if [[ -z "${CHECKM}" || ! -x "${CHECKM}" ]]; then
    echo "error: checkm executable not found; pass its path as the first argument" >&2
    exit 1
fi

if [[ -z "${CHECKM_DATA_PATH:-}" ]]; then
    echo "error: set CHECKM_DATA_PATH to a CheckM reference data directory first" >&2
    exit 1
fi

HELP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/src/commands/help"
mkdir -p "${HELP_DIR}"

# argparse wraps help text to the terminal width; fix it so the output does not
# depend on the terminal used to run this script
export COLUMNS=100

# Remove machine-specific content: trailing whitespace, the path of the checkm
# executable and the CheckM data directory (if they ever appear), and any
# timestamped log lines.
clean() {
    sed -e 's/[[:space:]]*$//' \
        -e "s#${CHECKM_DATA_PATH%/}#<checkm_data_dir>#g" \
        -e "s#$(dirname "${CHECKM}")/##g" \
        -e '/^\[[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\} [0-9:]\{8\}\]/d' \
    | cat -s
}

# strip leading and trailing blank lines
trim() {
    sed -e '/./,$!d' | sed -e ':a' -e '/^\n*$/{$d;N;ba' -e '}'
}

"${CHECKM}" -h < /dev/null 2>&1 | clean | trim > "${HELP_DIR}/checkm.txt"
echo "wrote ${HELP_DIR}/checkm.txt"

# subcommands are taken from the argparse definitions in the checkm script
COMMANDS=$(grep -o "add_parser('[a-z_]*'" "${CHECKM}" | sed "s/add_parser('\(.*\)'/\1/")
if [[ -z "${COMMANDS}" ]]; then
    echo "error: no subcommands found in ${CHECKM}" >&2
    exit 1
fi

for cmd in ${COMMANDS}; do
    "${CHECKM}" "${cmd}" -h < /dev/null 2>&1 | clean | trim > "${HELP_DIR}/${cmd}.txt"
    echo "wrote ${HELP_DIR}/${cmd}.txt"
done
