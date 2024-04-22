#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -euo pipefail

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "The script '$(readlink -f $0)' must be sourced, not executed."
    echo "Please run 'source $(readlink -f $0)' instead."
    exit 1
fi

###
# Directory configuration
OUTPUT_DIR="/srv/data/seiler/part_simulation"

####
# Parameters supporting suffixes (e.g. 4G, 1M, 1024)
# The following suffixes are supported (case-sensitive): <none>, K, M, G, T, P
# The values must not contain spaces, i.e., "4 G" is not allowed.
# The suffixes are interpreted as binary suffixes, i.e., 1K = 1024.
REFERENCE_LENGTH=128M
NUMBER_OF_PARTS="64 512 1024"
NUMBER_OF_HAPLOTYPES=16
READ_COUNT=1M
# The total sequence size will be REFERENCE_LENGTH * NUMBER_OF_HAPLOTYPES

###
# Other parameters
READ_ERRORS=2
READ_LENGTH=250
SEED=42
THREADS=64

### Internal
expand_suffix() {
    local VAR=$1
    local VALUE=${!1}

    if [[ ${VALUE} =~ " " ]]; then
        echo "The given value must not contain spaces: '${VALUE}'"
        return 1
    fi

    eval ${VAR}=$(numfmt --from=iec "${VALUE}")
}

to_iec() {
    numfmt --to=iec $1 --format '%.0f'
}

expand_suffix "REFERENCE_LENGTH"
expand_suffix "READ_COUNT"
expand_suffix "NUMBER_OF_HAPLOTYPES"
