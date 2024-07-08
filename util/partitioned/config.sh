#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -euo pipefail

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "The script '$(readlink -f $0)' must be sourced, not executed."
    echo "There is potentially another script that sources this script."
    echo "If not, please run 'source $(readlink -f $0)' instead."
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
NUMBER_OF_PARTS="64 512 1024 655360"
NUMBER_OF_HAPLOTYPES=16
READ_COUNT=10M
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

check_config() {
    if [[ ${REFERENCE_LENGTH} -lt ${READ_LENGTH} ]]; then
        echo "The reference length (${REFERENCE_LENGTH}) must be greater than the read length (${READ_LENGTH})."
        return 1
    fi

    if [[ ${READ_ERRORS} -gt ${READ_LENGTH} ]]; then
        echo "The number of read errors (${READ_ERRORS}) must be less than the read length (${READ_LENGTH})."
        return 1
    fi

    for PART in ${NUMBER_OF_PARTS}; do
        if [[ ${PART} -gt ${READ_COUNT} ]]; then
            echo "The number of parts (${PART}) must be less than the number of reads (${READ_COUNT})."
            return 1
        fi

        if [[ $((${READ_COUNT} % ${PART})) -ne 0 ]]; then
            echo "The number of reads (${READ_COUNT}) must be divisible by the number of parts (${PART})."
            return 1
        fi
    done
}

expand_suffix "REFERENCE_LENGTH"
expand_suffix "READ_COUNT"
expand_suffix "NUMBER_OF_HAPLOTYPES"

check_config
