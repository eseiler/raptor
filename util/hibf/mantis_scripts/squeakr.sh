#!/usr/bin/env bash
# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

set -e
SCRIPT_ROOT=$(dirname $(readlink -f $0))
source $SCRIPT_ROOT/variables.sh

run_squeakr () {
    local FASTQ="$1"
    $SQUEAKR_BINARY count \
        -e \
        -k $K \
        -c 1 \
        -n \
        -t 1 \
        -o $SQUEAKR_DIRECTORY/$(basename ${FASTQ} .fastq).squeakr \
        ${FASTQ}
}

while IFS= read -r FASTA; do
        run_squeakr "$WORKING_DIRECTORY/bins/$(basename "${FASTA}" .fna).fastq"
done < <(cat "$INSERTION_FILES" "$EXISTING_FILES")