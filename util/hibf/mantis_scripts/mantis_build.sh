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
> $WORKING_DIRECTORY/bin_list.tmp
while IFS= read -r FASTA; do
    echo "$SQUEAKR_DIRECTORY/$(basename "${FASTA}" .fna).squeakr"
done < "$EXISTING_FILES" >> "$WORKING_DIRECTORY/bin_list.tmp"

ulimit -Sn $(ulimit -Hn)
$MANTIS_BINARY build \
    -t $THREADS \
    -s 34 \
    -i $WORKING_DIRECTORY/bin_list.tmp \
    -o $MANTIS_INDEX
rm $WORKING_DIRECTORY/bin_list.tmp
