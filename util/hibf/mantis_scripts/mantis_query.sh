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

# Use bin_list.tmp to create query file
$READS_BINARY --output "$WORKING_DIRECTORY/tmp_queries.fastq" \
                 --errors "$ERRORS" \
                 --read_length "$READ_LENGTH" \
                 --number_of_reads "$NUMBER_OF_READS" \
                 "$WORKING_DIRECTORY/tmp_queries.txt"
> $READ_FILE
awk 'NR%4==2' "$WORKING_DIRECTORY/tmp_queries.fastq" > $READ_FILE #this select every 2nd line. Although, from the github page, it seems as if fasta files are also allowed as input, so I don't know if this is still nesessary?

$MANTIS_BINARY query -b \
    -p $MANTIS_INDEX \
    -k $K \
    -o $QUERY_OUT \
    $READ_FILE

rm $WORKING_DIRECTORY/tmp_queries.fastq
