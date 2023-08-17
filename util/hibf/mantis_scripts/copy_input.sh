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
> "$WORKING_DIRECTORY/tmp_queries.txt" # Clear the content of the query file if it exists already.
> $READ_FILE

while IFS= read -r FASTA || [[ -n "$FASTA" ]]; do
         $HELPER_BINARY --input "${FASTA}" \
                        --output $WORKING_DIRECTORY/bins/$(basename "${FASTA}" .fna).fastq
done < <(cat "$INSERTION_FILES" "$EXISTING_FILES")

cp "$EXISTING_FILES" "$WORKING_DIRECTORY/tmp_queries.txt"

#while IFS= read -r FASTA || [[ -n "$FASTA" ]]; do
#  for ((i=1; i<=$NUMBER_OF_QUERIES; i++)); do
#          head -n 1 "$FASTA" >> $QUERIES_FILE   # Take the first line from the FASTA file
#          sampled_line=$(grep -v '^>' "$FASTA" | shuf -n 1)   # Sample a line from the FASTA file TODO sample more lines, more queries.
#          echo "$sampled_line" >> "$QUERIES_FILE"    # Append the sampled line to the $QUERIES_FILE
#  done
#done <<< "$(cat "$EXISTING_FILES")"

awk 'NR%4==2' "$WORKING_DIRECTORY/tmp_queries.txt" > $READ_FILE #this select every 2nd line. Although, from the github page, it seems as if fasta files are also allowed as input, so I don't know if this is still nesessary?
