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
> $QUERIES_FILE # Clear the content of the query file if it exists already.

#for FASTA in $(seq -f "$INPUT_DIR/$BIN_NUMBER/bins/bin_%0${#BIN_NUMBER}.0f.fasta" 0 1 $((BIN_NUMBER-1))); do
#    $HELPER_BINARY --input ${FASTA} \
#                   --output $WORKING_DIRECTORY/bins/$(basename ${FASTA} .fasta).fastq
#done
## CREAT THE INITIAL QUERY FILE
#  for FASTA in $(seq -f "$INPUT_DIR/$BIN_NUMBER/bins/bin_%0${#BIN_NUMBER}.0f.fasta" 0 1 $((CURRENT_BIN_NUMBER-1))); do
#    head -n 1 "$FASTA" >> $QUERIES_FILE   # Take the first line from the FASTA file
#    sampled_line=$(grep -v '^>' "$FASTA" | shuf -n 1)   # Sample a line from the FASTA file
#    echo "$sampled_line" >> "$QUERIES_FILE"    # Append the sampled line to the $QUERIES_FILE
#  done

while read -r FASTA; do
    if [[ $counter -ge 0 && $counter -le $((BIN_NUMBER-1))]]; then
         $HELPER_BINARY --input ${FASTA} \
                        --output $WORKING_DIRECTORY/bins/$(basename ${FASTA} .fasta).fastq

        if [[ $counter -ge 0 && $counter -le $((CURRENT_BIN_NUMBER-1))]]; then
          head -n 1 "$FASTA" >> $QUERIES_FILE   # Take the first line from the FASTA file
          sampled_line=$(grep -v '^>' "$FASTA" | shuf -n 1)   # Sample a line from the FASTA file TODO sample more lines, more queries.
          echo "$sampled_line" >> "$QUERIES_FILE"    # Append the sampled line to the $QUERIES_FILE
        fi
    fi

    ((counter++))     # Increment the counter
done < "$ALL_FILENAMES"

awk 'NR%4==2' $QUERIES_FILE > $READ_FILE #this select every 2nd line. Although, from the github page, it seems as if fasta files are also allowed as input, so I don't know if this is still nesessary?
