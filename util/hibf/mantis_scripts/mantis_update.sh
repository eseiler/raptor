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

# convert the files in the bin list to their squeakr filenames.
while IFS= read -r FASTA; do
      echo $SQUEAKR_DIRECTORY/$(basename ${FASTA} .fna).squeakr > $WORKING_DIRECTORY/bin_list_copy.tmp
done <<< "$(cat "$WORKING_DIRECTORY/bin_list.tmp")"
mv $WORKING_DIRECTORY/bin_list_copy.tmp $WORKING_DIRECTORY/bin_list.tmp

ulimit -Sn $(ulimit -Hn)
$MANTIS_BINARY build \
    -t $THREADS \
    -s 34 \
    -i $WORKING_DIRECTORY/bin_list.tmp \
    -o $WORKING_DIRECTORY/mantis_tmp # mantis gives an error if bin_list.tmp is empty  2834345 Gleitkomma-Ausnahme
$MANTIS_BINARY mst \
    -p $WORKING_DIRECTORY/mantis_tmp \
    -t $THREADS \
    -d

$SCRIPT_ROOT/mantis_merge.sh

mv $WORKING_DIRECTORY/mantis_merged $MANTIS_INDEX # move the folder of the new mantis to the orignal one.
rm $WORKING_DIRECTORY/bin_list.tmp


#[2023-06-30 13:28:51.029] [mantis_console] [info] Running kruskal
#ERROR in selfParent => idx > vector size: 4 4; later in the code.