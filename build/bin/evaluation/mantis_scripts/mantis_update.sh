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
CURRENT_BIN_NUMBER=11 # adapt this variable in the root.

seq -f "$SQUEAKR_DIRECTORY/bin_%0${#BIN_NUMBER}.0f.squeakr" $CURRENT_BIN_NUMBER 1 $(($CURRENT_BIN_NUMBER + $BATCH_SIZE -1)) > $WORKING_DIRECTORY/bin_list.tmp

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
