#!/usr/bin/env bash
# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

set -e

READ_LENGTH=100
ERRORS=2
K=20
THREADS=32
BIN_NUMBER=1024
CURRENT_BIN_NUMBER=10 # The number of bins that should be used for building the initial index.
BATCH_SIZE=10
# /group/ag_abi/myrthew00/dynamic_mantis/util/Genome_Biology/mantis_scripts

SQUEAKR_BINARY="/scratch/dynamic_mantis/bin/squeakr"
MANTIS_BINARY="/scratch/dynamic_mantis/bin/mantis"
HELPER_BINARY="/scratch/dynamic_mantis/bin/fasta_to_fastq"
INPUT_DIR="/group/ag_abi/myrthew00/raptor/build/bin/example_data" #"/scratch/dynamic_mantis/big_dataset"
BENCHMARK_DIR="/group/ag_abi/myrthew00/dynamic_mantis/mantis_bench"

WORKING_DIRECTORY=$BENCHMARK_DIR/$BIN_NUMBER
mkdir -p $WORKING_DIRECTORY/bins/
mkdir -p $WORKING_DIRECTORY/reads/

SQUEAKR_DIRECTORY=$WORKING_DIRECTORY/squeakr
mkdir -p $SQUEAKR_DIRECTORY

MANTIS_INDEX=$WORKING_DIRECTORY/mantis/
mkdir -p $MANTIS_INDEX

QUERIES_FILE=/group/ag_abi/myrthew00/raptor/build/bin/evaluation/tmp/queries.fasta
READ_FILE=$WORKING_DIRECTORY/reads/all_10.fastq
QUERY_OUT=$WORKING_DIRECTORY/result.txt
PYTHON_OUTPUT=$WORKING_DIRECTORY/mantis_result.txt

copy_time=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_copy.time
copy_log=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_copy.log
squeakr_time=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_squeakr.time
squeakr_log=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_squeakr.log
mantis_build_time=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_mantis_build.time
mantis_build_log=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_mantis_build.log
mantis_mst_time=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_mantis_mst.time
mantis_mst_log=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_mantis_mst.log
mantis_query_time=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_mantis_query.time
mantis_query_log=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_mantis_query.log
mantis_update_log=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_mantis_update.log
mantis_update_time=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_mantis_update.time

