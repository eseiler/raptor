#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -euo pipefail

SCRIPT_ROOT=$(dirname $(readlink -f $0))
source ${SCRIPT_ROOT}/config.sh

trap_script() {
    exec 2>&4 1>&3

    if [ $status -ne 0 ]; then
        echo "[ERROR] The log file can be found at ${__LOG_FILE}"
        echo "[ERROR] Last 10 lines of the log file follow"
        tail -n 10 ${__LOG_FILE} | xargs -d '\n' -L 1 echo "    "
    else
        echo "[Success] The log file can be found at ${__LOG_FILE}"
        echo "[Success] The output can be found in ${OUTPUT_DIR}"
    fi

    if [[ -v __SIMULATION_TMP_DIR && -d ${__SIMULATION_TMP_DIR} ]]; then
        rm -fdr ${__SIMULATION_TMP_DIR}
    fi
}
trap 'status=$?; set +x; trap_script; exit $status' EXIT INT
trap 'status=$?; set +x; exit $status' ABRT HUP INT PIPE QUIT TERM

quiet_loop() {
    if [[ $- =~ x ]]; then
        set +x
    fi
}

__LOG_DIR=${OUTPUT_DIR}/logs
__BUILD_DIR=${OUTPUT_DIR}/build
__BINARY_DIR=${__BUILD_DIR}/bin

mkdir -p ${__LOG_DIR}
mkdir -p ${__BUILD_DIR}
mkdir -p ${__BINARY_DIR}

__LOG_FILE=${__LOG_DIR}/$(date +"%Y-%m-%d_%H-%M-%S").log
touch ${__LOG_FILE}
echo "## Log file can be found at ${__LOG_FILE}"
exec 3>&1 4>&2 1>>${__LOG_FILE} 2>&1

set -x

echo "## [$(date +"%Y-%m-%d %T")] Building Dependencies" | tee /dev/fd/3
cd ${__BUILD_DIR}
if [[ ! -d ${__BUILD_DIR}/raptor ]]; then
    git clone https://github.com/seqan/raptor.git
fi

cmake ${__BUILD_DIR}/raptor/util/iScience -DCMAKE_BUILD_TYPE=Release \
                                          -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
                                          -DRAPTOR_UTILITY_BUILD_MASON=ON \
                                          -Wno-dev \
                                          -DINSTALL_RAPTOR=OFF

make -j${THREADS} --no-print-directory install

set +x

for __NUMBER_OF_PARTS in ${NUMBER_OF_PARTS}; do
    expand_suffix "__NUMBER_OF_PARTS"

    __SIMULATION_DIR=${OUTPUT_DIR}/${__NUMBER_OF_PARTS}
    __SIMULATION_BIN_DIR=${__SIMULATION_DIR}/bins
    __SIMULATION_INFO_DIR=${__SIMULATION_DIR}/info
    __SIMULATION_TMP_DIR=${__SIMULATION_DIR}/TMP_$(date +"%Y-%m-%d_%H-%M-%S")

    mkdir -p ${__SIMULATION_DIR}
    mkdir -p ${__SIMULATION_BIN_DIR}
    mkdir -p ${__SIMULATION_INFO_DIR}
    mkdir -p ${__SIMULATION_TMP_DIR}

    set -x
    echo "## [$(date +"%Y-%m-%d %T")] Number of parts: ${__NUMBER_OF_PARTS}" | tee /dev/fd/3

    echo "## [$(date +"%Y-%m-%d %T")] Simulating ${__NUMBER_OF_PARTS} parts with a total length of $(to_iec ${REFERENCE_LENGTH})" | tee /dev/fd/3
    # Simulate reference
    ${__BINARY_DIR}/mason_genome \
        --contig-length ${REFERENCE_LENGTH} \
        --out-file ${__SIMULATION_TMP_DIR}/ref.fasta \
        --seed ${SEED} \
        --quiet

    # Evenly distribute it over bins
    ${__BINARY_DIR}/split_sequence \
        --input ${__SIMULATION_TMP_DIR}/ref.fasta \
        --parts ${__NUMBER_OF_PARTS}

    # We do not need the reference anymore
    rm ${__SIMULATION_TMP_DIR}/ref.fasta

    echo "## [$(date +"%Y-%m-%d %T")] Simulating $(to_iec ${NUMBER_OF_HAPLOTYPES}) haplotypes for each part" | tee /dev/fd/3
    # Simulate haplotypes for each bin
    for file in ${__SIMULATION_TMP_DIR}/*.fa
    do
        (
            ${__BINARY_DIR}/mason_variator \
                --in-reference ${file} \
                --num-haplotypes ${NUMBER_OF_HAPLOTYPES} \
                --out-fasta ${__SIMULATION_TMP_DIR}/$(basename ${file} .fa).fasta \
                --out-vcf ${__SIMULATION_INFO_DIR}/$(basename ${file} .fa).vcf \
                --quiet &>/dev/null
            gzip --force ${__SIMULATION_INFO_DIR}/$(basename ${file} .fa).vcf
        ) &

        if [[ $(jobs -r -p | wc -l) -ge ${THREADS} ]]; then
            wait -n
        fi
        quiet_loop
    done
    wait
    set -x

    echo "## [$(date +"%Y-%m-%d %T")] Simulating $(to_iec ${READ_COUNT}) reads of length ${READ_LENGTH} with ${READ_ERRORS} errors" | tee /dev/fd/3
    __LIST_FILE=${__SIMULATION_TMP_DIR}/all_bins.txt
    find ${__SIMULATION_TMP_DIR} -type f -name "*.fasta" | sort -V > ${__LIST_FILE}
    ${__BINARY_DIR}/generate_reads \
        --output ${__SIMULATION_TMP_DIR} \
        --errors ${READ_ERRORS} \
        --number_of_reads ${READ_COUNT} \
        --read_length ${READ_LENGTH} \
        --number_of_haplotypes ${NUMBER_OF_HAPLOTYPES} \
        ${__LIST_FILE}
    find ${__SIMULATION_TMP_DIR} -type f -name "*.fastq" | sort -V | xargs cat > ${__SIMULATION_DIR}/reads_e${READ_ERRORS}_${READ_LENGTH}.fastq

    echo "## [$(date +"%Y-%m-%d %T")] Splitting $(to_iec ${NUMBER_OF_HAPLOTYPES}) haplotypes into individual files" | tee /dev/fd/3
    cd ${__SIMULATION_BIN_DIR}
    for file in ${__SIMULATION_TMP_DIR}/*.fasta
    do
        (
            csplit \
                --prefix=$(basename ${file} .fasta)_ \
                --digits=${#__NUMBER_OF_PARTS} \
                --suffix-format="%0${#NUMBER_OF_HAPLOTYPES}d.fasta" \
                --elide-empty-files \
                --silent \
                ${file} \
                "/^>/" {*}
        ) &

        if [[ $(jobs -r -p | wc -l) -ge ${THREADS} ]]; then
            wait -n
        fi
        quiet_loop
    done
    wait
    set -x

    rm -fdr ${__SIMULATION_TMP_DIR}

    find ${__SIMULATION_BIN_DIR} -type f -name "*.fasta" | sort -V > ${__SIMULATION_DIR}/filenames.txt
    cp --no-preserve=mode ${SCRIPT_ROOT}/config.sh ${__SIMULATION_INFO_DIR}
    cp --no-preserve=mode $(readlink -f $0) ${__SIMULATION_INFO_DIR}

    set +x
done
