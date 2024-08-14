#!/usr/bin/env bash
# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

#https://raw.githubusercontent.com/seqan/raptor/main/util/Genome_Biology/refseq_truth/refseq.sample.ids => 25320	files.

set -e
SCRIPT_ROOT=$(dirname $(readlink -f $0))
source $SCRIPT_ROOT/variables.sh

########### FUNCTIONS
function convertElapsedTime() {     # Splitting the time string into hours, minutes, and seconds
    time_str="$1"
    IFS=':' read -r -a time_parts <<< "$time_str"

    if [ "${#time_parts[@]}" -eq 3 ]; then
        hours="${time_parts[0]}"
        minutes="${time_parts[1]}"
        seconds="${time_parts[2]}"

    elif [ "${#time_parts[@]}" -eq 2 ]; then
        hours=0
        minutes="${time_parts[0]}"
        seconds="${time_parts[1]}"
    else
        echo "Invalid time format: $time_str"
        return 1
    fi
    # Converting hours, minutes, and seconds to seconds
    total_seconds=$(echo "($hours * 3600) + ($minutes * 60) + $seconds" | bc)
    echo "$total_seconds"
}

join_by_comma() { # Function to join array elements with commas
  local IFS=","
  echo "$*"
}



INSERTION_TIMINGS=()  # Declare an empty array to store the timings
INSERTION_MEMORY=()  # Declare an empty array to store the memory values
QUERY_TIMINGS=()
QUERY_MEMORY=()
INDEX_SIZES=()
ERRORS=()
number_of_existing_bins=$(wc -l < "$EXISTING_FILES")
CURRENT_BIN_NUMBER=1
NUMBER_OF_FILES=$(( $(wc -l < "$INSERTION_FILES") + CURRENT_BIN_NUMBER ))
save_results() {
    ERRORS+=($((CURRENT_BIN_NUMBER + number_of_existing_bins)))     # Note down that the error occurs
    echo "Saving results..."
    echo "time_insertion = [$(join_by_comma "${INSERTION_TIMINGS[@]}")]" > $PYTHON_OUTPUT
    echo "time_query = [$(join_by_comma "${QUERY_TIMINGS[@]}")]" >> $PYTHON_OUTPUT
    echo "memory_insertion = [$(join_by_comma "${INSERTION_MEMORY[@]}")]" >> $PYTHON_OUTPUT
    echo "memory_query = [$(join_by_comma "${QUERY_MEMORY[@]}")]" >> $PYTHON_OUTPUT
    echo "size_index = [$(join_by_comma "${INDEX_SIZES[@]}")]" >> $PYTHON_OUTPUT
    echo "rebuild = []" >> $PYTHON_OUTPUT
    echo "errors = [] [$(join_by_comma "${ERRORS[@]}")]" >> $PYTHON_OUTPUT
}
##################


trap 'save_results; echo; echo "## [$(date +"%Y-%m-%d %T")] ERROR ##"' ERR

echo "## [$(date +"%Y-%m-%d %T")] Start ##"

echo -n "[$(date +"%Y-%m-%d %T")] Copying input..."
    /usr/bin/time -o $copy_time -v \
        $SCRIPT_ROOT/copy_input.sh \
            &>$copy_log
echo "Done."

echo -n "[$(date +"%Y-%m-%d %T")] Running squeakr..."
    /usr/bin/time -o $squeakr_time -v \
        $SCRIPT_ROOT/squeakr.sh \
            &>$squeakr_log
echo "Done."

echo -n "[$(date +"%Y-%m-%d %T")] Running mantis build..."
    /usr/bin/time -o $mantis_build_time -v \
        $SCRIPT_ROOT/mantis_build.sh \
            &>$mantis_build_log
echo "Done."

echo -n "[$(date +"%Y-%m-%d %T")] Running mantis mst..."
    /usr/bin/time -o $mantis_mst_time -v \
        $SCRIPT_ROOT/mantis_mst.sh \
            &>$mantis_mst_log
echo "Done."



while ((CURRENT_BIN_NUMBER <= NUMBER_OF_FILES)); do #looping over CURRENT_BIN_NUMBER to BIN_NUMBER, with stepsize = BATCH_SIZE. # inside loop: update, query, measure time. Update the query file
    echo "Current bin number: "
    echo $((CURRENT_BIN_NUMBER + number_of_existing_bins))
    #write the new batch of files from CURRENT_BIN_NUMBER to CURRENT_BIN_NUMBER + BATCH_SIZE to a query file temporary file .
    end_line=$((CURRENT_BIN_NUMBER + BATCH_SIZE - 1))
    sed -n "${CURRENT_BIN_NUMBER},${end_line}p" "$INSERTION_FILES" > "$WORKING_DIRECTORY/bin_list.tmp"
    sed -n "${CURRENT_BIN_NUMBER},${end_line}p" "$INSERTION_FILES" >> "$WORKING_DIRECTORY/tmp_queries.txt" #append the list of existing files
    if [[ ! -s "$WORKING_DIRECTORY/bin_list.tmp" ]]; then      # File is empty or contains only an empty line
        echo "File is empty or contains only an empty line"
        break
    fi

    folder_size=$(du -s "$MANTIS_INDEX" | cut -f1)
    INDEX_SIZES+=($folder_size)

    # QUERY
    echo -n "[$(date +"%Y-%m-%d %T")] Running mantis query..."
    SCRIPT_ROOT=$(dirname $(readlink -f $0))
        /usr/bin/time -o $mantis_query_time -v \
            $SCRIPT_ROOT/mantis_query.sh \
                &>$mantis_query_log || true
    echo "Done."
    elapsed_time=$(grep Elapsed ${mantis_query_time} | cut -d ' ' -f 8)
    elapsed_seconds=$(convertElapsedTime "$elapsed_time")
    QUERY_TIMINGS+=("elapsed_seconds")
    memory_str=$(grep Max ${mantis_query_time} | cut -d ' ' -f 6)
    QUERY_MEMORY+=("$memory_str")

    # UPDATE
    echo -n "[$(date +"%Y-%m-%d %T")] Running mantis update..."
        /usr/bin/time -o $mantis_update_time -v \
            $SCRIPT_ROOT/mantis_update.sh \
                &>$mantis_update_log || true
    echo "Done." # TODO check if the new user bins are appended to the query file.
    elapsed_time=$(grep Elapsed ${mantis_update_time} | cut -d ' ' -f 8)     # extract time and memory from  mantis_update_log, and append those to the python files.
    elapsed_seconds=$(convertElapsedTime "$elapsed_time")
    INSERTION_TIMINGS+=("$elapsed_seconds")
    memory_str=$(grep Max ${mantis_update_time} | cut -d ' ' -f 6) #| numfmt --to iec --format "%8.2f" --from-unit Ki --round nearest)
    INSERTION_MEMORY+=("$memory_str")  # Append the converted memory value to the array

    # Increment the current bin number by the batch size
    CURRENT_BIN_NUMBER=$((CURRENT_BIN_NUMBER + BATCH_SIZE))
    if (( CURRENT_BIN_NUMBER % 50 == 0 )); then
      save_results
    fi
done

save_results


echo "## [$(date +"%Y-%m-%d %T")] End ##"
