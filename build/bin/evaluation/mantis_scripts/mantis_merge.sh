set -e
SCRIPT_ROOT=$(dirname $(readlink -f $0))
source $SCRIPT_ROOT/variables.sh
#mkdir $WORKING_DIRECTORY/mantis_merged

#seq -f "$SQUEAKR_DIRECTORY/bin_%0${#BIN_NUMBER}.0f.squeakr" 0 1 $((BIN_NUMBER-1)) > $WORKING_DIRECTORY/bin_list.tmp
#ulimit -Sn $(ulimit -Hn) #i2 will be the smaller index containing the batch.
$MANTIS_BINARY merge -rm \
	    -i1 $MANTIS_INDEX \
	    -i2 $WORKING_DIRECTORY/mantis_tmp \
	    -t $THREADS \
	    -o $WORKING_DIRECTORY/mantis_merged # $MANTIS_INDEX



    #mantis merge -i1 <input_dir_1> -i2 <input_dir_2> [-t <num_threads>] -o <merge_output>
