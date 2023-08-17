set -e
SCRIPT_ROOT=$(dirname $(readlink -f $0))
source $SCRIPT_ROOT/variables.sh

$MANTIS_BINARY merge -rm \
	    -i1 $MANTIS_INDEX \
	    -i2 $WORKING_DIRECTORY/mantis_tmp \
	    -t $THREADS \
	    -o $WORKING_DIRECTORY/mantis_merged # $MANTIS_INDEX

