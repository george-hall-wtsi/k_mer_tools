#! /bin/bash -x

REFERENCE=$1
REFERENCE_NAME=${REFERENCE##*/}
REFERENCE_NAME=${REFERENCE_NAME%*.*}

REPEATS=$2
REPEATS_NAME=${REPEATS##*/}
REPEATS_NAME=${REPEATS_NAME%*.*}

PEAK_NUM=$3 
MAIN_LOC=$4

SMALT_BIN=$MAIN_LOC"/../bin/smalt-0.7.4"

HASH_NAME=$REFERENCE_NAME".hash"
HASH_LOCATION=$PWD"/"$HASH_NAME

WORKING_DIR=$PWD"/"$REPEATS_NAME"_reads"

K_SIZE=31
NUM_PROCESSORS=20

cd $WORKING_DIR"/peak_"$PEAK_NUM

# Generate hash of reference if requried (but hopefully will already be there)
if [ ! -f $HASH_NAME".smi" ] || [ ! -f $HASH_NAME".sma" ]; then
	$MAIN_LOC"/scripts/generate_hash.sh" $HASH_LOCATION $REFERENCE $MAIN_LOC
fi

$SMALT_BIN map -m 200 -f ssaha -n $NUM_PROCESSORS -O -d 0 $HASH_LOCATION "k"$K_SIZE"-2.fastq" > "peak_"$PEAK_NUM"_map"
grep "alignment:S:00" "peak_"$PEAK_NUM"_map" > "grepped"
mv "grepped" "peak_"$PEAK_NUM"_map"

cd ../..
