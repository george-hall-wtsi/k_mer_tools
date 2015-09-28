#! /bin/bash -x

# ref path, file name, peak number, location

CONTIG_MASK_BIN=$4"/../bin/contig_mask"

WORKING_DIR=$2
WORKING_DIR=${WORKING_DIR##*/}
WORKING_DIR=${WORKING_DIR%*.*}
WORKING_DIR=$WORKING_DIR"_reads"

PEAK_NUM=$3

cd $WORKING_DIR

if [ ! -d "Masked Repeats" ]; then
	mkdir "Masked Repeats"
fi

cd "peak_"$PEAK_NUM
$CONTIG_MASK_BIN $1 "peak_"$PEAK_NUM"_map" "peak_"$PEAK_NUM"_mask"
mv "peak_"$PEAK_NUM"_mask" ../"Masked Repeats"

cd ../..
