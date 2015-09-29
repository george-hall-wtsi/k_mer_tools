#! /bin/bash -x

REFERENCE=$1
WORKING_DIR=$2
CONTIG_MASK_BIN=$3"/../bin/contig_mask"
IN_PATH=$4

IN_NAME=${IN_PATH##*/}
IN_NAME=${IN_NAME%*.*}

cd $WORKING_DIR

if [ ! -d "Masked Repeats" ]; then
	mkdir "Masked Repeats"
fi

$CONTIG_MASK_BIN $REFERENCE $IN_PATH $WORKING_DIR"/Masked Repeats/"$IN_NAME"_mask"

cd ..
