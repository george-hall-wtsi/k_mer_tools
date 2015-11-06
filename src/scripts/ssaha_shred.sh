#!/bin/bash -x


################################################################################
# Copyright (c) 2015 Genome Research Ltd. 
# 
# Author: George Hall gh10@sanger.ac.uk 
# 
# This program is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software 
# Foundation; either version 3 of the License, or (at your option) any later 
# version. 
# 
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
# details. 
# 
# You should have received a copy of the GNU General Public License along with 
# this program. If not, see <http://www.gnu.org/licenses/>. 
################################################################################


REFERENCE=$1
REFERENCE_NAME=${REFERENCE##*/}
REFERENCE_NAME=${REFERENCE_NAME%*.*}

FILE_NAME=$2
MAIN_LOC=$3

HASH_NAME=$REFERENCE_NAME".hash"
HASH_LOCATION=$PWD"/"$HASH_NAME

WORKING_DIR=$FILE_NAME"_reads"
cd $WORKING_DIR

SMALT_BIN=$MAIN_LOC"/../bin/smalt-0.7.4"
SSAHA_SHRED_BIN=$MAIN_LOC"/../bin/ssaha_shred"

NUM_CPUS=20

SHRED_SIZE=100
$SSAHA_SHRED_BIN -rlength $SHRED_SIZE $REFERENCE $REFERENCE_NAME"-shred-"$SHRED_SIZE"bp.fasta"

# Generate hash of reference if requried (but hopefully will already be there)
if [ ! -f $HASH_LOCATION".smi" ] || [ ! -f $HASH_LOCATION".sma" ]; then
	$MAIN_LOC"/scripts/generate_hash.sh" $HASH_LOCATION $REFERENCE $MAIN_LOC
fi

$SMALT_BIN map -m 20 -f ssaha -n $NUM_CPUS -O -d 0 \
	$HASH_LOCATION $REFERENCE_NAME"-shred-"$SHRED_SIZE"bp.fasta" > "shred_map"

cd ..

