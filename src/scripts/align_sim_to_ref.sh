#! /bin/bash -x

PEAK_NUM=$1
NAME=$2
REFERENCE=$3

REFERENCE=${REFERENCE##*/ }
echo $REFERENCE

WORKING_DIR=${NAME%*.*}
WORKING_DIR=${WORKING_DIR##*/}"_reads"

HASH_LOCATION=${REFERENCE%*.*}".hash"

# Generate hash of reference if requried (but hopefully will already be there)
if [ ! -f $HASH_LOCATION".smi" ] || [ ! -f $HASH_LOCATION".sma" ]; then
	/nfs/users/nfs_g/gh10/Documents/Repositories/k_mer_tools/src/scripts/generate_hash.sh $HASH_LOCATION $REFERENCE
fi

/software/hpag/bin/smalt-0.7.4 map -m 20 -f ssaha -n 4 -O -d -0 $HASH_LOCATION $WORKING_DIR"/peak_"$PEAK_NUM"_k_mers-read.fastq" > $WORKING_DIR"/peak_"$PEAK_NUM"_occs.ssaha" 
sed -i "s|$| $PEAK_NUM|" $WORKING_DIR"/peak_"$PEAK_NUM"_occs.ssaha"  
