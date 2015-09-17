#! /bin/bash -x

PEAK_NUM=$1
NAME=$2
REFERENCE=$3

REFERENCE=${REFERENCE##*/ }
echo $REFERENCE


HASH_LOCATION=$REFERENCE".hash"

# Generate hash of reference if requried (but hopefully will already be there)
if [ ! -f $HASH_LOCATION".smi" ] || [ ! -f $HASH_LOCATION".sma" ]; then
	/software/hpag/bin/smalt-0.7.4 index -k 17 -s 17 $HASH_LOCATION $REFERENCE
fi

/software/hpag/bin/smalt-0.7.4 map -m 20 -f ssaha -n 4 -O -d -0 $HASH_LOCATION $NAME"_reads/peak_"$PEAK_NUM"_k_mers-read.fastq" > $NAME"_reads/peak_"$PEAK_NUM"_occs.ssaha" 
sed -i "s|$| $PEAK_NUM|" $NAME"_reads/peak_"$PEAK_NUM"_occs.ssaha"  

