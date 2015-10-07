#######
####### NEEDS TO BE PASSED ABSOLUTE PATH TO REFERENCE AND REPEAT FILE
#######

REFERENCE=$1
REFERENCE_NAME=${REFERENCE##*/}
REFERENCE_NAME=${REFERENCE_NAME%*.*}

REPEATS=$2
REPEATS_NAME=${REPEATS##*/}
REPEATS_NAME=${REPEATS_NAME%*.*}

PEAK_NUM=$3

MAIN_LOC=$4
RENAME_FASTQ_BIN=$MAIN_LOC"/../bin/rename_fastq"
SMALT_BIN=$MAIN_LOC"/../bin/smalt-0.7.4"
SOAP_BIN=$MAIN_LOC"/../bin/SOAPdenovo-63mer"
GAP_CLOSER_BIN=$MAIN_LOC"/../bin/GapCloser"
SPADES_BIN="/nfs/users/nfs_z/zn1/src/SPAdes/SPAdes-3.5.0-Linux/bin/spades.py"

ASSEMBLER=$5 # Can be either 'soap' or 'spades'

WORKING_DIR=$PWD"/"$REPEATS_NAME"_reads"

HASH_NAME=$REFERENCE_NAME".hash"
HASH_LOCATION=$PWD"/"$HASH_NAME

# Generate hash of reference if requried (but hopefully will already be there)
if [ ! -f $HASH_NAME".smi" ] || [ ! -f $HASH_NAME".sma" ]; then
	$MAIN_LOC"/scripts/generate_hash.sh" $HASH_LOCATION $REFERENCE $MAIN_LOC
fi

ASSEMBLY_CONFIG_LOCATION=$MAIN_LOC"/scripts/assembly_config"

cd $WORKING_DIR

# Change assembly_config and sh.run-soap (maybe make k = 37)
K_SIZE=31
NUM_PROCESSORS=20

if [ $ASSEMBLER = "soap" ]; then
	$SOAP_BIN all -s $ASSEMBLY_CONFIG_LOCATION -K $K_SIZE -k $K_SIZE -o "k"$K_SIZE -p $NUM_PROCESSORS > "k"$K_SIZE".all.err"
	$GAP_CLOSER_BIN -o "k"$K_SIZE".fasta" -t $NUM_PROCESSORS -b $ASSEMBLY_CONFIG_LOCATION -a "k"$K_SIZE".scafSeq" > "k"$K_SIZE"-new.gf.err"
fi

if [ $ASSEMBLER = "spades" ]; then
	$SPADES_BIN --s1 "peak_"$PEAK_NUM"_k_mers-read.fastq" -t $NUM_PROCESSORS \
		-o "out-spades"
	mv "out-spades/contigs.fasta" "k"$K_SIZE".fasta"
fi

$RENAME_FASTQ_BIN -name contig -len 200 "k"$K_SIZE".fasta" "k"$K_SIZE"-2.fastq"

$SMALT_BIN map -m 200 -f ssaha -n $NUM_PROCESSORS -O -d 0 $HASH_LOCATION "k"$K_SIZE"-2.fastq" > "peak_"$PEAK_NUM"_map"
grep "alignment:S:00" "peak_"$PEAK_NUM"_map" > "grepped"
mv "grepped" "peak_"$PEAK_NUM"_map"

mkdir "peak_"$PEAK_NUM
find . -maxdepth 1 -type f -exec mv {} ./"peak_"$PEAK_NUM/ \;

cd ..
