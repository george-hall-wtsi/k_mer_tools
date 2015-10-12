#######
####### NEEDS TO BE PASSED ABSOLUTE PATH TO REFERENCE AND REPEAT FILE
#######

REPEATS=$1
REPEATS_NAME=${REPEATS##*/}
REPEATS_NAME=${REPEATS_NAME%*.*}

PEAK_NUM=$2

MAIN_LOC=$3
RENAME_FASTQ_BIN=$MAIN_LOC"/../bin/rename_fastq"
SOAP_BIN=$MAIN_LOC"/../bin/SOAPdenovo-63mer"
GAP_CLOSER_BIN=$MAIN_LOC"/../bin/GapCloser"
SPADES_BIN="/nfs/users/nfs_z/zn1/src/SPAdes/SPAdes-3.5.0-Linux/bin/spades.py"

ASSEMBLER=$4 # Can be either 'soap' or 'spades'

WORKING_DIR=$PWD"/"$REPEATS_NAME"_reads"

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

mkdir "peak_"$PEAK_NUM
find . -maxdepth 1 -type f -exec mv {} ./"peak_"$PEAK_NUM/ \;

cd ..
