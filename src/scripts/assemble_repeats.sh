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

RENAME_FASTQ_BIN="/software/hpag/icas/0.61/icas/bin/rename_fastq"
SMALT_BIN="/software/hpag/bin/smalt-0.7.4"
SOAP_BIN="/nfs/users/nfs_z/zn1/src/soap/SOAPdenovo2-bin-LINUX-generic-r240/SOAPdenovo-63mer"
GAP_CLOSER_BIN="/software/hpag/icas/0.61/icas/bin/GapCloser"

WORKING_DIR=$PWD"/"$REPEATS_NAME"_reads"

HASH_NAME=$REFERENCE_NAME".hash"
HASH_LOCATION=$PWD"/"$HASH_NAME

# Generate hash of reference if requried (but hopefully will already be there)
if [ ! -f $HASH_NAME".smi" ] || [ ! -f $HASH_NAME".sma" ]; then
	/nfs/users/nfs_g/gh10/Documents/Repositories/k_mer_tools/src/scripts/generate_hash.sh $HASH_LOCATION $REFERENCE
fi

ASSEMBLY_CONFIG_LOCATION="/nfs/users/nfs_g/gh10/Documents/Repositories/k_mer_tools/src/scripts/assembly_config"

cd $WORKING_DIR

# Change assembly_config and sh.run-soap (maybe make k = 37)
K_SIZE=31
NUM_PROCESSORS=20

### USED TO BE SEPERATE FILE CALLED sh.run-soap
$SOAP_BIN all -s $ASSEMBLY_CONFIG_LOCATION -K $K_SIZE -k $K_SIZE -o "k"$K_SIZE -p $NUM_PROCESSORS > "k"$K_SIZE".all.err"
$GAP_CLOSER_BIN -o "k"$K_SIZE".fasta" -t $NUM_PROCESSORS -b $ASSEMBLY_CONFIG_LOCATION -a "k"$K_SIZE".scafSeq" > "k"$K_SIZE"-new.gf.err"
### END OF OLD sh.run-soap

$RENAME_FASTQ_BIN -name contig -len 200 "k"$K_SIZE".fasta" "k"$K_SIZE"-2.fastq"
#$RENAME_FASTQ_BIN -name contig -len 100 "k"$K_SIZE"-2.fastq" "k"$K_SIZE"-2.fastq"

$SMALT_BIN map -m 20 -f ssaha -n $NUM_PROCESSORS -O -d 0 $HASH_LOCATION "k"$K_SIZE"-2.fastq" > "peak_"$PEAK_NUM"_map"

mkdir "peak_"$PEAK_NUM
find . -type f -maxdepth 1 -exec mv {} ./"peak_"$PEAK_NUM/ \;

cd ..
