REFERENCE=$1
PEAK_NUM=$2
NAME=${REFERENCE##*/}
NAME=${NAME%*.*}
echo $NAME

RENAME_FASTQ_BIN="/software/hpag/icas/0.61/icas/bin/rename_fastq"
SMALT_BIN="/software/hpag/bin/smalt-0.7.4"
SOAP_BIN="/nfs/users/nfs_z/zn1/src/soap/SOAPdenovo2-bin-LINUX-generic-r240/SOAPdenovo-63mer"
GAP_CLOSER_BIN="/software/hpag/icas/0.61/icas/bin/GapCloser"
SSAHA_SHRED_BIN="/nfs/users/nfs_z/zn1/bin/ssaha_shred"

SUBDIR=$NAME"_k_mer_words_locations"

ASSEMBLY_CONFIG_LOCATION="/nfs/users/nfs_g/gh10/Documents/Repositories/k_mer_tools/src/scripts/assembly_config"

# Change assembly_config and sh.run-soap (maybe make k = 37)
K_SIZE=33

cd $SUBDIR

### USED TO BE SEPERATE FILE CALLED sh.run-soap
$SOAP_BIN all -s $ASSEMBLY_CONFIG_LOCATION -K $K_SIZE -k $K_SIZE -o "k"$K_SIZE -p 20 > "k"$K_SIZE".all.err"
$GAP_CLOSER_BIN -o "k"$K_SIZE".fasta" -t 20 -b $ASSEMBLY_CONFIG_LOCATION -a "k"$K_SIZE".scafSeq" > "k"$K_SIZE"-new.gf.err"
### END OF OLD sh.run-soap

$RENAME_FASTQ_BIN -name "contig" -len 100 -num 2 "k"$K_SIZE".fasta" "k"$K_SIZE"-2.fastq"

cd ..

if [ ! -f $NAME"_hash_file"* ]; then
	$SMALT_BIN index -k 17 -s 17 $NAME"_hash_file" $REFERENCE
fi

HASH_LOC=$PWD"/"$NAME"_hash_file"
echo $HASH_LOC

cd $SUBDIR

$RENAME_FASTQ_BIN -name "contig" -len 100 -num 2 "k"$K_SIZE"-2.fasta" "k"$K_SIZE"-2.fastq"

$SMALT_BIN map -m 20 -f ssaha -n 4 -O -d 0 $HASH_LOC "k"$K_SIZE"-2.fastq"

SHRED_SIZE=500
$SSAHA_SHRED_BIN -rlength $SHRED_SIZE $REFERENCE $NAME"-shred-"$SHRED_SIZE"bp.fasta"
$SMALT_BIN -m 20 -f ssaha -n 4 -O -d 0 $HASH_LOC $NAME"-shred-"$SHRED_SIZE"bp.fasta"
