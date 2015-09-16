NAME=$1
REFERENCE=$2
PEAK_NUM=$3

RENAME_FASTQ_BIN="~zn1/bin/rename_fastq"
SMALT_BIN="/software/hpag/bin/smalt-0.7.4"
SOAP_BIN="/nfs/users/nfs_z/zn1/src/soap/SOAPdenovo2-bin-LINUX-generic-r240/SOAPdenovo-63mer"
GAP_CLOSER_BIN="~zn1/bin/GapCloser"
SSAHA_SHRED_BIN="/nfs/users/nfs_z/zn1/bin/ssaha_shred"

cd $NAME"_k_mer_words_locations"
"k_"$PEAK_NUM"_k_mers.fasta" | awk '{print ">reads \n" $1}' > "peak_"$PEAK_NUM"_k_mers-read.fasta"
$RENAME_FASTQ_BIN -name kmer-read -len 10 "peak_"$PEAK_NUM"_k_mers-read.fasta" "peak_"$PEAK_NUM"_k_mers-read.fastq"

# Change assembly_config and sh.run-soap (maybe make k = 37)

### USED TO BE SEPERATE FILE CALLED sh.run-soap
$SOAP_BIN all -s assembly_config -K 33 -k 33 -o k33 -p 20 > k33.all.err
$GAP_CLOSER_BIN -o k33.fasta -t 20 -b assembly_config -a k33.scafSeq > k33-new.gf.err 
### END OF OLD sh.run-soap

$RENAME_FASTQ_BIN -name contig -len 100 -num 2 k37.fasta k37-2.fastq

if [ ! -f $NAME"_hash_file"* ]; then
	$SMALT_BIN index -k 17 -s 17 $NAME"_hash_file" $REFERENCE
fi

$RENAME_FASTQ_BIN -name contig -len 100 -num 2 k37-2.fasta k37-2.fastq

$SMALT_BIN map -m 20 -f ssaha -n 4 -O -d 0 $NAME"_hash_file" k37-2.fastq

$SSAHA_SHRED_BIN -rlength 500 $REFERENCE $NAME"-shred-500bp.fasta"

$SMALT_BIN -m 20 -f ssaha -n 4 -O -d 0 $NAME"_hash_file" $NAME"-shred-500bp.fasta"
