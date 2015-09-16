NAME=$1
REFERENCE=$2
PEAK_NUM=$3
RENAME_FASTQ_BIN="~zn1/bin/rename_fastq"

cd $NAME"_k_mer_words_locations"
"k_"$PEAK_NUM"_k_mers.fasta" | awk '{print ">reads \n" $1}' > "peak_"$PEAK_NUM"_k_mers-read.fasta"
$RENAME_FASTQ_BIN -name kmer-read -len 10 "peak_"$PEAK_NUM"_k_mers-read.fasta" "peak_"$PEAK_NUM"_k_mers-read.fastq"

# Change the file of assembly_config
# Change the file of sh.run-soap
# Maybe make k = 37

sh sh.run-soap

$RENAME_FASTQ_BIN -name contig -len 100 -num 2 k37.fasta k37-2.fastq

if [ ! -f $NAME"_hash_file"* ]; then
	/software/hpag/bin/smalt-0.7.4 index -k 17 -s 17 $NAME"_hash_file" $REFERENCE
fi

$RENAME_FASTQ_BIN -name contig -len 100 -num 2 k37-2.fasta k37-2.fastq

/software/hpag/bin/smalt-0.7.4 map -m 20 -f ssaha -n 4 -O -d 0 $NAME"_hash_file" k37-2.fastq

/nfs/users/nfs_z/zn1/bin/ssaha_shred -rlength 500 $REFERENCE $NAME"-shred-500bp.fasta"

/software/hpag/bin/smalt-0.7.4 map -m 20 -f ssaha -n 4 -O -d 0 $NAME"_hash_file" $NAME"-shred-500bp.fasta"

