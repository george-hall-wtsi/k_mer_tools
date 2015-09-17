#! /bin/bash -x

# USING THE MOST RECENT mer_counts.jf, this script generates all k-mers occurring with the specified number. Jellyfish then dumps these
# into a .fasta format, and then they are located on the reference genome using Smalt.
# This script should be run in the directory containing the reference 

NAME=$1
LOWER_LIM=$2
UPPER_LIM=$3
PEAK_NUM=$4

WITHOUT_EXTENSION=${NAME%.*}
RENAME_FASTQ_BIN="/software/hpag/icas/0.61/icas/bin/rename_fastq"

if [ ! -d $NAME"_reads" ]; then
	mkdir $NAME"_reads"
fi

/nfs/users/nfs_g/gh10/src/jellyfish-2.2.3/bin/jellyfish dump -L $LOWER_LIM -U $UPPER_LIM -ct $WITHOUT_EXTENSION"_mer_counts.jf" > $NAME"_reads/"$PEAK_NUM"_words.tmp.fasta"

cd $NAME"_reads"
cat $PEAK_NUM"_words.tmp.fasta" | awk '{print ">reads \n" $1}' > "peak_"$PEAK_NUM"_k_mers-read.fasta"
$RENAME_FASTQ_BIN -name kmer-read -len 10 "peak_"$PEAK_NUM"_k_mers-read.fasta" "peak_"$PEAK_NUM"_k_mers-read.fastq"
rm $PEAK_NUM"_words.tmp.fasta"

cd ..
