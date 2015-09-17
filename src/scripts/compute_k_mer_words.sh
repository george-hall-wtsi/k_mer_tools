#! /bin/bash -x

# USING THE MOST RECENT mer_counts.jf, this script generates all k-mers occurring with the specified number. Jellyfish then dumps these
# into a .fasta format, and then they are located on the reference genome using Smalt.
# This script should be run in the directory containing the reference 

NAME=$1
LOWER_LIM=$2
UPPER_LIM=$3
PEAK_NUM=$4

WITHOUT_EXTENSION=${NAME%*.*}
WITHOUT_EXTENSION=${WITHOUT_EXTENSION##*/}

RENAME_FASTQ_BIN="/software/hpag/icas/0.61/icas/bin/rename_fastq"

if [ ! -d $WITHOUT_EXTENSION"_reads" ]; then
	mkdir $WITHOUT_EXTENSION"_reads"
fi

/nfs/users/nfs_g/gh10/src/jellyfish-2.2.3/bin/jellyfish dump -L $LOWER_LIM -U $UPPER_LIM -ct $WITHOUT_EXTENSION"_mer_counts.jf" > $WITHOUT_EXTENSION"_reads/"$PEAK_NUM"_words.tmp.fasta"

cd $WITHOUT_EXTENSION"_reads"
cat $PEAK_NUM"_words.tmp.fasta" | awk '{print ">reads \n" $1}' > "peak_"$PEAK_NUM"_k_mers-read.fasta"
$RENAME_FASTQ_BIN -name kmer-read -len 10 "peak_"$PEAK_NUM"_k_mers-read.fasta" "peak_"$PEAK_NUM"_k_mers-read.fastq"
rm $PEAK_NUM"_words.tmp.fasta"

cd ..
