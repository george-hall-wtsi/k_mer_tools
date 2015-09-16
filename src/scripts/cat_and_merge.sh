#! /bin/bash -x

PEAK_NUM=$1
NAME=$2
RENAME_FASTQ_BIN="/software/hpag/icas/0.61/icas/bin/rename_fastq"

if [ ! -d $NAME"_k_mer_words_locations" ]; then
	echo "Creating k_mer_locations directory locally"
	mkdir $NAME"_k_mer_words_locations"
fi

TEMP_DIR=$NAME"_temp"
cd $TEMP_DIR
cat *.tmp.dump.fasta > "peak_"$PEAK_NUM"_k_mers.fasta"
mv "peak_"$PEAK_NUM"_k_mers.fasta" ../$NAME"_k_mer_words_locations"
rm *.tmp.dump.fasta
cat *.tmp.ssaha > "peak_"$PEAK_NUM".ssaha"
mv "peak_"$PEAK_NUM".ssaha" ../$NAME"_k_mer_words_locations"
rm *.tmp.ssaha
cd ..
rm -rf $TEMP_DIR
echo "Finished concatenating for peak number "$PEAK_NUM

echo "Renaming for peak number "$PEAK_NUM
cd $NAME"_k_mer_words_locations"

cat "peak_"$PEAK_NUM"_k_mers.fasta" | awk '{print ">reads \n" $1}' > "peak_"$PEAK_NUM"_k_mers-read.fasta"
$RENAME_FASTQ_BIN -name kmer-read -len 10 "peak_"$PEAK_NUM"_k_mers-read.fasta" "peak_"$PEAK_NUM"_k_mers-read.fastq"
echo "Finished renaming for peak number "$PEAK_NUM

cd ..
