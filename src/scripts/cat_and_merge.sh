#! /bin/bash -x

PEAK_NUM=$1
NAME=$2

if [ ! -d $NAME"_k_mer_words_locations" ]; then
	echo "Creating k_mer_locations directory locally"
	mkdir $NAME"_k_mer_words_locations"
fi

TEMP_FILE=$NAME"_temp"
cat $TEMP_FILE"/*.tmp.dump.fasta" > $NAME"_k_mer_locations_and_words/peak_"$PEAK_NUM"_k_mers.fasta"
rm $TEMP_FILE"/*.tmp.dump.fasta"
cat $TEMP_FILE"/*.tmp.ssaha" > $NAME"_k_mer_locations_and_words/peak_"$PEAK_NUM".ssaha"
rm $TEMP_FILE"/*.tmp.ssaha"
echo "Finished concatenating for peak number "$PEAK_NUM
