#! /bin/bash -x

PEAK_NUM=$1
cat *.tmp.dump.fasta > "./k_mer_locations_and_words/peak_"$PEAK_NUM"_k_mers.fasta"
rm *.tmp.dump.fasta
cat ./k_mer_locations_and_words/*.tmp.ssaha > "./k_mer_locations_and_words/peak_"$PEAK_NUM".ssaha"
rm ./k_mer_locations_and_words/*.tmp.ssaha
echo "Finished concatenating for peak number "$PEAK_NUM
