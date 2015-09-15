#! /bin/bash -x

PEAK_NUM=$1
cat *.tmp.dump.fasta > "/lustre/scratch110/sanger/gh10/Data/k_mer_locations_and_words/peak_"$PEAK_NUM"_k_mers.fasta"
rm *.tmp.dump.fasta
cat *.tmp.ssaha > "/lustre/scratch110/sanger/gh10/Data/k_mer_locations_and_words/peak_"$PEAK_NUM".ssaha"
rm *.tmp.ssaha
echo "Finished concatenating for peak number "$PEAK_NUM
