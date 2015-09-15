#! /bin/bash -x

PEAK_NUM=$1
cat /lustre/scratch110/sanger/gh10/Data/tmp/*.tmp.dump.fasta > "/lustre/scratch110/sanger/gh10/Data/k_mer_locations_and_words/peak_"$PEAK_NUM"_k_mers.fasta"
rm /lustre/scratch110/sanger/gh10/Data/tmp/*.tmp.dump.fasta
cat /lustre/scratch110/sanger/gh10/Data/tmp/*.tmp.ssaha > "/lustre/scratch110/sanger/gh10/Data/k_mer_locations_and_words/peak_"$PEAK_NUM".ssaha"
rm /lustre/scratch110/sanger/gh10/Data/tmp/*.tmp.ssaha
echo "Finished concatenating for peak number "$PEAK_NUM
