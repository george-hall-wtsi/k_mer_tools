################################################################################
# Copyright (c) 2015 Genome Research Ltd. 
# 
# Author: George Hall gh10@sanger.ac.uk 
# 
# This program is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software 
# Foundation; either version 3 of the License, or (at your option) any later 
# version. 
# 
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
# details. 
# 
# You should have received a copy of the GNU General Public License along with 
# this program. If not, see <http://www.gnu.org/licenses/>. 
################################################################################


REPEATS=$1
REPEATS_NAME=${REPEATS##*/}
REPEATS_NAME=${REPEATS_NAME%*.*}

PEAK_NUM=$2

MAIN_LOC=$3
RENAME_FASTQ_BIN=$MAIN_LOC"/../bin/rename_fastq"

ASSEMBLER=$4 # Can be either 'soap' or 'spades'

WORKING_DIR=$PWD"/"$REPEATS_NAME"_reads"

ASSEMBLY_CONFIG_LOCATION=$MAIN_LOC"/scripts/assembly_config"

cd $WORKING_DIR

K_SIZE=$5
NUM_PROCESSORS=$6

ASSEMBLER_BIN=$7
GAP_CLOSER_BIN=$8

if [ $ASSEMBLER = "soap" ]; then
	$ASSEMBLER_BIN all -s $ASSEMBLY_CONFIG_LOCATION -K $K_SIZE -k $K_SIZE -o "k"$K_SIZE \
		-p $NUM_PROCESSORS > "k"$K_SIZE".all.err"
	$GAP_CLOSER_BIN -o "k"$K_SIZE".fasta" -t $NUM_PROCESSORS -b $ASSEMBLY_CONFIG_LOCATION \
		-a "k"$K_SIZE".scafSeq" > "k"$K_SIZE"-new.gf.err"
	rm -f "k"$K_SIZE".bubbleInScaff" "k"$K_SIZE".contigPosInscaff" "k"$K_SIZE".fasta.fill" \
		"k"$K_SIZE".links" "k"$K_SIZE".peGrads" "k"$K_SIZE".readInGap.gz" \
		"k"$K_SIZE".scaf_gap" "k"$K_SIZE".updated.edge" "k"$K_SIZE".all.err" \
		"k"$K_SIZE".contig" "k"$K_SIZE".edge.gz" "k"$K_SIZE".gapSeq" \
		"k"$K_SIZE".newContigIndex" "k"$K_SIZE".preArc" "k"$K_SIZE".readOnContig.gz" \
		"k"$K_SIZE".scafSeq" "k"$K_SIZE".vertex" "k"$K_SIZE".Arc" "k"$K_SIZE".ContigIndex" \
		"k"$K_SIZE".kmerFreq" "k"$K_SIZE"-new.gf.err" "k"$K_SIZE".preGraphBasic" \
		"k"$K_SIZE".scaf" "k"$K_SIZE".scafStatistics"  

fi

if [ $ASSEMBLER = "spades" ]; then
	$ASSEMBLER_BIN --s1 "peak_"$PEAK_NUM"_k_mers-read.fastq" -t $NUM_PROCESSORS \
		-o "out-spades"
	mv "out-spades/contigs.fasta" "k"$K_SIZE".fasta"
	rm -rf "out-spades"
fi

rm "peak_"$PEAK_NUM"_k_mers-read.fasta"

$RENAME_FASTQ_BIN -name contig -len 200 "k"$K_SIZE".fasta" "contigs.fastq"
rm "k"$K_SIZE".fasta"

mkdir "peak_"$PEAK_NUM
find . -maxdepth 1 -type f -exec mv {} ./"peak_"$PEAK_NUM/ \;

cd ..
