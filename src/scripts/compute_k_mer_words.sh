#! /bin/bash -x


################################################################################
# Copyright (c) 2015 Genome Research Ltd. 
#  
# Author: George Hall <gh10@sanger.ac.uk> 
# 
# This file is part of K-mer Toolkit. 
# 
# K-mer Toolkit is free software: you can redistribute it and/or modify it under 
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


# Using the relevant mer_counts_*.jf, this script generates all k-mers occurring with the 
# specified number. Jellyfish then dumps these into a .fasta format, and then they are located 
# on the reference genome using Smalt. This script should be run in the directory containing 
# the reference 

NAME=$1
LOWER_LIM=$2
UPPER_LIM=$3
PEAK_NUM=$4
MAIN_LOC=$5
K_SIZE=$6

WITHOUT_EXTENSION=${NAME%*.*}
WITHOUT_EXTENSION=${WITHOUT_EXTENSION##*/}

RENAME_FASTQ_BIN=$MAIN_LOC"/../bin/rename_fastq"

if [ ! -d $WITHOUT_EXTENSION"_reads" ]; then
	mkdir $WITHOUT_EXTENSION"_reads"
fi

$MAIN_LOC"/../bin/jellyfish" dump -L $LOWER_LIM -U $UPPER_LIM \
	-ct $WITHOUT_EXTENSION"_mer_counts_"$K_SIZE".jf" > \
	$WITHOUT_EXTENSION"_reads/"$PEAK_NUM"_words.tmp.fasta"

cd $WITHOUT_EXTENSION"_reads"
cat $PEAK_NUM"_words.tmp.fasta" | awk '{print ">reads \n" $1}' > \
	"peak_"$PEAK_NUM"_k_mers-read.fasta"
$RENAME_FASTQ_BIN -name kmer-read "peak_"$PEAK_NUM"_k_mers-read.fasta" "peak_"$PEAK_NUM"_k_mers-read.fastq"

rm $PEAK_NUM"_words.tmp.fasta"

cd ..
