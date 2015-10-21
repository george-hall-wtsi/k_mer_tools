#! /bin/bash -x


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


# Takes .fasta file, desired coverage, read lengthn isnert size, and desired output name as 
# inputs and returns two reads, which simulate 100x coverage reads of length 100 in both 
# forward and backward directions. Assumes that original .fasta reads are of length 250
# Call as script in directory containing .fasta file. Finally combines both reads into one file

INPUT_FILE=$1
COVERAGE=$2
READ_LENGTH=$3
INSERT_SIZE=$4
NAME=$5
MAIN_LOC=$6

$MAIN_LOC"/../bin/simulation_reads-randam2" -rlength $READ_LENGTH -cover $COVERAGE \
	-insert $INSERT_SIZE $INPUT_FILE "temp-simu-random.fastq"
cat temp-simu-random.fastq_*.fastq > temp-simu-random.fastq
rm temp-simu-random.fastq_*.fastq
$MAIN_LOC"/../bin/ssaha_reads" -file 22 "temp-simu-random.fastq" \
	$NAME-simu-random_1.fastq $NAME-simu-random_2.fastq
rm temp-simu-random.fastq
cat $NAME-simu-random_1.fastq $NAME-simu-random_2.fastq > $NAME-simu-random_both.fastq

echo "Finished simulating reads"
