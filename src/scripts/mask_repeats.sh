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


REFERENCE=$1
WORKING_DIR=$2
CONTIG_MASK_BIN=$3"/../bin/contig_mask"
IN_PATH=$4

IN_NAME=${IN_PATH##*/}
IN_NAME=${IN_NAME%*.*}

cd $WORKING_DIR

if [ ! -d "Masked Repeats" ]; then
	mkdir "Masked Repeats"
fi

$CONTIG_MASK_BIN $REFERENCE $IN_PATH $WORKING_DIR"/Masked Repeats/"$IN_NAME"_mask"

cd ..
