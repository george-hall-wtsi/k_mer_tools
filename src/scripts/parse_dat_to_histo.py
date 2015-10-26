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


def parse(input_file_path, k_mer_size):

	"""
	Takes file stored at 'input_file_path' in .dat format (i.e. formatted such as:)
	
	hist:     0            0
	hist:     1      7919008
	hist:     2       404988
	hist:     3       109952
	
	and saves the data in the same format as is returned by compute_histogram(),
	for example if the above file were passed, it would be saved in the format:
	
	1 7919008
	2 404988
	3 109952
	
	A file passed in as foo.* will be saved in '/lustre/scratch110/sanger/gh10/hgram_data'
	as foo_nmer.hgram, where n = k_mer_size.
	
	It is worth noting that, as the k-mers corresponding to the file which is being read have 
	already been counted, the k_mer_size variable is used merely to catalogue the k-mer size 
	used. It is the user's responsibility to  record this number. 

	"""
	
	file_name = input_file_path.split("/")[-1].split(".")[0] + "_" + str(k_mer_size) + "mer"
	
	with open("/".join(input_file_path.split("/")[:-1]) + "/" + file_name + ".hgram", 'w') \
		as to_write: 
	
		with open(input_file_path,"r") as f:
			file_lines = f.readlines()
	
		line_count = sum(1 for line in file_lines)
	
		iCount = 0
		for line in file_lines:
			storeStr = ""
			iCount += 1
			
			for word in line.split(" "):
				if len(word) == 0 or str(word) in ['0', '0\n']:
					pass
				elif word[0].isdigit() and storeStr == "":
					storeStr += str(word) + " "
				elif word[0].isdigit():
					storeStr += str(word)[:-1]
			if iCount == line_count:
				storeStr += "0" 

			if len(storeStr) > 0:
				to_write.write(storeStr + "\n")

