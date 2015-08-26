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
	
	with open("/lustre/scratch110/sanger/gh10/Code/k_mer_scripts/hgram_data/" + file_name \
	+ ".hgram","w") as to_write:
	
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

