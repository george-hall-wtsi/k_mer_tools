def check_equal(lst):
		return lst[1:] == lst[:-1]

def assess(file_path):
	file_path = file_path.split(".")[0]
	with open(file_path + "_reads/" + "shred_map", "r+") as f:
		data = [line.split() for line in f.readlines()]

	new_data = [item[2] for item in data]

	store_dict = {} 
	tmp = {}
	for n in xrange(2,15):
		store_lst = []
		for i  in xrange(1, len(new_data) - n):
			if check_equal(new_data[i:i+n]) and (new_data[i] != new_data[i-1]) and (new_data[i] != new_data[i+n]):
				store_lst.append(new_data[i])

		store_dict[n] = []
		for repeat in store_lst:
			new_suffix = str(int(repeat.split("_X")[1]) + 200)
			second = repeat[:-len(new_suffix)] + new_suffix
			store_dict[n].append((repeat, second))

		tmp[n] = []

		start = None
		for (first, last) in sorted(store_dict[n]):
			if start is None:
				start = first
			elif prev_end != first:
				tmp[n].append((start, prev_end))
				start = first

			prev_end = last

	for key in tmp.keys():
		print "Count:" , key , tmp[key] , "\n"

	for key in tmp.keys():
		print key, len(tmp[key])

assess("/lustre/scratch110/sanger/gh10/c_elegans_kmers/c_elegans_ref-simu-random_both")
