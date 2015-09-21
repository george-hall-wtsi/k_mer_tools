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

		store_dict[n] = [(int(repeat.split("00000000_X")[1]), int(repeat.split("00000000_X")[1]) + 100) for repeat in store_lst]
		tmp[n] = []

		start = None
		for (first, last) in sorted(store_dict[n]):
			if start is None:
				start = first
			elif prev_end != first:
				tmp[n].append((start, prev_end))
				start = first

			prev_end = last
	print sorted(store_dict[2])
	for key in tmp.keys():
		print tmp[key] 


assess("/lustre/scratch110/sanger/gh10/ecoli/Escherichiacoli-K-12-simu-random_both")
