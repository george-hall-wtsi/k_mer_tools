import subprocess32

with open("temp", "r") as f:
	data = f.readlines()
temp = [line.split() for line in data]
data = temp
prev = data[0]
nxt = data[2]
new_data = [prev]

for line in data[1:]:
	if line[2] == prev[2] and line[3] == prev[3] and (int(prev[6]) < int(line[6]) < int(prev[7])):
		pass
#	if line[8] == 'C' and line[2] == next[2] and line[3] == next[3] and 
	else:
		new_data.append(line)
		prev = line

with open("outfile", "w") as out:	
	for line in sorted(new_data, key = lambda x: int(x[6])):
		out.write(" ".join(x for x in line[:3]) + " " + " ".join(str(x).rjust(10) for x in line[3:8]) + " " + " ".join(str(x) for x in line[8:]) + "\n")

store_dict = {}
with open("outfile", "r") as f:
	for line in f.readlines():
		line = line.split()
		try:
			store_dict[line[2]].append(int(line[7]))
		except KeyError:
			store_dict[line[2]] = [int(line[7])]

	for x in sorted(store_dict.keys()):
		y = store_dict[x]
		new_suffix = str(int(x.split("_X")[1]) + 200)
		second = x[:-len(new_suffix)] + new_suffix
		try:
			next_shred = store_dict[second]
			for i in y:
				if any(abs(i - j) <= 200 for j in next_shred):
					y.remove(i)	

		except KeyError:
			pass


	out_dict = {}
	for i in store_dict.keys():
		try:
			out_dict[len(store_dict[i])] += 1
		except KeyError:
			out_dict[len(store_dict[i])] = 1

	print out_dict


print "Done!"
