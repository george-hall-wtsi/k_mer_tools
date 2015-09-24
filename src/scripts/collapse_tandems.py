import collections

with open("temp", "r") as f:
	data = f.readlines()
data = [line.split() for line in data]
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

new_data.sort(key = lambda x: int(x[6]))

store_dict = {}
for line in new_data:
	try:
		store_dict[line[2]].append(int(line[7]))
	except KeyError:
		store_dict[line[2]] = [int(line[7])]

for (k, v) in collections.OrderedDict(sorted(store_dict.items())).iteritems():
	new_suffix = str(int(k.split("_X")[1]) + 200)
	follower = k[:-len(new_suffix)] + new_suffix
	try:
		next_shred = store_dict[follower]
		temp = v[:]
		for i in temp:
			if any(abs(i - j) <= 200 for j in next_shred):
				#v.remove(i)	
				pass
	except KeyError:
		pass

count_dict = {}

for key in store_dict.keys():
	try:
		count_dict[len(store_dict[key])] += 1
	except KeyError:
		count_dict[len(store_dict[key])] = 1

for k,v in count_dict.items():
	print k,v

print "Done!"

###	out_list.append(" ".join(x for x in line[:3]) + " " + " ".join(str(x).rjust(10) for x in line[3:8]) + " " + " ".join(str(x) for x in line[8:])) # will need to append new line character eventually
