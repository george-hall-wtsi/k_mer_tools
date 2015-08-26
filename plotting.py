#!/nfs/users/nfs_g/gh10/Documents/Code/Python/virtualenvs/venv/bin/python

import os.path
import subprocess32
import random

import matplotlib
import matplotlib.pyplot as plt

import parse_dat_to_histo as parse_data
import graph_settings


def compute_hist_from_fast(input_file_path, k_size):
	
	"""
	Uses Jellyfish to count k-mers of length k_size from the file stored at 
	'input_file_path'. This data is stored as a .hgram file in
	/lustre/scratch110/sanger/gh10/Code/k_mer_scripts/hgram_data/. 
	"""

	print "Computing hgram data for k = " + str(k_size) + " for first time"

	print "Reading data for k = " + str(k_size)

	# Counts occurences of k-mers of size "k-size" in "file_input":  
	subprocess32.call(["/nfs/users/nfs_g/gh10/src/jellyfish-2.2.3/bin/jellyfish", "count", 
	("-m " + str(k_size)), "-s 100M", "-t 20", "-C", input_file_path])

	print "Processing histogram for k = " + str(k_size)
	
	file_name = str(input_file_path.split("/")[-1].split(".")[0]) + "_" + str(k_size) + "mer"
	
	with open("/lustre/scratch110/sanger/gh10/Code/k_mer_scripts/hgram_data/" + 
	file_name + ".hgram","w") as out_file:
		# Computes histogram data and stores in "out_file"
		subprocess32.call(["/nfs/users/nfs_g/gh10/src/jellyfish-2.2.3/bin/jellyfish", "histo", 
		"mer_counts.jf"], stdout=out_file)
	
	print "Finished for k = " + str(k_size)

		
def find_extrema(hist_dict, window_size, alternate = False):
	
	"""
	Returns a dict with 2 keys (max and min) with the values for each of these keys being 
	tuples which correspond (occurrence, frequency) pairs which are either a maximum or a 
	minimum.
	"""
	
	i,j,k = (window_size), (window_size+1), (window_size+2)
	store_dict = {'Max' : [], 'Min' : []}

	if alternate == False:

		while k < (hist_dict.keys()[-1] - window_size + 1):
			
			if (all(hist_dict[i - x] < hist_dict[j] > hist_dict[k + x] \
			for x in xrange(0, (window_size + 1)))):
				store_dict['Max'].append((j + 1, hist_dict[j]))
			
			if (all(hist_dict[i - x] > hist_dict[j] < hist_dict[k + x] \
			for x in xrange(0, (window_size + 1)))):
				store_dict['Min'].append((j + 1, hist_dict[j]))

			i += 1
			j += 1
			k += 1
			
	else:	
		prev_min, prev_max = False, False
		
		while k < (hist_dict.keys()[-1] - window_size + 1):
			if prev_max == False and (all(hist_dict[i - x] < hist_dict[j] > hist_dict[k + x] \
			for x in xrange(0, (window_size + 1)))):
				store_dict['Max'].append((j + 1,hist_dict[j]))
				prev_min, prev_max = False, True
			
			elif prev_min == False and (all(hist_dict[i - x] > hist_dict[j] < hist_dict[k + x] \
			for x in xrange(0, (window_size + 1)))):			
				store_dict['Min'].append((j + 1,hist_dict[j]))
				prev_min, prev_max = True, False
		
			i += 1
			j += 1
			k += 1
		
	return store_dict
	
	
def calculate_modes(hist_dict, n):
	"""Takes a hist_dict as input and returns a list containing its first n modes. """
	hist_dict_augmented = format_data(hist_dict)
	modes = []
	for x in reversed(sorted(find_extrema(hist_dict_augmented, window_size = 15, 
	alternate = False)['Max'], key = lambda pair: pair[1])):
		if len(modes) < n:
			modes.append(x)
		else:
			break

	return modes
	
	
def format_data(hist_dict):

	"""
	This function is required when, for example, simulated data is being used, as frequency 
	values are not generated for all occurrence values.	That is, not all points on the 
	x-axis will be used when the graph is plotted. This causes problems when trying to use 
	data points as if they are spaced at unit length along the x-axis. To combat this problem, 
	this function returns a dict which contains frequency values for all x values (most of 
	which could well be 0). This allows the data to be used in the correct manner. 
	"""
	
	for i in xrange(sorted(hist_dict.keys())[0], sorted(hist_dict.keys())[-1]):
		hist_dict.setdefault(i,0)	
	return hist_dict


def compute_num_kmer_words(hist_dict):

	total_kmer_words = 0
	total_kmer_words = sum(occurrence * hist_dict[occurrence] for occurrence in hist_dict)
	
	return total_kmer_words


def generate_histogram(input_file_path, k_mer_size):
	
	"""
	Effectively ensures that a .hgram file exists and is stored at the correct location for
	the file stored at 'input_file_path'. 
	"""

	file_name = input_file_path.split("/")[-1].split(".")[0]
	extension = input_file_path.split("/")[-1].split(".")[-1]
	
	if os.path.isfile("/lustre/scratch110/sanger/gh10/Code/k_mer_scripts/hgram_data/" + 
	file_name + "_" + str(k_mer_size) + "mer.hgram"):
		return
	
	elif extension in ["data","dat"]:
		parse_data.parse(input_file_path, k_mer_size)
		
	elif extension in ["fasta","fa","fsa","fastq"]:
		if k_mer_size == []:
			raise Exception("Cannot use an empty k-mer size list when trying to create \
			histograms")
		compute_hist_from_fast(input_file_path, k_mer_size)
		
	else:
		raise Exception("Incorrect file extension.")
		

def compute_genome_size(hists_dict):

	genome_size_list = []
	for size in hists_dict.keys():
		hist_dict_augmented = format_data(hists_dict[size])		
		modes = calculate_modes(hist_dict_augmented, 1)
		
 		# genome_size = total num of k-mer words / first mode of occurences
		genome_size = compute_num_kmer_words(hists_dict[size]) / modes[0][0]
		genome_size_list.append((size, genome_size))
	
	return genome_size_list


def plot_graph(hists_dict):

	k_mer_sizes = hists_dict.keys()
	for size in k_mer_sizes:
		plt.plot(hists_dict[size].keys(), hists_dict[size].values())
	
	reload(graph_settings)
	settings = graph_settings.generate_settings() 
	
	plt.xlim(settings['x_lower'], settings['x_upper'])
	plt.ylim(settings['y_lower'], settings['y_upper'])
	plt.xscale(settings['x_scale'])
	plt.yscale(settings['y_scale'])
	plt.xlabel(settings['x_label'])
	plt.ylabel(settings['y_label'])
	
	plt.legend(hists_dict.keys())
	desired_title = raw_input("Enter graph's title: ")
	plt.title(desired_title)
	plt.tick_params(labelright = True)

	plt.show()
	
	return
	
	
def calculate_hist_dict(input_file_path, k_size):

	"""
	Returns dictionary consisting of keys corresponding to occurrences and values corresponding 
	to frequencies.
	"""
	
	generate_histogram(input_file_path, k_size)
		
	file_name = str(input_file_path.split("/")[-1].split(".")[0]) + "_" + str(k_size) + "mer" 
	extension = str(input_file_path.split("/")[-1].split(".")[-1])
	
	if extension == "hgram":
		hgram_name = input_file_path
	else:
		hgram_name = "/lustre/scratch110/sanger/gh10/Code/k_mer_scripts/hgram_data/" + \
		file_name + ".hgram"
		
	with open(hgram_name, 'r') as hgram_data:

		store_dict = {}
		for line in hgram_data.readlines():
			occ_and_freq = line.split(" ")
			store_dict[int(occ_and_freq[0])] = int(occ_and_freq[1])
				
	return store_dict

def generate_sample(hist_dict, sample_size):
	
	"""
	Generates and returns a sample of size 'sample_size' from 'hist_dict'. Each iteration 
	samples a single occurrence/frequency pair, and stores them in dict 'sample'. This dict 
	is then returned in the same format as a hist_dict.
	""" 
	
	sample = {}
	total_kmer_reads = compute_num_kmer_words(hist_dict)
	
	for i in xrange(sample_size):
		x = random.randint(1,total_kmer_reads)
		iCount = hist_dict[1]
		j = 1
		for j in hist_dict.keys():
			iCount += j * hist_dict[j]
			if iCount > x:
				break
			j += 1
		
		# Increment occurrence frequency count by 1, add occurrence to sample with value 1
		# if previously unobserved
		sample[j] = sample.get(j,0) + 1 
	
	return sample


def get_path_k_size():
	
	"""
	Prompts the user to input the location of the data they want to use, and a list containing
	all the k-mer siaes which they wish to use.
	"""
	
	input_file_path = raw_input("Enter the location of data to be used: ") 
	# Often stored at /lustre/scratch110/sanger/gh10/jellyfish/
	
	if input_file_path == "":
		raise Exception("File location cannot be empty!")
	
	extension = input_file_path.split("/")[-1].split(".")[-1]
	if extension != "hgram":
		k_mer_sizes = input("\nEnter k-mer sizes: ") # List containing all desired k-mer sizes
		if k_mer_sizes == []:
			raise Exception("k-mer size array cannot be empty!")
	else:
		# k-mer size of .hgram file
		k_mer_size = input_file_path.split("/")[-1].split(".")[0].split("_")[-1][:-3]
		k_mer_sizes = [k_mer_size]	
	
	# Dict in which to store k-mer size as key, and hist_dict for that k-mer size as value:
	hists_dict = {} 
	
	for size in k_mer_sizes:
		hists_dict[size] = calculate_hist_dict(input_file_path, size)
		
	return hists_dict

		
def main():

	hists_dict = get_path_k_size()

	while True:
	
		choice = raw_input("\n1 = Calculate Size; 2 = Plot Graph; 0 = Change Data; ENTER = Exit; ")
		
		if choice == "1":
			# Calculate genome size
			print ""
			for size in compute_genome_size(hists_dict):
				print "Size calculated to be " + str(size[1]) + " base pairs (using " + \
				str(size[0]) + "mers)"
		
		elif choice == "2":
			# Plot graphs
			plot_graph(hists_dict)	
			
		elif choice == "0":
			# Get new input
			hists_dict = get_path_k_size()		
			
		else:
			print "\nExiting...\n"
			return		
				
								
if __name__ == "__main__":
	main()

