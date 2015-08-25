#!/nfs/users/nfs_g/gh10/Documents/Code/Python/virtualenvs/venv/bin/python

import os.path
import subprocess32
import random

import matplotlib
import matplotlib.pyplot as plt

import parse_dat_to_histo as parse_data


def compute_hist_from_fasta(input_file_path, k_size):

	print "Computing hgram data for k = " + str(k_size) + " for first time"

	print "Reading data for k = " + str(k_size)
	file_location = input_file_path
	# Counts occurences of k-mers of size "k-size" in "file_input":  
	subprocess32.call(["/nfs/users/nfs_g/gh10/src/jellyfish-2.2.3/bin/jellyfish", "count", 
	("-m " + str(k_size)), "-s 100M", "-t 20", "-C", file_location])

	print "Processing histogram for k = " + str(k_size)
	
	file_name = str(input_file_path.split("/")[-1].split(".")[0]) + "_" + str(k_size) + "mer"
	
	with open("/lustre/scratch110/sanger/gh10/Code/k_mer_scripts/hgram_data/" + 
	file_name + ".hgram","w") as out_file:
		# Computes histogram data and stores in "out_file"
		subprocess32.call(["/nfs/users/nfs_g/gh10/src/jellyfish-2.2.3/bin/jellyfish", "histo", 
		"mer_counts.jf"], stdout=out_file)
	
	print "Finished for k = " + str(k_size)

		
def find_extrema(dct, window_size, alternate = False):
	
	"""
	Returns a dict with 2 keys (max, min) with the values corresponding to these keys 
	denoting occurrence, frequency pairs.
	"""
	
	i,j,k = (window_size), (window_size+1), (window_size+2)
	store_dict = {'Max' : [], 'Min' : []}

	if alternate == False:

		while k < (dct.keys()[-1] - window_size + 1):
			
			if (all(dct[i - x] < dct[j] > dct[k + x] for x in xrange(0, (window_size + 1)))):
				store_dict['Max'].append((j + 1, dct[j]))
			
			if (all(dct[i - x] > dct[j] < dct[k + x] for x in xrange(0, (window_size + 1)))):
				store_dict['Min'].append((j + 1, dct[j]))

			i += 1
			j += 1
			k += 1
			
	else:	
		prev_min, prev_max = False, False
		
		while k < (dct.keys()[-1] - window_size + 1):
			if prev_max == False and (all(dct[i - x] < dct[j] > dct[k + x] for x in xrange(
			0,(window_size + 1)))):
			
				store_dict['Max'].append((j + 1,dct[j]))
				prev_min, prev_max = False, True
			
			elif prev_min == False and (all(dct[i - x] > dct[j] < dct[k + x] for x in xrange(
			0,(window_size + 1)))):
			
				store_dict['Min'].append((j + 1,dct[j]))
				prev_min, prev_max = True, False
		
			i += 1
			j += 1
			k += 1
		
	return store_dict
	
	
def calculate_modes(hist_dict, n):
	"""Takes a hist_dict as input and returns its first n modes. """
	hist_dict_augmented = format_data(hist_dict)
	modes = []
	for x in reversed(sorted(find_extrema(hist_dict_augmented, window_size = 15, 
	alternate = False)['Max'], key = lambda pair: pair[1])):
		if len(modes) < n:
			modes.append(x)
		else:
			break

	return modes
	
	
def format_data(dct):

	"""
	This function is required for example when simulated data is being used, as frequency 
	values are not generated for all occurrence values.	That is, not all points on the 
	x-axis are used. This causes problems when trying to use data points as if they are 
	constantly spaced along the x-axis. Therefore this function returns a dict which 
	contains frequency values for all x values (most of which could well be 0). This 
	allows the data to be used in the correct manner. 
	"""
	
	for i in xrange(sorted(dct.keys())[0], sorted(dct.keys())[-1]):
		dct.setdefault(i,0)	
	return dct


def compute_num_kmer_words(hist_dict):

	total_kmer_words = 0
	total_kmer_words = sum(occurrence * hist_dict[occurrence] for occurrence in hist_dict)
	
	return total_kmer_words


def generate_histogram(input_file_path, k_mer_size): 

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
		compute_hist_from_fasta(input_file_path, k_mer_size)
		
	elif extension != "hgram":
		raise Exception("Incorrect file extension.")
		

def compute_genome_size(input_file_path, k_size_list):

	genome_size_list = []
	
	for k_size in k_size_list:
	
		hist_dict = calculate_hist_dict(input_file_path, k_size)
		hist_dict_augmented = format_data(hist_dict)
		
		modes = calculate_modes(hist_dict_augmented, 1)
 		# genome_size = total num of k-mer words / first mode of occurences
 		print compute_num_kmer_words(hist_dict) , modes[0][0]
		genome_size = compute_num_kmer_words(hist_dict) / modes[0][0]
		genome_size_list.append((k_size, genome_size))
	
	return genome_size_list


def plot_graph(input_file_path, k_mer_sizes):

	for k_mer in k_mer_sizes:
		hist_dict = calculate_hist_dict(input_file_path, k_mer)	
		plt.plot(hist_dict.keys(),hist_dict.values())
		
	plt.xlim(0,5000)
	#plt.ylim(0,2400)
	plt.yscale('log')
	plt.xlabel("Occurrences")
	plt.ylabel("Frequencies")
	plt.legend(k_mer_sizes)
	plt.title(input_file_path)
	plt.tick_params(labelright=True)

	plt.show()
	
	
def calculate_hist_dict(input_file_path, k_size):

	"""
	Returns dictionary consisting of keys corresponding to occurrences and values corresponding 
	to frequency.
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
	is then returned in the same format as hist_dict.
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
		
	return (input_file_path, k_mer_sizes)

		
def main():

	input_file_path, k_mer_sizes = get_path_k_size()

	while True:
	
		choice = raw_input("\n1 = Calculate Size; 2 = Plot Graph; 0 = Change Data; ENTER = Exit; ")
		
		if choice == "1":
			print ""
			for size in compute_genome_size(input_file_path, k_mer_sizes):
				print "Size calculated to be " + str(size[1]) + " base pairs"
		
		elif choice == "2":
			plot_graph(input_file_path, k_mer_sizes)	
			
		elif choice == "0":
			input_file_path, k_mer_sizes = get_path_k_size()
			
		else:
			print "\nExiting...\n"
			return		
				
								
if __name__ == "__main__":
	main()

