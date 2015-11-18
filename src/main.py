#! /usr/bin/env python


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


import os.path

# Activate virtualenv
file_location = os.path.dirname(__file__)
activate_this = os.path.join(file_location, 'virtualenv/bin/activate_this.py')
execfile(activate_this, dict(__file__ = activate_this))
# Virtualenv should now be activated

import sys
import subprocess32
import random
import argparse
import math

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as spysig

import scripts.parse_dat_to_histo as parse_data
import settings.all_settings as all_settings


def simulate_reads(reference, coverage, read_length, insert_size, name):

	"""
	Cuts reads at random locations on a reference genome and stores them in a .fasta file.
	"""

	subprocess32.call(['sh', os.path.join(os.path.dirname(__file__), "scripts/sim_reads.sh"), 
		reference, str(coverage), str(read_length), str(insert_size), name, 
		os.path.dirname(__file__)])

	return


def update_assembly_config(new_location):
	config_location = os.path.join(os.path.dirname(__file__), "scripts/assembly_config")
	with open(config_location, 'r') as assembly_config:
		lines = assembly_config.readlines()
	lines[10] = new_location
	with open(config_location , 'w') as assembly_config:
		assembly_config.writelines(lines)

	return


def process_peak(file_path, file_name, lower_limit, upper_limit, peak_number, reference_path, 
	assembler, k_size, assembler_k, processors):

	"""
	Takes a file and computes k-mer words present in the section of the k-mer spectrum graph 
	between lower_limit and upper_limit. These k-mer words are then assembled into contigs. If 
	the reference sequence has been provided (i.e. the reads have been simulated from a 
	reference for error checking), these contigs are mapped against it.
	"""

	if assembler_k >= k_size:
		raise Exception("Assembler k-mer size must be smaller than overall k-mer size")

	subprocess32.call(['sh', os.path.join(os.path.dirname(__file__), 
		"scripts/compute_k_mer_words.sh"), file_name, str(lower_limit), str(upper_limit), 
		str(peak_number), os.path.dirname(__file__), str(k_size)]) 

	if assembler == 'soap':
		new_location = "q=" + os.path.abspath(file_path).split(".")[0] + "_reads/peak_" + \
			str(peak_number) + "_k_mers-read.fastq\n"
		update_assembly_config(new_location)

	subprocess32.call(['sh', os.path.join(os.path.dirname(__file__), 
		"scripts/assemble_repeats.sh"), os.path.abspath(file_path), str(peak_number), 
		os.path.dirname(__file__), assembler, str(assembler_k), str(processors)])
	
	if reference_path != "":
		subprocess32.call(['sh', os.path.join(os.path.dirname(__file__), 
			"scripts/align_sim_to_ref.sh"), os.path.abspath(reference_path), 
			os.path.abspath(file_path), str(peak_number), os.path.dirname(__file__), 
			str(assembler_k), str(processors)])

	return


def ranges_from_extrema(extrema):
	minima = extrema['Min']
	intervals = [(y - x) for (x, y) in zip([m for m in minima], [m for m in minima[1:]])] 

	peak_ranges = zip(minima, minima[1:])
	peak_widths = [(j - i) for (i, j) in peak_ranges]
	new_ranges = []
	for i in xrange(len(peak_ranges)):
		new_ranges.append([0, 0])
		reload(all_settings)
		settings = all_settings.generate_settings()
		desired_border = settings['desired_border']

		new_ranges[i][0] = peak_ranges[i][0] + (desired_border * peak_widths[i])
		new_ranges[i][1] = peak_ranges[i][1] - (desired_border * peak_widths[i])
	peak_ranges = [(int(i), int(j)) for (i, j) in new_ranges][1:]

	return peak_ranges


def calculate_peak_ranges(hist_dict, max_peak):
	extrema = find_extrema(hist_dict, max_peak)
	
	return ranges_from_extrema(extrema)


def find_repeats(hist_dict, file_path, max_peak, assembler, k_size, assembler_k, 
	processors, reference_path = ""):
	
	"""
	Finds distinct peaks of k-mer spectrum, then uses Smalt to discover k-mer words associated
	with each peak (i.e. which occur within an interval half the width of the peak either side
	of the peak. If the optinal reference sequence has been provided, it is shredded and mapped
	against itself, to discover sequence of length 500 or more which are repetitive. This is 
	used to test the de novo repetition detection. 
	"""
	
	file_path = os.path.abspath(file_path)
	if reference_path:
		reference_path = os.path.abspath(reference_path)
	src = os.path.dirname(__file__)
	working_dir = file_path.split(".")[0] + "_reads"

	file_name = file_path.split("/")[-1]

	peak_ranges = calculate_peak_ranges(hist_dict, max_peak)

	for (peak_number, (lower_limit, upper_limit)) in enumerate(peak_ranges, 2):
		print "Started processing peak" , peak_number
		process_peak(file_path, file_name, lower_limit, upper_limit, peak_number, 
			reference_path, assembler, k_size, assembler_k, processors)
		
		if reference_path != "":
			# Mask repeats found in each peak (replace their loci with Xs on a copy of 
			# the reference fasta)
			subprocess32.call(['sh', os.path.join(src, "scripts/mask_repeats.sh"), 
				reference_path, working_dir, src, 
				(working_dir + "/peak_" + str(peak_number) +"/peak_" + str(peak_number) + \
				"_map")])

		print "Finished processing peak number" , peak_number

	if reference_path != "":	
		# 'Shred' reference and map to itself (to find all repeats for testing purposes):
		update_assembly_config("q=" + reference_path + "\n")
		subprocess32.call(['sh', os.path.join(src, "scripts/ssaha_shred.sh"), 
			reference_path, file_name.split(".")[0], src])

		# Mask repeated regions from each mode in shredded reads
		subprocess32.call(['grep', ':00', working_dir + "/shred_map"], 
			stdout = open(working_dir + "/shred_grep", "w"))
		
		with open(working_dir + "/shred_grep", "r") as f:
			data = [line.split() for line in f.readlines()]

		for n in xrange(2, max_peak + 1):
			print "Masking repeats occuring " + str(n) + " times"
			iCount = 0
			for i  in xrange(1, len(data) - n):
				if all(x[2] == data[i][2] for x in data[i+1: i+n]) and \
					(data[i][2] != data[i-1][2]) and (data[i][2] != data[i+n][2]):
					
					for line in data[i:i+n]:
						with open(working_dir + "/shred_" + str(n) + "_repeats", "a") as out:
							out.write(" ".join(x for x in line[:3]) + " " + \
								" ".join(str(x).rjust(10) for x in line[3:8]) + " " + \
								" ".join(str(x) for x in line[8:] ) + "\n")
		
			subprocess32.call(['sh', os.path.join(src, "scripts/mask_repeats.sh"), 
			reference_path, working_dir, src, os.path.abspath(working_dir + "/shred_" + \
				str(n) + "_repeats")])

	return 


def calculate_ex_score(ex_dict):

	"""
	Returns a score describing how well the predicted extrema match with where we expect them 
	to occur. Basically assesses how periodic the extrema are, and returns lower (that is, 
	better) scores for those which display (non-trivially) periodic behaviour. 
	"""
	
	diff_list = []

	# Don't allow trivially 'periodic' extrema (normally very crowded around the origin):
	if ex_dict['Max'][0] < 3:
		return float("inf")

	# Check that max is within peak range	
	peak_ranges = ranges_from_extrema(ex_dict)
	for (mx, peak_range) in zip(ex_dict['Max'][1:], peak_ranges):
		if (peak_range[0] < mx < peak_range[1]) == False:
			return float("inf")

	for x in xrange(1, len(ex_dict['Max'])):
		diff_list.append(abs(ex_dict['Max'][x] - (ex_dict['Max'][0] * (x + 1))))

	score = sum((float(i) / ex_dict['Max'][0]) for i in diff_list) / len(diff_list)
	score = abs(score * (ex_dict['Max'][0] - (0.5 * ex_dict['Max'][1])))
	average_width = (sum((y - x) for (x, y) in zip(ex_dict['Min'], ex_dict['Min'][1:])) \
		/ float(len(ex_dict['Min']) - 1))
	score = score * abs((ex_dict['Min'][1] - ex_dict['Min'][0]) - average_width)

	return score


def estimate_extrema(hist_dict, window_size, order_num, num_peaks_desired):
 
	"""
	Smooths the data using a moving average. Uses Scipy.signals.argrelextrema to then detect
	which points correspond to extrema. 
	"""

	window = np.ones(int(window_size))/float(window_size)
	moving_average = np.convolve(hist_dict.values(), window, 'same')
	smoothed_data = dict(zip(hist_dict.keys(), [int(x) for x in moving_average]))

	store_dict = {'Min': [], 'Max': []}

	min_list = spysig.argrelextrema(np.array(smoothed_data.values()), np.less_equal, 
		order = order_num)[0].tolist()[1:]
	max_list = spysig.argrelextrema(np.array(smoothed_data.values()), np.greater_equal, 
		order = order_num)[0].tolist()[1:]

	store_dict['Min'].append(min_list[0])
	iCount = 0
	min_index = 1
	max_index = 0
	
	while iCount < num_peaks_desired:
		
		if max_index >= len(max_list) or min_index >= len(min_list):
			break
		
		while max_list[max_index] < store_dict['Min'][-1]:
			max_index += 1
			if max_index == len(max_list):
				break
		else:
			store_dict['Max'].append(max_list[max_index])
		while min_list[min_index] < store_dict['Max'][-1]:
			min_index += 1
			if min_index == len(max_list):
				break
		else:
			store_dict['Min'].append(min_list[min_index])
		iCount += 1

	return store_dict

		
def find_extrema(hist_dict, num_peaks_desired):
	
	"""
	Returns a dict with 2 keys (Max and Min) with the values for each of these keys being 
	tuples which correspond to (occurrence, frequency) pairs which are either a maximum or a 
	minimum.
	"""

	hist_dict = pad_data(hist_dict)
	score_list = []
	for i in xrange(max(2, num_peaks_desired), num_peaks_desired + 4):
		score_list.append([find_extrema_main(hist_dict, i), i]) 
	sorted_scores = sorted(score_list, key = lambda x: x[0][1])

	extrema = sorted_scores[0][0][0]

	if extrema['Max'] == [] or extrema['Min'] == []:
		raise Exception("Could not find suitable extrema")

	return {'Max': extrema['Max'][:num_peaks_desired], 'Min': extrema['Min'][:num_peaks_desired + 1]}


def find_extrema_main(hist_dict, num_peaks_desired):
	(window_size, order_num) = (10, 10)

	while True:
		score_list = []
		for (w, o) in [(window_size, order_num), (window_size + 5, order_num), 
			(window_size, order_num + 1), (window_size - 5, order_num), 
			(window_size, order_num - 1)]:

			if (w > 0) and (o > 0):
				score_list.append(((w, o), 
					calculate_ex_score(estimate_extrema(hist_dict, w, o, num_peaks_desired))))

		sort_scores = sorted(score_list, key = lambda x: x[1])

		if sort_scores[0][1] == float("inf"):
			if window_size > 1000:
				return [{'Max': [], 'Min': []}, float("inf")]
			(window_size, order_num) = (window_size + 10, order_num + 10)
			continue
		
		# Current estimate is the best we can do
		elif sort_scores[0][0] == (window_size, order_num):
			extrema = estimate_extrema(hist_dict, window_size, order_num, num_peaks_desired)
			return [extrema, sort_scores[0][1]]

		# Perfect score, so return
		elif sort_scores[0][1] == 0.0:
			extrema = estimate_extrema(hist_dict, sort_scores[0][0][0], sort_scores[0][0][1], 
				num_peaks_desired)
			return [extrema, sort_scores[0][1]]

		else:
			(window_size, order_num) = sort_scores[0][0]




		
def pad_data(hist_dict):

	"""
	This function is required when, for example, simulated data is being used, as frequency 
	values are not generated for all occurrence values.	That is, not all points on the 
	x-axis will be used when the graph is plotted. This causes problems when trying to use 
	data points as if they are spaced at unit length along the x-axis. To combat this problem, 
	this function returns a dict which contains frequency values for all x values (most of 
	which could well be 0). This allows the data to be used in the correct manner. 
	"""

	for i in xrange(1, sorted(hist_dict.keys())[-1]):
		hist_dict.setdefault(i,0)	

	return hist_dict


def compute_genome_size(hists_dict):

	genome_size_list = []
	for size in hists_dict.keys():
		# Calculate more than the first extremum in order to more accurately estimate the peaks
		mode = find_extrema(hists_dict[size], 3)['Max'][0]
 		# Genome Size = total num of k-mer words / first mode of occurences
		genome_size = compute_num_kmer_words(hists_dict[size]) / mode
		genome_size_list.append((size, genome_size))
	
	return genome_size_list


def plot_graph(hists_dict, graph_title, use_dots, max_peak = None):

	k_mer_sizes = hists_dict.keys()
	for size in k_mer_sizes:
		padded_data = pad_data(hists_dict[size])

		if use_dots:
			plt.plot(padded_data.keys(), padded_data.values(), 'o')
		else:
			plt.plot(padded_data.keys(), padded_data.values())

		if max_peak is not None:

			if max_peak <= 0:
				raise Exception("Maximum desired peak must be a positive integer")

			extrema = find_extrema(hists_dict[size], max_peak)
			for (extremum, ordinates) in extrema.items():
				if extremum == 'Max':
					peak_ranges = ranges_from_extrema(extrema)
					for (lower, upper) in peak_ranges:
						plt.axvspan(lower, upper, color = 'green', alpha = 0.5)
					for x in ordinates:
						plt.axvline(x, c = 'b')
				else:
					for x in ordinates:
						plt.axvline(x, c = 'r')
				
				
	reload(all_settings)
	settings = all_settings.generate_settings() 
	
	plt.xlim(settings['x_lower'], settings['x_upper'])
	plt.ylim(settings['y_lower'], settings['y_upper'])
	plt.xscale(settings['x_scale'])
	plt.yscale(settings['y_scale'])
	plt.xlabel(settings['x_label'])
	plt.ylabel(settings['y_label'])

	plt.title(graph_title)
	plt.legend(hists_dict.keys())
	plt.tick_params(labelright = True)

	plt.show()
	
	return


def compute_num_kmer_words(hist_dict):

	total_kmer_words = 0
	total_kmer_words = sum(occurrence * hist_dict[occurrence] for occurrence in hist_dict)
	
	return total_kmer_words
	

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


def compute_hist_from_fast(input_file_path, k_size, processors, hash_size):
	
	"""
	Uses Jellyfish to count k-mers of length k_size from input file. 
	"""

	if (processors == 1) and (hash_size == 1):
		print "Number of processors used and hash size have both been left at their default"+ \
			"\nvalues. This is not a problem, but may not be what you wanted."

	print "Computing histogram data for k = " + str(k_size) + " for first time"
	print "Reading data for k = " + str(k_size)

	file_name = input_file_path.split("/")[-1].split(".")[0]
	mer_count_file = file_name + "_mer_counts_" + str(k_size) + ".jf"
	current_dir = os.path.dirname(__file__)

	# Count occurences of k-mers of size "k_size" in input file  

	s = 100000000 

	if k_size > 21: 
		# Otherwise memory usage prediction fails and starts predicting that no memory is going
		# to be used by hash table
		l = math.ceil(math.log(s, 2))
		mem_used = ((2**(l - 33)) * ((2*k_size) - l + 7))

		while mem_used < hash_size:
			s = s * 1.1
			l = math.ceil(math.log(s, 2))
			mem_used = ((2**(l - 33)) * ((2*k_size) - l + 7))

			if mem_used <= 0.0:
				break

		if s != 100000000:
			s = s / 1.1
		
	hash_size = int(s)

	subprocess32.call([os.path.join(current_dir, "../bin/jellyfish"), 
		"count", ("-m " + str(k_size)), "-s " + str(hash_size), "-t " + str(processors), "-C", 
		input_file_path, '-o', mer_count_file])

	print "Processing histogram for k = " + str(k_size)
	
	file_name = str(input_file_path.split("/")[-1].split(".")[0]) + "_" + str(k_size) + "mer"
	
	with open(file_name + ".hgram","w") as out_file:
		# Computes histogram data and stores in "out_file"
		subprocess32.call([os.path.join(current_dir, "../bin/jellyfish"), 
			"histo", mer_count_file], stdout = out_file)
	
	print "Finished for k = " + str(k_size)
	

def generate_histogram(input_file_path, k_mer_size, processors, hash_size, force_jellyfish):
	
	"""
	Essentially ensures that a .hgram file exists and is stored at the correct location for
	the file stored at 'input_file_path'. 
	"""
	
	file_name = input_file_path.split("/")[-1].split(".")[0]
	extension = input_file_path.split("/")[-1].split(".")[-1]
	
	if os.path.isfile(file_name + "_" + str(k_mer_size) + "mer.hgram") and not force_jellyfish:
		return
	
	elif extension in ["data","dat"]:
		parse_data.parse(os.path.abspath(input_file_path), k_mer_size)
		
	elif extension == "hgram":
		if str(k_mer_size) != file_name[-len(str(k_mer_size)) - 3:-3]:
			raise Exception("Incompatible k-mer size and .hgram file. ")
		else:
			return
	
	else:
		compute_hist_from_fast(input_file_path, k_mer_size, processors, hash_size)

def calculate_hist_dict(input_file_path, k_size, processors, hash_size, force_jellyfish):

	"""
	Returns dictionary consisting of keys corresponding to occurrences and values 
	corresponding to frequencies.
	"""
	
	generate_histogram(input_file_path, k_size, processors, hash_size, force_jellyfish)
		
	file_name = str(input_file_path.split("/")[-1].split(".")[0]) + "_" + str(k_size) + "mer" 
	extension = str(input_file_path.split("/")[-1].split(".")[-1])
	
	if extension == "hgram":
		hgram_name = input_file_path
	else:
		hgram_name = file_name + ".hgram"
		
	with open(hgram_name, 'r') as hgram_data:

		store_dict = {}
		for line in hgram_data.readlines():
			occ_and_freq = line.split(" ")
			store_dict[int(occ_and_freq[0])] = int(occ_and_freq[1])
	
	return store_dict


def argument_parsing():
	
	"""
	Uses argparse module to create an argument parser. Its first argument is the function which 
	the user wishes to execute.  
	"""

	# Most basic parser - all it asks for is path to some data
	base_parser = argparse.ArgumentParser(add_help = False,
		description = "A tool for computing genomic characteristics using k-mers")

	base_parser.add_argument("path", type = str, help = "location at which the data is stored")

	# Adds options to choose numbers of CPUs used and size of Jellyfish's hash table (could be 
	# relevant in all fucntions other than simulate)
	basic_options = argparse.ArgumentParser(add_help = False, parents = [base_parser])
	basic_options.add_argument("-p", "--processors", 
		help = "maximum number of CPUs used (default: 1)", default = 1, type = int)
	basic_options.add_argument("-s", "--hash-size", 
		help = "maximum size of Jellyfish's hash table in memory (in Gb). Only relevant if \
			Jellyfish has to count k-mers (default: 1)", 
		default = 1, type = int)
	basic_options.add_argument("-f", "--force-jellyfish", help =  "force Jellyfish to be run on\
			new data even if k-mers already appear to have been counted", action = "store_true")
	
	# Actual parser which is used
	parser = argparse.ArgumentParser()

	# For functions which must only have one value for k
	single_k_required = argparse.ArgumentParser(add_help = False, parents = [basic_options])
	single_k_required.add_argument("k", help = "k value to use", type = int, nargs = 1)

	# For functions which are capable of being passed multiple values of k
	multiple_k_possible = argparse.ArgumentParser(add_help = False, parents = [basic_options])
	multiple_k_possible.add_argument("k", help = "k value(s) to use (seperate with spaces)", 
		type = int, nargs = '+')

	# For functions which are in some way related to finding repetitive sequence
	some_repeats = argparse.ArgumentParser(add_help = False, parents = [single_k_required])
	some_repeats.add_argument("-r", "--reference", 
		help = "location of reference if reads are simulated", type = str, default = "")
	some_repeats.add_argument("-a", "--assembler", 
		help = "If SOAPdenovo is to be used, instead of SPAdes", type = str, 
		default = "spades", choices = ["soap", "spades"])
	some_repeats.add_argument("-d", "--assembler_k",  
		help = "k-mer size for assembler (must be smaller than overall k-mer size)",
		type = int, default = 31)

	subparsers = parser.add_subparsers(help = "select which function to execute")

	plot_subparser = subparsers.add_parser("plot", help = "plot k-mer spectra", 
		parents = [multiple_k_possible])
	size_subparser = subparsers.add_parser("size", help = "estimate genome size", 
		parents = [multiple_k_possible])
	repeats_subparser = subparsers.add_parser("repeats", 
		help = "find repetitive k-mer words, and align repetitive contigs to reference", 
		parents = [some_repeats])
	indiv_repeats_subparser = subparsers.add_parser("indiv-repeats", 
		help = "find repetitive k-mer words, and align repetitive contigs to reference for a \
		specified range", parents = [some_repeats])
	simulate_subparser = subparsers.add_parser("simulate", 
		help = "simulate random reads from reference genome", parents = [base_parser])

	plot_subparser.add_argument("-o", "--dots", 
		help = "plot the histogram using red dots", action = "store_true")
	plot_subparser.add_argument("-l", "--lines", 
		help = "draw lines to split graph into peaks up to peak number l", type = int)
	plot_subparser.add_argument("-t", "--title", help = "specify the title for the graph", 
		type = str, default = "")
	plot_subparser.set_defaults(func = "plot")

	size_subparser.set_defaults(func = "size")

	repeats_subparser.add_argument("max_peak", 
		help = "highest peak number to consider", type = int)
	repeats_subparser.set_defaults(func = "repeats")

	indiv_repeats_subparser.add_argument("peak_name", type = str, 
		help = "name of peak to be calulated")
	indiv_repeats_subparser.add_argument("l_lim", type = int, help = "lower limit of range")
	indiv_repeats_subparser.add_argument("u_lim", type = int, help = "upper limit of range")
	indiv_repeats_subparser.set_defaults(func = "indiv-repeats")

	simulate_subparser.add_argument("coverage", type = float, 
		help = "desired simulated coverage")
	simulate_subparser.add_argument("length", type = int, 
		help = "desired simulated read length")
	simulate_subparser.add_argument("insert", type = int, 
		help = "desired simulated insert size")
	simulate_subparser.add_argument("-n", "--name", type = str, 
		help = "name for simulated data", default = "")
	simulate_subparser.set_defaults(func = "simulate_reads")

	args = parser.parse_args()

	return args

		
def main():

	args = argument_parsing()

	# Calculate hists_dict if k_mer_words have been supplied
	if args.func in ['plot', 'size', 'repeats', 'indiv-repeats']:

		# Dict in which to store k-mer size as key, and hist_dict for that k-mer size as value:
		hists_dict = {}

		for size in args.k:
			hists_dict[size] = calculate_hist_dict(args.path, size, args.processors, 
				args.hash_size, args.force_jellyfish)

	if args.func == "plot":
		graph_title = args.title or args.path # If user has entered title then set title
		plot_graph(hists_dict, graph_title, args.dots, args.lines)

	if args.func == "size":
		for size in compute_genome_size(hists_dict):
			print "Size calculated to be " + str(size[1]) + " base pairs (using " + \
				str(size[0]) + "mers)"

	if args.func == "repeats":

		extension = ".".join(args.path.split("/")[-1].split(".")[1:])

		if extension not in ["fasta", "fastq"]:
			raise Exception("Incorrect file extension: file must be either .fasta or .fastq")

		for size in hists_dict.keys():
			file_name = args.path.split("/")[-1].split(".")[0]
			if not os.path.isfile(file_name + "_mer_counts_" + str(size) + ".jf"):
				compute_hist_from_fast(args.path, size, args.processors, args.hash_size)
			find_repeats(hists_dict[size], args.path, args.max_peak, args.assembler, size, 
				args.assembler_k, args.processors, args.reference)
			print "Finished finding repeats"

	if args.func == "indiv-repeats":
		extension = ".".join(args.path.split("/")[-1].split(".")[1:])

		if extension not in ["fasta", "fastq"]:
			raise Exception("Incorrect file extension: file must be either .fasta or .fastq")

		for size in hists_dict.keys():
			file_name = args.path.split("/")[-1].split(".")[0]
			if not os.path.isfile(file_name + "_mer_counts_" + str(size) + ".jf"):
				compute_hist_from_fast(args.path, size, args.processors, args.hash_size)

			process_peak(args.path, file_name, args.l_lim, args.u_lim, args.peak_name, 
				args.reference, args.assembler, size, args.assembler_k, args.processors)
			print "Finished finding repeats"

	if args.func == "simulate_reads":
		name = (args.name or args.path.split("/")[-1].split(".")[0])
		simulate_reads(args.path, args.coverage, args.length, args.insert, name)

	return


if __name__ == "__main__":
	main()

