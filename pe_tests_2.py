#!/nfs/users/nfs_g/gh10/Documents/Code/Python/virtualenvs/venv/bin/python


################################################################################
# Copyright (c) 2015 Genome Research Ltd. 
# 
# Author: George Hall gh10@sanger.ac.uk 
# 
# This program is free software: you can redistribute it and/or modify it under 
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


import parameter_estimation
import plotting
import matplotlib
import matplotlib.pyplot as plt

def tests():

	file_path = "/lustre/scratch110/sanger/gh10/Data/yeast-miseq_1_0002_uncorrected_untrimmed.fasta"
	
	num_parameters = 6
	k_size = 31
	
	hist_dict = plotting.calculate_hist_dict(file_path, k_size)
	new_hist_dict = {}
	for key in hist_dict:
		if 40 <= key <= 400:
			new_hist_dict[key-40] = hist_dict[key]
			
			
	hist_dict = new_hist_dict
	
	lambdas, mixture_proportions = parameter_estimation.learn_mixture_parameters(hist_dict, k_size, num_parameters)
	
	print "Lambda values:" , lambdas
	print "Mixture proportions:" , mixture_proportions

	total_num_kmers = plotting.compute_num_kmer_words(hist_dict)
	total_kmers_without_singles = sum(i * hist_dict[i] for i in hist_dict.keys()) - hist_dict[1]

	first_mode = plotting.calculate_modes(hist_dict, 1)[0][0]

	print "Not corrected GS:" , int(total_num_kmers / first_mode) , "base pairs"
	print "Disregarding single occurrences:" , int(total_kmers_without_singles / first_mode) , \
	"base pairs"
	
	g_size = int(total_num_kmers / first_mode)


	# GRAPH TESTS:

	number_iterations = 1000
	y_list = []
		
	for j in xrange(1, number_iterations):
		y_list.append(g_size * parameter_estimation.calculate_poisson_prob(lambdas, 
		mixture_proportions, j))
	plt.plot(xrange(1, number_iterations), [y for y in y_list])
	
	for curve in xrange(len(lambdas)):
		print curve
		y_list = []
		for j in xrange(1, number_iterations):
			y_list.append(g_size * parameter_estimation.calculate_poisson_prob([lambdas[curve]], 
			[mixture_proportions[curve]], j))
		plt.plot(xrange(1, number_iterations), [y for y in y_list])

	new_hists = plotting.calculate_hist_dict(file_path, k_size)
	
	new_hist_dict = {}
	for key in new_hists:
		if 40 <= key <= 10000:
			new_hist_dict[key-40] = new_hists[key]
	new_hists = new_hist_dict
	plotting.plot_graph({k_size: new_hists})
	
	
	
if __name__ == "__main__":
	tests()
