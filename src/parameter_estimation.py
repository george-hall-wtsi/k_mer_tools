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


from __future__ import division																						

import sys
import time
import itertools 

import math
import matplotlib
import matplotlib.pyplot as plt
import scipy.stats

import main

def generate_initial_parameters(hist_dict, num_parameters):	

	"""
	Generates the initial lists of occurrence rates (lambda_list) and the list of mixture 
	proportions	(mixture_proportions) which are required when using EM to generate new
	parameters. They are generated according to the specification laid out in 'Exploring 
	Genome Characteristics and Sequence Quality Without a Reference'. 
	"""
	
	first_mode = main.calculate_modes(hist_dict, 1)[0][0]
	
	lambda_list = []
	if num_parameters == 1:
		print "Warning: Only one lambda value"
		return lambda_list
	elif num_parameters > 1: 
		for i in xrange(0, num_parameters):
			lambda_list.append(first_mode * (i + 1) + (40 * (i)))
	else:
		raise Exception("num_parameters must be >= 0")

	parameters = [lambda_list,[(1 / num_parameters) for i in xrange(0, num_parameters)]] 
	return parameters
	
	
def log_poisson(k,m):
	return (k * math.log(m) - m - math.log(math.factorial(k)))
	
	
def calculate_poisson_prob(lambda_vals, weights, c_val):
	
	"""Returns value of zero-truncated Poisson Distribution formula as given in paper. """
	
	return sum(((scipy.stats.poisson(lambda_vals[i]).pmf(c_val) / 						\
	(1 - math.exp(-lambda_vals[i]))) * weights[i]) for i in xrange(len(weights)))
	
	
def learn_mixture_parameters(original, k_size, num_parameters):
	
	"""
	Attempts to approximate k-mer count data using a mixture model of Poisson distributions. 
	Uses the probability formula as set out in Jared Simpson's paper 'Exploring Genome 
	Characteristics and Sequence Quality Without a Reference' and the Expectation Maximisation 
	algorithm. This function is largely based off Jared's code, which can be found at: 
	https://github.com/jts/sga/blob/master/src/SGA/preqc.cpp#L429
	"""
	
	lambdas, mixture_proportions = generate_initial_parameters(original, num_parameters)
	#sample_hist_dict = main.generate_sample(original, sample_size = 50000)
	
	sample_hist_dict = original
	
	max_copy_number = lambdas[2] * (len(lambdas) - 1)
	count_vector = {}
	for key, value in sorted(sample_hist_dict.iteritems()):
		if len(count_vector) < max_copy_number:
			count_vector[key] = value
		else:
			break

	error_sum = 0
	error_n = 0
	prev_ll = -float(sys.maxint) 
	
	max_iterations = 30
	iteration = 0
	
	number_of_components = len(lambdas)
	
	while iteration < max_iterations:
		print "Beginning iteration " + str(iteration)
		iteration += 1
		
		# Array storing number of k-mers which have been softly classified as each state for 
		# the current iteration
		mixture_counts = [0 for x in xrange(number_of_components)]
		
		for occurrence, frequency in count_vector.iteritems():
			
			# Calculate logs of P(c) and P(c|m_i)P(m_i) by summing over all mixture states
			mixture_log_p = [0 for x in xrange(number_of_components)]
			log_p_c = 0.0
			
			for i in xrange(num_parameters):	
				
				# Calculate the zero-truncation term:
				log_zero_trunc = math.log(1 / (1 - math.exp(-lambdas[i])))
				
				# P(c| m_i) * P(m_i)
				mixture_log_p[i] = log_poisson(occurrence, lambdas[i]) + log_zero_trunc + \
				math.log(mixture_proportions[i])
				
				if (i==0):
					log_p_c = mixture_log_p[i]
				else:
					# Implements addLogs inline function from JS's C++ program
					log_p_c = (lambda a,b: a + (math.log(1.0 + math.exp(b - a))) if (a > b) 
					else b + (math.log(1.0 + math.exp(a - b))))(log_p_c, mixture_log_p[i])
							
			# Calculate P(m_i|c): 
			mixture_log_p = [(probability - log_p_c) for probability in mixture_log_p]
			#mixture_counts = [mixture_counts[i] + (math.exp(mixture_log_p[i]) 
			#* frequency) for i in xrange(number_of_components)]
			for i in xrange(number_of_components):
				mixture_counts[i] += (math.exp(mixture_log_p[i]) * frequency)
	
			error_sum += occurrence * frequency * math.exp(mixture_log_p[0])
			error_n += frequency * math.exp(mixture_log_p[0])
				
		# Recalculate mixture proportions:
		sum_of_mixture_counts = sum(i for i in mixture_counts)
		for i in xrange(number_of_components):
			mixture_proportions[i] = mixture_counts[i] / sum_of_mixture_counts
			
		# Recalculate lambda_0:
		lambdas[0] = error_sum / error_n
		
		# Calculate log-likelihood:
		sum_ll = 0.0
		for occurrence, frequency in count_vector.iteritems():
			partial_ll = 0.0
			for i in xrange(number_of_components):
				# Calculate zero-truncation term
				log_zero_trunc = math.log(1 / (1 - math.exp(-lambdas[i])))
				ll = log_poisson(occurrence, lambdas[i]) + log_zero_trunc + \
				math.log(mixture_proportions[i])
				if (i == 0):
					partial_ll = ll
				else:
					# Implements addLogs inline function from JS's C++ program
					partial_ll = (lambda a,b: a + (math.log(1.0 + math.exp(b - a))) \
					if (a > b) else b + (math.log(1.0 + math.exp(a - b))))(partial_ll, ll)
					
			partial_ll = partial_ll * frequency
			sum_ll += partial_ll
			
		difference = sum_ll - prev_ll
		ll_improvement = difference / abs(sum_ll)
		if (ll_improvement < 0.00001):
			break
		
		prev_ll = sum_ll
		
	if (iteration >= max_iterations):
		print "Warning: Coverage model failed to converge"
		
	return [lambdas, mixture_proportions]
	

