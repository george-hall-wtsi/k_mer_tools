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


def generate_settings():	

	settings = {}

	# GRAPH SETTINGS: 
	
	settings['x_lower'] = 1
	settings['x_upper'] = 500
	settings['y_lower'] = 1
	settings['y_upper'] = 10**8
	
	# Scales can be linear or log
	settings['x_scale'] = 'linear'
	settings['y_scale'] = 'log'
	settings['x_label'] = 'k-mer Coverage'
	settings['y_label'] = 'k-mer Count Frequency'




	# REPEATS SETTINGS:

	settings['desired_border'] = 0.1 # Desire middle 80%


	return settings

