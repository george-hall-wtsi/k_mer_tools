Set of k-mer tools used primarily for quality control. Nowhere near finished. Will also only work on my local machine.

The algorithm used in parameter_estimation.py is based off Jared Simpson's paper 'Exploring Genome Characteristics and Sequence Quality Without a Reference'. The code used in the function 'learn_mixture_parameters' is largely based off code used in his implementation of the algorithm, which can be found at https://github.com/jts/sga/blob/master/src/SGA/preqc.cpp#L429.






Copyright (c) 2015 Genome Research Ltd. 

Author: George Hall gh10@sanger.ac.uk 

This program is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software 
Foundation; either version 3 of the License, or (at your option) any later 
version. 

This program is distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
details. 

You should have received a copy of the GNU General Public License along with 
this program. If not, see <http://www.gnu.org/licenses/>. 
