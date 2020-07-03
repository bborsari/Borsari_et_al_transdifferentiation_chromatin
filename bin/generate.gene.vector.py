#!/usr/bin/env python

# takes in input the output of NW.bidirectional.matches.py and returns for each gene a 
# binary vector of size = 552 (46*12) or an integer vector of size = 12


#************
# LIBRARIES *
#************

import sys
from optparse import OptionParser



#******************
# OPTION PARSING *
#******************

parser = OptionParser()
parser.add_option("-i", "--input", dest="input", default="stdin")
parser.add_option("-t", "--type", dest="type", default="binary")
options, args = parser.parse_args()

open_input = sys.stdin if options.input == "stdin" else open(options.input)

# type can be either 'binary' or 'integer'


#********
# BEGIN *
#********


d = { 'expression_000h': 0,
	'expression_003h': -1,
	'expression_006h': -2,
	'expression_009h': -3,
	'expression_012h': -4,
	'expression_018h': -5,
	'expression_024h': -6,
	'expression_036h': -7,
	'expression_048h': -8,
	'expression_072h': -9,
	'expression_120h': -10,
	'expression_168h': -11,
	'mark_000h': 0,
	'mark_003h': 1,
	'mark_006h': 2,
	'mark_009h': 3,
	'mark_012h': 4,
	'mark_018h': 5,
	'mark_024h': 6,
	'mark_036h': 7,
	'mark_048h': 8,
	'mark_072h': 9,
	'mark_120h': 10,
	'mark_168h': 11 }

genes = {}

expression_cols = ['expression_000h',
	'expression_003h',
	'expression_006h',
	'expression_009h',
	'expression_012h',
	'expression_018h',
	'expression_024h',
	'expression_036h',
	'expression_048h',
	'expression_072h',
	'expression_120h',
	'expression_168h']

positions = { -11: 0,
	-10: 1,
	-9: 2,
	-8: 3,
	-7: 4,
	-6: 5,
	-5: 6,
	-4: 7,
	-3: 8,
	-2: 9,
	-1: 10,
	0: 11,
	1: 12,
	2: 13,
	3: 14,
	4: 15,
	5: 16,
	6: 17,
	7: 18,
	8: 19,
	9: 20,
	10: 21,
	11: 22 }



i=1



for line in open_input.readlines():
	if (i > 1):
	 	id, expression_time_point, mark_time_point, distance, step_type = line.strip().split('\t')

		genes[id] = genes.get(id, {})
		if (	(genes[id].get(expression_time_point) == None)
			or
				(	(genes[id].get(expression_time_point)[1] == 'insert')
					 and	( step_type == 'match' )
				)
		):
			genes[id][expression_time_point] = [mark_time_point, step_type]
	i += 1


for id in genes.keys():
	total_matrix = []
	for expression_time_point in expression_cols:
		mark_time_point = genes[id][expression_time_point][0]
		n = d[expression_time_point] + d[mark_time_point] 
		
		if ( type == 'binary' ):
			current_matrix = [[0]*23, [0]*23]
			pos = positions[n]
			if (genes[id][expression_time_point][1] == 'match'):
				current_matrix[0][pos] = 1
			else:
				current_matrix[1][pos] = 1 	

			current_matrix = sum(current_matrix, [])
			current_matrix = [str(i) for i in current_matrix]
			total_matrix.append(('\t'.join(current_matrix)))

		else:
			total_matrix.append((str(n)))			
	
	total_matrix = '\t'.join(total_matrix)
	print '\t'.join((id, total_matrix))


