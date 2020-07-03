#!/usr/bin/env python

# takes in input the output of NW.alignment.path.R and, for each pair, it reports
# 'match' if it's the best bidirectional match involving the two time points, otherwise 'insertion'


#************
# LIBRARIES *
#************

import sys
from optparse import OptionParser
from decimal import *


#*****************
# OPTION PARSING *
#*****************

parser = OptionParser()
parser.add_option("-i", "--input", dest="input", default="stdin")
options, args = parser.parse_args()

open_input = sys.stdin if options.input == "stdin" else open(options.input)



#********
# BEGIN *
#********

i = 1


expression_cols = ['expression_000h','expression_003h','expression_006h','expression_009h','expression_012h','expression_018h','expression_024h','expression_036h','expression_048h','expression_072h','expression_120h','expression_168h']
mark_cols = ['mark_000h','mark_003h','mark_006h','mark_009h','mark_012h','mark_018h','mark_024h','mark_036h','mark_048h','mark_072h','mark_120h','mark_168h']



alignments = {}

for line in open_input.readlines():
	if (i > 1):
	 	id, step_type, score, expression_time_point, mark_time_point, expression_z, mark_z, distance = line.strip().split('\t')
		alignments[id] = alignments.get(id, {})
		float_score = Decimal(distance)

		if (	(alignments[id].get(expression_time_point) == None)
			or
			  	(float_score < alignments[id][expression_time_point][0])
				or
					(	(float_score == alignments[id][expression_time_point][0])
					  and	(mark_time_point < alignments[id][expression_time_point][1])
					)
		):
			alignments[id][expression_time_point] = [float_score, mark_time_point]
					
			
                if (    (alignments[id].get(mark_time_point) == None)
                       	or
                                (float_score < alignments[id][mark_time_point][0])
                               	or
                                       	(       (float_score == alignments[id][mark_time_point][0])
                                       	  and   (expression_time_point < alignments[id][mark_time_point][1])
                                       	)
               	):
                        alignments[id][mark_time_point] = [float_score, expression_time_point]
	i += 1

print '\t'.join(['gene_id', 'expression_time_point', 'mark_time_point', 'distance', 'step-type'])

for id in alignments.keys():
	lines = []
	for time_point in expression_cols:
		corresponding_time_point = alignments[id][time_point][1]
		if (alignments[id][corresponding_time_point][1] == time_point):
			type = "match"
		else:
			type = "insert"
		l = [id, time_point, corresponding_time_point, str(alignments[id][time_point][0]), type]
		lines.append(('\t'.join(l)))
	for time_point in mark_cols:
                corresponding_time_point = alignments[id][time_point][1]
                if (alignments[id][corresponding_time_point][1] == time_point):
                        type = "match"
                else:
                     	type = "insert"
                l = [id, corresponding_time_point, time_point, str(alignments[id][time_point][0]), type]
                if (('\t'.join(l)) not in lines):
			lines.append(('\t'.join(l)))
	lines.sort()
	for l in lines:
		print l



 
