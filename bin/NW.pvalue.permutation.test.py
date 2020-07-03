#!/usr/bin/env python

# takes in input the output of NW.alignment.permutations.R and computes the p-value for each gene


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
options, args = parser.parse_args()

open_input = sys.stdin if options.input == "stdin" else open(options.input)



#********
# BEGIN *
#********

i=1


print '\t'.join(('score', 'p-value'))


for line in open_input.readlines():
	if (i > 1):
		pvalue = None
	 	id, score, permutations = line.strip().split('\t')
		p = permutations.split(',')
		counts=1
		for val in p:
			if (float(val) <= float(score)):
				counts+=1
			else:
				pvalue = str((float(counts)/(len(p)+1)))
				break
		if (pvalue is None):
			pvalue = str((len(p)+1)/(len(p)+1))
		print '\t'.join([id, score, pvalue])
	i+=1


 
