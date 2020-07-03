#!/usr/bin/env python


#************
# LIBRARIES *
#************

import sys
from optparse import OptionParser


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

d={}
vals=[]

for line in open_input.readlines():
	name, val = line.strip().split('\t')
	d[name] = d.get(name, {})
	d[name][val] = d[name].get(val, 0) +1
	vals.append(val)

vals = set(vals)
 
for key in d.keys():
	mylist = []
	for val in vals:
		if (val not in d[key].keys()):
			d[key][val] = 0
		mystring = ":".join([val, str(d[key][val])])
		mylist.append(mystring)
	print key, "\t", ";".join(mylist)
