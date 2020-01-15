#!/usr/bin/env python

#takes in input chromosome name, length of the chromosome and window length, and returns a file with genomic windows

#usage is: create_genomic_windows.py chr_name chr_length window

import sys

chr_name = sys.argv[1] 
chr_length = int(sys.argv[2])
window = int(sys.argv[3])

end = window
start = 1
while (end <= chr_length):
	print chr_name,"\t",start,"\t",end
	start = (end +1)
	end = (start + window -1)
