#!/usr/bin/env python

# parses a metadata file from Encode based on the number of audit categories
# discards experiments with at least one audit error notification
# for each experimental target selects the experiment with lowest number of audit categories
# in case two experiments have equal number of audit categories,
# it selects the most recent one


#************
# LIBRARIES *
#************

import sys
from optparse import OptionParser


#*****************
# OPTION PARSING *
#*****************


parser = OptionParser()
parser.add_option("-I", "--input", dest="input", default="stdin")
options, args = parser.parse_args()

open_input = sys.stdin if options.input == "stdin" else open(options.input)



#********
# BEGIN *
#********


d1={}
d2={}

for line in open_input.readlines():
	splitline = line.strip("\r\n").split("\t")
	file_id = splitline[0]
	experiment_id = splitline[3]
        target = splitline[18]
        date = splitline[25]
	replicate = splitline[30]

	if splitline[49] == "":
		audit_warning = 0
	else:
		audit_warning = len(splitline[49].split(","))
		
	if splitline[50] == "":
		audit_internal_action = 0
	else:
		audit_internal_action = len(splitline[50].split(","))

	if splitline[51] == "":
		audit_not_compliant = 0
	else:
		audit_not_compliant = len(splitline[51].split(","))

	if splitline[52] == "":
		audit_error = 0
	else:
		audit_error = len(splitline[52].split(","))

	total_audit = audit_warning + audit_internal_action + audit_not_compliant
	if (audit_error == 0):
		d1[target] = d1.get(target, {})
		d1[target][experiment_id] = [total_audit, date]
		d2[experiment_id] = d2.get(experiment_id, {})
		d2[experiment_id][file_id] = replicate
		

for el in d1.keys():
	for i in range(0,len(d1[el].keys())):
		if (i == 0):
			min = int(d1[el][d1[el].keys()[i]][0])
			pos = 0
		else:	
			if (int(d1[el][d1[el].keys()[i]][0]) < min):
				min = d1[el][d1[el].keys()[i]][0]
				pos = i
			elif (d1[el][d1[el].keys()[i]][0] == min):
				
				date1 = d1[el][d1[el].keys()[pos]][1].split("-")
				date2 = d1[el][d1[el].keys()[i]][1].split("-")
				
				# check for year of release
				if (int(date1[0]) != int(date2[0])):
					if (max(int(date1[0]), int(date2[0])) == int(date1[0])):
						min = d1[el][d1[el].keys()[pos]][0]
						pos = pos
					else:
						min = d1[el][d1[el].keys()[i]][0]
						pos = i
				else:
					# check for month of release
					if (int(date1[1]) != int(date2[1])):
						if (max(int(date1[1]), int(date2[1])) == int(date1[1])):
							min = d1[el][d1[el].keys()[pos]][0]
							pos = pos
						else:
							min = d1[el][d1[el].keys()[i]][0]
							pos = i
					else:
						# check for day of release
						if (int(date1[2]) != int(date2[2])):
							if (max(int(date1[2]), int(date2[2])) == int(date1[2])):
								min = d1[el][d1[el].keys()[pos]][0]
								pos = pos
							else:
								min = d1[el][d1[el].keys()[i]][0]
								pos = i
						else:
							min = d1[el][d1[el].keys()[pos]][0]
							pos = pos
	
	for j in d2[d1[el].keys()[pos]]:
		l = []
		rep = (d2[d1[el].keys()[pos]][j]).replace(" ", "").replace("_", "")
		l.append(el)
		l.append(d1[el].keys()[pos])
		l.append(j) 
		l.append(rep)
		print "\t".join(l)

		


