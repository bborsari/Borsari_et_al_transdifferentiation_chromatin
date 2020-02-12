#!/bin/bash

#adds a header to a file; the header is R-type (length of the header = NF-1)

[ $# -ge 1 -a -f "$1" ] && input="$1" || input="-"

awk 'BEGIN{FS=OFS="\t"}{if(NR==1){
				for(i=1;i<=(NF-2);i++)
				{printf "V"i"\t"};
				printf "V"NF-1"\n";
				print $0}
			else
				{print $0}}' $input
