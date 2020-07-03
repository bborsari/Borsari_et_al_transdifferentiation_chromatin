#!/bin/bash


#********
# USAGE *
#********

display_usage() { 
	echo -e "DESCRIPTION: parses a list of experimental targets and corresponding bed/bigBed identifier(s)\n"
	echo -e "\t--input <.txt file with target_id, experiment_id, file_id, replicate (optional)>\n"
	echo -e "\t--dir <directory of bed/bigBed files>\n"
	echo -e "\t--ext <bed/bigBed> (file extension of peaks; default: 'bigBed')\n"
	echo -e "\t--verbose <yes/no> (default: no)\n"
} 


if [[  $1 == "--help" ||  $1 == "-h" ]]
then
    	display_usage
        exit 0
fi


if [  $# -le 3  ]
then
	echo "ERROR: insufficient number of arguments"
    	display_usage
        exit 1
fi


#******************
# READING OPTIONS *
#******************


while [[ $# -gt 1 ]]; do

	key="$1"
	
	case $key in

	--input)
	input="$2"
	shift
	;;

	--dir)
	dir="$2"
	shift
	;;
    	
	--ext)
	ext="$2"
	shift
	;;
	
	--verbose)
	verbose="$2"
	;;
	*)
	
	;;
	esac
	shift
done


: ${ext:="bigBed"}


#******************************
# PARSING LIST OF IDENTIFIERS *
#******************************

if [[ "$verbose" == "yes" ]]
then
	echo -e "Parsing list of identifiers..\n"
fi

tail -n+2 $input | sed 's/\-human//g' |  awk -v ext="$ext" -v dir="$dir" '
BEGIN {
	FS="\t"
	} 
{ 
if (NR==1) {
	target=$1
	i=1
	exp_id=$2
	b[i]=dir"/"$3"."ext 
	}
else {
	if ($2==exp_id) {
		i++; 
		b[i]=dir"/"$3"."ext 
		}  
	else {
		printf "%s\t", target 
		for(j=1; j<i; j++) {
			printf "%s,", b[j] 
			} 
		print b[i] 
		delete b
		i=1
		target=$1
		exp_id=$2
		b[i]=dir"/"$3"."ext
		}
	}
}
END {
	printf "%s\t", target 
	for(j=1; j<i; j++){ 
		printf "%s,", b[j] 
		} 
	print b[i] 
	}'
