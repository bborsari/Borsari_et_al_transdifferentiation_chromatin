#!/bin/bash

#********
# USAGE *
#********

display_usage() { 
	echo -e "DESCRIPTION: provided an encode metadata file, returns a list of 'File accession' Ids, one for each experimental target\n"
	echo -e "\t--input <Input metadata file>\n"
	echo -e "\t--format <File format>\n"
	echo -e "\t--type <Output type>\n"
	echo -e "\t--assembly\n"
        echo -e	"\t--status <File status>\n"
	echo -e "\t--assay\n"
	echo -e "\t--cell_line <Biosample term name>\n"
	echo -e "\t--outFile (default: 'parsed.metadata.tsv')\n"
	echo -e "\t--verbose <yes/no> (default: no)"
} 


if [[  $1 == "--help" ||  $1 == "-h" ]]
then
    	display_usage
        exit 0
fi


if [  $# -le 13  ]
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
    	
	--format)
    	format="$2"
    	shift
    	;;
    	
	--type)
    	type="$2"
    	shift
    	;;
    	
	--assembly)
    	assembly="$2"
    	shift
    	;;

	--status)
	status="$2"
	shift
	;;

	--assay)
	assay="$2"
	shift
	;;

	--cell_line)
	cell_line="$2"
	shift
	;;

	--outFile)
	outFile="$2"
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


: ${outFile:="parsed.metadata.tsv"}
: ${verbose:="no"}


if [[ "$verbose" == "yes" ]]
then
	echo "Reading options .."
	echo "Input metadata file =" "${input}"
	echo "File format =" "${format}"
	echo "Output type =" "${type}"
	echo "Assembly =" "${assembly}"
	echo "File status =" "${status}"
	echo "Assay =" "${assay}"
	echo -e "Biosample term name =" "${cell_line}\n"
fi





#**************************
# PARSE AUDIT CATEGORIES *
#**************************


if [[ "$verbose" == "yes" ]]
then
	echo -e "Parsing audit categories ..\n"
fi

awk -v format="$format" -v type="$type" -v assembly="$assembly" -v status="$status" -v assay="$assay" -v cell_line="$cell_line" '
BEGIN{FS="\t"}
$2==format && $3==type && $5==assay && $7==cell_line && $44==assembly && $48==status
' "$input" | /no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/bin/parse.metadata.audit.categories.py  > "$outFile"




#**************
# HEADERIZING *
#**************

if [[ "$verbose" == "yes" ]]
then
	echo "Headerizing .."
fi

sed -i '1iTarget\tExperiment_id\tFile_id\tReplicate(s)' "$outFile"
