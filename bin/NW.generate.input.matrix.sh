#!/bin/bash


#********
# USAGE *
#********

display_usage() { 
	echo -e "DESCRIPTION: provided matrices of expression and chromatin profiles, pastes the corresponding z-score matrices one next to the other\n"
	echo -e "\t--expression <matrix of expression profiles>\n"
	echo -e "\t--mark <matrix of chromatin profiles>\n"
	echo -e "\t--outFolder <output folder> (default: cwd)\n"
	echo -e "\t--outFile <output file> (default: 'NW.input.matrix.tsv')\n"
	echo -e "\t--keep <yes/no> (default: no)\n"
} 


if [[  $1 == "--help" ||  $1 == "-h" ]]
then
    	display_usage
        exit 0
fi


if [  $# -le 1  ]
then
	echo -e "ERROR: insufficient number of arguments\n"
    	display_usage
        exit 1
fi



#******************
# READING OPTIONS *
#******************

while [[ $# -gt 1 ]]; do

	key="$1"
	
	case $key in
    	
	--expression)
    	expression="$2"
    	shift # past argument
    	;;
    	
	--mark)
    	mark="$2"
    	shift # past argument
    	;;
    	    	
	--outFolder)
    	outFolder="$2"
    	shift # past argument
    	;;

	--outFile)
	outFile="$2"
	shift # past argument
	;;

	--keep)
	keep="$2"
	;;
	*)
	
	;;
	esac
	shift
done


: ${outFolder:="."}
: ${outFile:="NW.input.matrix.tsv"}
: ${keep:="no"}





#********
# BEGIN *
#********



echo -e "expression_000h\texpression_003h\texpression_006h\texpression_009h\texpression_012h\texpression_018h\texpression_024h\texpression_036h\texpression_048h\texpression_072h\texpression_120h\texpression_168h\tmark_000h\tmark_003h\tmark_006h\tmark_009h\tmark_012h\tmark_018h\tmark_024h\tmark_036h\tmark_048h\tmark_072h\tmark_120h\tmark_168h" > "$outFolder"/"$outFile"

~/software/R-3.5.2/bin/Rscript ~abreschi/R/scripts/scale_matrix.R -i $expression -C -S -n 0 -r > "$outFolder"/tmp_expression.tsv


~/software/R-3.5.2/bin/Rscript ~abreschi/R/scripts/scale_matrix.R -i $mark -C -S -n 0 -r > "$outFolder"/tmp_mark.tsv


~/bin/join.py -b <(tail -n+2 "$outFolder"/tmp_expression.tsv) -a <(tail -n+2 "$outFolder"/tmp_mark.tsv) | grep -v "NA" >> "$outFolder"/"$outFile"


if [[ "$keep" == "no" ]]
then
	rm "$outFolder"/tmp_expression.tsv
	rm "$outFolder"/tmp_mark.tsv
fi

