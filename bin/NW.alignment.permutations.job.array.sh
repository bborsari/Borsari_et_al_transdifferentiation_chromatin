#!/bin/bash
#$ -N NW.alignment.permutations.job.array
#$ -cwd
#$ -q rg-el7
#$ -j y
#$ -l h_rt=200:00:00
#$ -M beatrice.borsari@crg.eu -m a
#$ -o ../../logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -e ../../errors/$JOB_NAME.$JOB_ID.$TASK_ID.error
#$ -t 1-162


#********
# USAGE *
#********

display_usage() { 
	echo -e "DESCRIPTION: provided an input file and a sample replicate, runs a job array over 'NW.alignment.permutations.R' with default options (-d, -p, -b, --t1, --t2 flags)\n"
	echo -e "\t--input <NW input matrix (see '/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/bin/NW.generate.input.matrix.sh')>\n"
	echo -e "\t--rep <sample replicate> (default: leave it empty)\n"
	echo -e "\t--outFolder <output folder> (default: 'out.NW')\n"
} 


if [[  $1 == "--help" ||  $1 == "-h" ]]
then
    	display_usage
        exit 0
fi


if [  $# -lt 1  ]
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
    	
	--input)
    	input="$2"
    	shift
    	;;

	--outFolder)
	outFolder="$2"
	shift
	;;
    	
	--rep)
    	rep="$2"
    	;;
	*)    	    	
	
	;;
	esac
	shift
done


: ${rep:=""}
: ${outFolder:="out.NW"}




#********
# BEGIN *
#********


if [[  "$rep" == ""  ]]
then
	~/software/R-3.5.2/bin/Rscript /no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/bin/NW.alignment.permutations.R -i $input -b $SGE_TASK_ID -d $outFolder --t1 12 --t2 12 -p 100000
else
	~/software/R-3.5.2/bin/Rscript /no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/bin/NW.alignment.permutations.R -i $input -b $SGE_TASK_ID -d $outFolder --t1 12 --t2 12 -p 100000 -r "$rep"
fi


