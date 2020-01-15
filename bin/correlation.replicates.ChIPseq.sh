#!/bin/bash


#********
# USAGE *
#********

display_usage() {
	echo -e "DESCRIPTION: computes correlation between 2 replicates of a ChIP-seq sample provided a bedfile of genomic regions\n"
	echo -e "\t--mark <ChIP-seq target>\n"
	echo -e "\t--sample <experiment Id>\n"
	echo -e "\t--db <dashboard_file> (i.e. chipseq.db)\n"
	echo -e "\t--loci <bedFile>\n"
	echo -e "\t--method <pearson/spearman> (default: pearson)\n"
	echo -e "\t--log <yes/no> (default: no)\n"
	echo -e "\t--window <500/1000/2000> (i.e. size of the genomic regions)\n"
	echo -e "\t--type <pileupSignal/pvalueSignal/fcSignal> (default: pileupSignal)\n"
	echo -e "\t--rep1 <rep1_mergedId> (default: samplemarkX1)\n"
	echo -e "\t--rep2 <rep2_mergedId> (default: samplemarkX2)\n"
	echo -e "\t--keep <yes/no> (default: no)\n"
	echo -e "\t--verbose <yes/no> (default: no)\n"
} 


if [[  $1 == "--help" ||  $1 == "-h" ]]
then
    	display_usage
        exit 0
fi




if [  $# -le 9 ]
then
	echo "ERROR: insufficient number of arguments\n"
    	display_usage
	exit 1
fi




#***************
# READ OPTIONS *
#***************

while [[ $# -gt 1 ]]; do

	key="$1"
	
	case $key in
    	
	--mark)
    	mark="$2"
    	shift # past argument
    	;;
    	
	--sample)
    	sample="$2"
    	shift # past argument
    	;;
    	
	--db)
    	db="$2"
    	shift # past argument
    	;;
    	
	--loci)
    	loci="$2"
    	shift # past argument
    	;;

	--method)
	method="$2"
	shift # past argument
	;;

	--log)
	log="$2"
	shift # past argument
	;;

	--window)
	window="$2"
	shift # past argument
	;;

	--type)
	type="$2"
	shift # past argument
	;;
	
	--rep1)
	rep1="$2"
	shift # past argument
	;;

	--rep2)
	rep2="$2"
	shift # past argument
	;;

	--keep)
	keep="$2"
	shift # past argument
	;;

	--verbose)
	verbose="$2"
	;;
	*)
	
	;;
	esac
	shift
done

: ${log:="no"}
: ${type:="pileupSignal"}
: ${method:="pearson"}
: ${rep1:="$sample""$mark""X1"}
: ${rep2:="$sample""$mark""X2"}
: ${keep:="no"}
: ${verbose:="no"}


if [[ "$method" != "pearson" && "$method" != "spearman" ]]
then
    	echo "ERROR: 'method' must be either <pearson> or <spearman>"
        exit 1
fi


if [[ "$log" != "no" && "$log" != "yes" ]]
then
    	echo "ERROR: 'log' must be either <no> or <yes>"
        exit 1
fi


if [[ "$window" != 500 && "$window" != 1000 && "$window" != 2000 ]]
then
	echo "ERROR: 'window' must be either <500>, <1000> or <2000>"
	exit 1
fi


if [[ "$type" != "pileupSignal" && "$type" != "pvalueSignal" && "$type" != "fcSignal" ]]
then
	echo "ERROR: 'type' must be a .bw file type"
	exit 1
fi 




if [[ "$verbose" == "yes" ]]
then

        echo "reading options .."
        echo "ChIP-seq target: " "${mark}"
        echo "ChIP-seq sample: " "${sample}"
        echo "dashboard file: " "${db}"
        echo "bedFile genomic regions: " "${loci}"
        echo "correlation method: " "${method}"
        echo "applying log: " "${log}"
        echo "genomic window: " "${window}"
        echo "bigWig data type: " "${type}"
        echo "merged ID rep.1: " "${rep1}"
        echo "merged ID rep.2: " "${rep2}"
	echo -e "keep tmp files: " "${keep}\n"

fi




#***********************
# COMPUTE CORRELATIONS *
#***********************


if [[ "$verbose" == "yes" ]]
then
	echo -e "computing correlations .."
fi

awk -v rep1="$rep1" -v type="$type" 'BEGIN{FS=OFS="\t"}
$1==rep1 && $5==type{print $2}' $db | while read line; do
	bwtool summary $loci -header -keep-bed $line "$rep1"."$(basename $db .db)"."$window"."$method"."$log".bwtool.summary.tsv
done

awk -v rep2="$rep2" -v type="$type" 'BEGIN{FS=OFS="\t"}
$1==rep2 && $5==type{print $2}' $db | while read line; do
	bwtool summary $loci -header -keep-bed $line "$rep2"."$(basename $db .db)"."$window"."$method"."$log".bwtool.summary.tsv
done


~abreschi/utils/join.py -b <(awk 'NR>1{print $1"_"$2"_"$3"\t"$8}' "$rep1"."$(basename $db .db)"."$window"."$method"."$log".bwtool.summary.tsv) -a <(awk 'NR>1{print $1"_"$2"_"$3"\t"$8}' "$rep2"."$(basename $db .db)"."$window"."$method"."$log".bwtool.summary.tsv) > "$sample""$mark"."$(basename $db .db)"."$window"."$method"."$log".bwtool.summary.tsv
sed -i '1irep1\trep2' "$sample""$mark"."$(basename $db .db)"."$window"."$method"."$log".bwtool.summary.tsv

if [[ "$log" == "yes" ]]
then
	coeff=$(~abreschi/R/scripts/matrix_to_dist.R -i "$sample""$mark"."$(basename $db .db)"."$window"."$method"."$log".bwtool.summary.tsv -c $method -l -p 0.0001 | awk 'NR==2{print $3}')
else
	coeff=$(~abreschi/R/scripts/matrix_to_dist.R -i "$sample""$mark"."$(basename $db .db)"."$window"."$method"."$log".bwtool.summary.tsv -c $method | awk 'NR==2{print $3}')
fi

echo -e "$mark\t$window\t$method\t$log\t$sample\t$coeff"


#*******************
# REMOVE TMP FILES *
#*******************


if [[ "$keep" == "no" ]]
then
	rm "$rep1"."$(basename $db .db)"."$window"."$method"."$log".bwtool.summary.tsv "$rep2"."$(basename $db .db)"."$window"."$method"."$log".bwtool.summary.tsv "$sample""$mark"."$(basename $db .db)"."$window"."$method"."$log".bwtool.summary.tsv
	if [[ "$verbose" == "yes" ]]
	then
		echo "removing tmp files .."
	fi
fi
