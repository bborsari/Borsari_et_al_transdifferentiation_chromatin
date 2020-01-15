#!/bin/bash
#$ -N job.array.correlation.chipseq.replicates
#$ -cwd
#$ -q rg-el7
#$ -j y
#$ -o /no_backup/rg/bborsari/outputs/$JOB_NAME.$TASK_ID.log
#$ -t 1-648

input="$1"

options=$(cat $input | sed -n ${SGE_TASK_ID}p)
mark=$(echo "$options" | cut -f1)
sample=$(echo "$options" | cut -f2)
db=$(echo "$options" | cut -f3)
loci=$(echo "$options" | cut -f4)
method=$(echo "$options" | cut -f5)
log=$(echo "$options" | cut -f6)
window=$(echo "$options" | cut -f7)
type=$(echo "$options" | cut -f8)
rep1=$(echo "$options" | cut -f9)
rep2=$(echo "$options" | cut -f10)
outFile=$(echo "$options" | cut -f11)

/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/bin/correlation.replicates.ChIPseq.sh --mark "$mark" --sample "$sample" --db "$db" --loci "$loci" --method "$method" --log "$log" --window "$window" --type "$type" --rep1 "$rep1" --rep2 "$rep2" > $outFile
