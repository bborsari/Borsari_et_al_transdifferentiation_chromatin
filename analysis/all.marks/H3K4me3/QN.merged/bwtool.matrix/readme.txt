#*****Feb 24th, 2020*****

# run cmnd.txt 1-3
qsub -q rg-el7 -N H3K4me3.QN.merged.bwtool.matrix.1-3 -m ea -M beatrice.borsari@crg.eu -cwd -o ../../logs/$JOB_NAME -e ../../errors/$JOB_NAME cmnd.txt


#*****Feb 25th, 2020*****

# run cmnd.txt 4-7
qsub -q rg-el7 -N H3K4me3.QN.merged.bwtool.matrix.4-7 -m ea -M beatrice.borsari@crg.eu -cwd -o ../../logs/$JOB_NAME -e ../../errors/$JOB_NAME cmnd.txt
