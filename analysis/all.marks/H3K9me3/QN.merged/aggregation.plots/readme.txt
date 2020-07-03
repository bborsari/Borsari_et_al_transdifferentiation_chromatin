#*****Feb 24th, 2020*****

# run cmnd.txt 1-3
qsub -q rg-el7 -N H3K9me3.QN.merged.aggregation.plots.1-3 -m ea -M beatrice.borsari@crg.eu -cwd -o ../../logs/$JOB_NAME -e ../../errors/$JOB_NAME cmnd.txt
