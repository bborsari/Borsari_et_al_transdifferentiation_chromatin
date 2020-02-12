#*****Jan 21st, 2020*****

# run 1-7 cmnd.txt
qsub -q rg-el7 -N H3K4me2.QN.merged.1-7 -m ea -M beatrice.borsari@crg.eu -cwd -o ../logs/$JOB_NAME -e ../errors/$JOB_NAME cmnd.txt


