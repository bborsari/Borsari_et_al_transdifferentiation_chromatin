#*****Jan 21st, 2020*****

# run 1-7 cmnd.txt
qsub -q rg-el7 -N H3K4me1.QN.merged.1-7 -m ea -M beatrice.borsari@crg.eu -cwd -o ../logs/$JOB_NAME -e ../errors/$JOB_NAME cmnd.txt


#*****Feb 13th, 2020*****

# run 8-11 cmnd.txt
qsub -q rg-el7 -N H3K4me1.QN.merged.8-11 -m ea -M beatrice.borsari@crg.eu -cwd -o ../logs/$JOB_NAME -e ../errors/$JOB_NAME cmnd.txt


#*****Feb 23rd, 2020*****

# run 12 cmnd.txt
qsub -q rg-el7 -N H3K4me1.QN.merged.12 -m ea -M beatrice.borsari@crg.eu -cwd -o ../logs/$JOB_NAME -e ../errors/$JOB_NAME cmnd.txt


#*****Apr 7th, 2020******

# re-run 10 cmnd.txt
qsub -q rg-el7 -N H3K4me1.QN.merged.10 -m ea -M beatrice.borsari@crg.eu -cwd -o ../logs/$JOB_NAME -e ../errors/$JOB_NAME cmnd.txt
