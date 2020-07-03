#*****Apr 23rd, 2020*****

# run 1-8 cmnd.txt
qsub -q rg-el7 -N cEBPa.QN.merged.1-8 -m ea -M beatrice.borsari@crg.eu -cwd -o ../logs/$JOB_NAME -e ../errors/$JOB_NAME cmnd.txt
