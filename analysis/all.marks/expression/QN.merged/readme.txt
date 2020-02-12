#*****Jan 21st, 2020*****

# run 1-4 cmnd.txt
qsub -q rg-el7 -N expression.QN.merged.1-4 -m ea -M beatrice.borsari@crg.eu -cwd -o ../logs/$JOB_NAME -e ../errors/$JOB_NAME cmnd.txt

# run 5 cmnd.txt manually (need to add and then remove gene_id)

# run 6-10 cmnd.txt
qsub -q rg-el7 -N expression.QN.merged.6-10 -m ea -M beatrice.borsari@crg.eu -cwd -o ../logs/$JOB_NAME -e ../errors/$JOB_NAME cmnd.txt


