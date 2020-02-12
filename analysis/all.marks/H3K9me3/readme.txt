#******Jan 21st, 2020*******

# run cmnd.txt 2-7
qsub -q rg-el7 -N H3K9me3.2-7 -m ea -M beatrice.borsari@crg.eu -cwd -o ./logs/$JOB_NAME -e ./errors/$JOB_NAME cmnd.txt
