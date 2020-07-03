#******Jan 21st, 2020*******

# run cmnd.txt 2-7
qsub -q rg-el7 -N H3K4me2.2-7 -m ea -M beatrice.borsari@crg.eu -cwd -o ./logs/$JOB_NAME -e ./errors/$JOB_NAME cmnd.txt


#******Feb 12th, 2020*******

# run cmnd.txt 8-9
qsub -q rg-el7 -N H3K4me2.8-9 -m ea -M beatrice.borsari@crg.eu -cwd -o ./logs/$JOB_NAME -e ./errors/$JOB_NAME cmnd.txt


#******Feb 23rd, 2020*******

# run cmnd.txt 10
qsub -q rg-el7 -l h_rt=18:00:00 -N H3K4me2.10 -m ea -M beatrice.borsari@crg.eu -cwd -o ./logs/$JOB_NAME -e ./errors/$JOB_NAME cmnd.txt
