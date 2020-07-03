#*****Feb 24th, 2020*****

# run cmnd.txt 1-3
qsub -q rg-el7 -N H3K4me2.QN.merged.bedfiles.1-3 -m ea -M beatrice.borsari@crg.eu -cwd -o ../../logs/$JOB_NAME -e ../../errors/$JOB_NAME cmnd.txt


#*****Feb 25th, 2020*****

# run cmnd.txt 4
qsub -q rg-el7 -N H3K4me2.QN.merged.bedfiles.4 -m ea -M beatrice.borsari@crg.eu -cwd -o ../../logs/$JOB_NAME -e ../../errors/$JOB_NAME cmnd.txt
