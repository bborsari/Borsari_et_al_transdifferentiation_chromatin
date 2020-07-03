#*****Feb 14th, 2020******

# run cmnd.txt 1-2 manually

#*****Feb 23rd, 2020******

# run cmnd.txt 3-8
qsub -l virtual_free=10G -q rg-el7 -N H3K9me3.QN.merged.NW.pr.3-8 -m ea -M beatrice.borsari@crg.eu -cwd -o ../../logs/$JOB_NAME -e ../../errors/$JOB_NAME cmnd.txt
