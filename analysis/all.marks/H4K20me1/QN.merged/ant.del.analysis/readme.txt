#******May 15th, 2020*******

# run cmnd.txt 1-12
qsub -q rg-el7 -N H4K20me1.QN.merged.ant.del.analysis.1-12 -m ea -M beatrice.borsari@crg.eu -cwd -o ../../logs/$JOB_NAME -e ../../errors/$JOB_NAME cmnd.txt
