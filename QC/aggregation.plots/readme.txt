#******Jan 13th, 2020*******

# run cmnd 1-12 manually

# run cmnd 13-16
qsub -q rg-el7 -N aggregation.plots.+-.5Kb -m bea -M beatrice.borsari@crg.eu -cwd -o ../logs/$JOB_NAME -e ../errors/$JOB_NAME cmnd.txt


#******Jan 27th, 2020*******

# run cmnd 17
qsub -q rg-el7 -N aggregation.plots.gene.body.+-.5Kb -m bea -M beatrice.borsari@crg.eu -cwd -o ../logs/$JOB_NAME -e ../errors/$JOB_NAME cmnd.txt
