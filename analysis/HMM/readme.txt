#******Jun 16th, 2020********
qsub -q rg-el7 -l h_rt=168:00:00 -l virtual_free=60G -N HMM -m ea -M beatrice.borsari@crg.eu -cwd -o logs/$JOB_NAME -e errors/$JOB_NAME cmnd.txt
