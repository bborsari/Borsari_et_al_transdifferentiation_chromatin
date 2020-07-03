#******May 7th, 2020********
qsub -q rg-el7 -N GM12878.TFs.gb -m ea -M beatrice.borsari@crg.eu -cwd -o /no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/ENCODE.data/logs/$JOB_NAME -e /no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/ENCODE.data/errors/$JOB_NAME cmnd.txt
