#!/bin/bash -l
#SBATCH --partition=macaque                    # write partition name(chimp;bonobo;macaque$i)
#SBATCH --job-name=vioplots                  # write your job name
#SBATCH --ntasks=1                              # number of tasks
#SBATCH --output=vioplots.txt
#SBATCH --error=vioplots.err




Rscript --verbose --vanilla violin_plots_figure1B.R /mnt/NAS/projects/2019_epipaleomix/epipaleo/out_shotgun/with_marchi/ratios_10



exit
