#!/bin/bash -l
#SBATCH --partition=macaque  			# write partition name(chimp;bonobo;macaque$i)
#SBATCH --job-name=epiPALEOMIX                 	# write your job name
#SBATCH --ntasks=1               		# number of tasks
#SBATCH --output=epipaleomix.txt
#SBATCH --error=epipaleomix.err

Rscript --verbose --vanilla join_matrix.R


exit
