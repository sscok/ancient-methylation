#!/bin/bash -l
#SBATCH --partition=macaque  			# write partition name(chimp;bonobo;macaque$i)
#SBATCH --job-name=makemeans                 	# write your job name
#SBATCH --ntasks=1               		# number of tasks
#SBATCH --output=mat.txt
#SBATCH --error=mat.err

Rscript --verbose --vanilla reshape_replicated_dataset.R


exit