#!/bin/bash -l
#SBATCH --partition=bonobo  			# write partition name(chimp;bonobo;macaque$i)
#SBATCH --job-name=ANOVA                 	# write your job name
#SBATCH --ntasks=1               		# number of tasks
#SBATCH --output=ANOVA.txt
#SBATCH --error=ANOVA.err


cd /mnt/NAS/projects/2019_epipaleomix/epipaleo/out_shotgun/epipaleo_out/anova_marchi/downsampled/ds_rep

for i in {1..20} ; do
	Rscript --verbose --vanilla ../../../../anova.R all_replicated_ds_sorted_a_$i /mnt/NAS/projects/2019_epipaleomix/epipaleo/out_shotgun/epipaleo_out/anova_marchi/downsampled/ds_rep/factor_Individual/anova_summaries_ind_$i /mnt/NAS/projects/2019_epipaleomix/epipaleo/out_shotgun/epipaleo_out/anova_marchi/downsampled/ds_rep/factor_Individual/genes_ind_$i
done

exit
