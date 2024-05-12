#!/bin/bash -l

intersectBed=/usr/local/sw/bedtools2/bin/intersectBed							# /path/to/bedtools/intersect
epipaleo=$1	# /path/to/epiPALEOMIX/outputs/
output=$2	# /output/directory/name/
refanno=$3	# /path/to/reference/annotation/file		
refcpg=$4	# /path/to/reference/cpg/positions
methcgi=$5	# /path/to/reference/cgi/shelf/shores
scripts=$6	# /path/to/scripts/

mkdir ${epipaleo}
# process the input files
Rscript --vanilla ${scripts}process_and_filter_epiPALEOMIX.R
# replicate the original files
bash ${scripts}replicated_dataset_generation.sh ${epipaleo} ${output} ${refanno} ${scripts}

# Monte-Carlo Method
Rscript --vanilla ${scripts}downsampled_dataset_generation.R ${epipaleo}/${output}

cd ${epipaleo}/${output}

# join all the replicated files
for i in {1..20} ; do
	cat *_replicated_ds_$i > replicated_ds_all_$i
done
# ANOVA.sh
for i in {1..20} ; do
	Rscript --vanilla ${scripts}ANOVA.R replicated_ds_all_$i anova_summaries_ind_$i genes_ind_$i
done

# process the anova results by p-values and turn them into tables

# reshape replicated dataset
Rscript --vanilla ${scripts}reshape_replicated_dataset.R
Rscript --vanilla ${scripts}downsampled_means_per_gene.R

# Matrix Generation

Rscript --vanilla ${scripts}join_matrix.R ${refcpg}

# Plot the figures
mkdir plots

Rscript --vanilla ${scripts}violin_plots_figure1B.R

Rscript --vanilla ${scripts}mean_MS_per_individual_on_genomic_areas.R ${methcgi}

Rscript --vanilla ${scripts}multi_dimensional_scaling.R ../../data/Shotgun_inds.tsv

Rscript --vanilla ${scripts}scatterplots_of_candidate_genes.R ../../data/genes/

exit
