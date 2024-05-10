#!/bin/bash -l

intersectBed=/usr/local/sw/bedtools2-2.25.0/bin/intersectBed							# /path/to/bedtools/intersect
epipaleo=$1	# /path/to/epiPALEOMIX/outputs/
output=$2	# /output/directory/name/
refanno=$3	# /path/to/reference/annotation/file		
refcpg=$4	# /path/to/reference/cpg/positions
methcgi=$5	# /path/to/reference/cgi/shelf/shores
scripts=$6	# /path/to/scripts/

# process the input files
R CMD BATCH ${scripts}process_and_filter_epiPALEOMIX.R

# replicate the original files
bash ${scripts}replicated_dataset_generation.sh ${epipaleo} ${output} ${refanno}

# Monte-Carlo Method
R CMD BATCH ${scripts}downsampled_dataset_generation.R ${epipaleo}/${output}

cd ${epipaleo}/${output}
# prepare files for ANOVA
for i in {1..20} ; do
	for file in *_ds_$i; do
		echo $file
		cut -f 1,2 $file > chr_pos
		cut -f 2 $file > pos
		cut -f 3,4 $file > values
		paste chr_pos pos values > raw
		$intersectBed -a raw -b ${refanno} -wa -wb -f 1.0 > intersect
		awk '!seen[$1,$2]++' intersect > processed 
		cut -f 1,2,4,5,11 processed > $file
	done
	rm *pos* values processed intersect raw 
	# replicate.sh; replicate the positions according to MS to produce reads (deaminated: 1, not deaminated 0)
	individuals=( UstIshim Loschbour prs002 prs009 prs013 prs016 R2 R3 R9 R7 R15 Motala12 K14 Sf12 cta016 irk025 irk061 kra001 LBK AKT16 Asp6 BAR25 Dil16 Ess7 Herx Klein7 LEPE48 LEPE52 Nea2 Nea3 STAR1 VC3-2 VLASA32 VLASA7 )
	output_names=( UstIshim_replicated_ds_$i Loschbour_replicated_ds_$i prs002_replicated_ds_$i prs009_replicated_ds_$i prs013_replicated_ds_$i prs016_replicated_ds_$i R2_replicated_ds_$i R3_replicated_ds_$i R9_replicated_ds_$i R7_replicated_ds_$i R15_replicated_ds_$i Motala12_replicated_ds_$i K14_replicated_ds_$i Sf12_replicated_ds_$i cta016_replicated_ds_$i irk025_replicated_ds_$i irk061_replicated_ds_$i kra001_replicated_ds_$i LBK_replicated_ds_$i AKT16_replicated_ds_$i Asp6_replicated_ds_$i BAR25_replicated_ds_$i Dil16_replicated_ds_$i Ess7_replicated_ds_$i Herx_replicated_ds_$i Klein7_replicated_ds_$i LEPE48_replicated_ds_$i LEPE52_replicated_ds_$i Nea2_replicated_ds_$i Nea3_replicated_ds_$i STAR1_replicated_ds_$i VC3-2_replicated_ds_$i VLASA32_replicated_ds_$i VLASA7_replicated_ds_$i )
	group=( HG HG NE NE NE NE NE NE NE HG HG HG HG HG HG HG HG HG NE NE NE NE NE NE NE NE NE NE NE NE NE NE HG HG )
	tissue=( Bone Tooth Tooth Tooth Tooth Tooth Bone Bone Bone Bone Bone Tooth Bone Bone Tooth Tooth Tooth Tooth Tooth Bone Bone Bone Bone Bone Bone Bone Bone Bone Bone Bone Bone Bone Bone Bone )
	genetic_sex=( XY XY XX XY XY XY XX XX XY XY XY XY XY XX XX XY XY XY XX XX XY XY XY XY XX XX XY XY XX XX XX XY XY XY )
	lab=( Fu Lazaridis Quinto Quinto Quinto Quinto Antonio Antonio Antonio Antonio Antonio Lazaridis Seguin Gunther Kilinc Kilinc Kilinc Kilinc Lazaridis Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi )

	for i in "${!individuals[@]}"
	do
		python3 ${scripts}replicate.py "${individuals[i]}" "${output_names[i]}" "${individuals[i]}" "${group[i]}" "${tissue[i]}" "${genetic_sex[i]}" "${lab[i]}"
	done

done
# join all the replicated files
cat *_replicated_ds_$i > all_replicated_ds_$i

# ANOVA.sh
for i in {1..20} ; do
	R CMD BATCH ${scripts}ANOVA.R all_replicated_ds_$i anova_summaries_ind_$i genes_ind_$i
done

# process the anova results by p-values and turn them into tables

# reshape replicated dataset
R CMD BATCH ${scripts}reshape_replicated_dataset.R
R CMD BATCH ${scripts}downsampled_means_per_gene.R

# Matrix Generation

R CMD BATCH ${scripts}join_matrix.R ${refcpg}

# Plot the figures
mkdir plots

R CMD BATCH ${scripts}violin_plots_figure1B.R

R CMD BATCH ${scripts}mean_MS_per_individual_on_genomic_areas.R ${methcgi}

R CMD BATCH ${scripts}multi_dimensional_scaling.R data/Shotgun_inds.tsv

R CMD BATCH ${scripts}scatterplots_of_candidate_genes.R data/genes/

exit