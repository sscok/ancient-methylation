#!/bin/bash -l
#SBATCH --partition=bonobo  			# write partition name(chimp;bonobo;macaque$i)
#SBATCH --job-name=replicate                 	# write your job name
#SBATCH --ntasks=1               		# number of tasks
#SBATCH --output=replicate.txt
#SBATCH --error=replicate.err


intersectBed=/usr/local/sw/bedtools2-2.25.0/bin/intersectBed							# /path/to/bedtools/intersect
epipaleo=$1	# /path/to/epiPALEOMIX/outputs/
output=$1	# /output/directory/name/
refanno=$3	# /path/to/reference/annotation/file		

cd /mnt/NAS/projects/2019_epipaleomix/epipaleo/out_shotgun/epipaleo_out/anova_marchi/downsampled/ds_rep


mkdir ${output}
# prepare files for ANOVA
for file in *; do    
	echo $file
	cut -f 1,2 $file > chr_pos
	cut -f 2 $file > pos
	cut -f 3,4 $file > values
	paste chr_pos pos values > raw
	$intersectBed -a raw -b ${refanno} -wa -wb -f 1.0 > intersect
	awk '!seen[$1,$2]++' intersect > processed 
	cut -f 1,2,4,5,11 processed > ${output}/$file
done
rm *pos* values processed intersect raw ${output}/${output} 
# replicate.sh; replicate the positions according to MS to produce reads (deaminated: 1, not deaminated 0)
individuals=( 14_ds_UstIshim 14_ds_Loschbour 14_ds_prs002 14_ds_prs009 14_ds_prs013 14_ds_prs016 14_ds_R2 14_ds_R3 14_ds_R9 14_ds_R7 14_ds_R15 14_ds_Motala12 14_ds_K14 14_ds_Sf12 14_ds_cta016 14_ds_irk025 14_ds_irk061 14_ds_kra001 14_ds_LBK 14_ds_AKT16 14_ds_Asp6 14_ds_BAR25 14_ds_Dil16 14_ds_Ess7 14_ds_Herx 14_ds_Klein7 14_ds_LEPE48 14_ds_LEPE52 14_ds_Nea2 14_ds_Nea3 14_ds_STAR1 14_ds_VC3-2 14_ds_VLASA32 14_ds_VLASA7 )
output_names=( UstIshim_replicated Loschbour_replicated prs002_replicated prs009_replicated prs013_replicated prs016_replicated R2_replicated R3_replicated R9_replicated R7_replicated R15_replicated Motala12_replicated K14_replicated Sf12_replicated cta016_replicated irk025_replicated irk061_replicated kra001_replicated LBK_replicated AKT16_replicated Asp6_replicated BAR25_replicated Dil16_replicated Ess7_replicated Herx_replicated Klein7_replicated LEPE48_replicated LEPE52_replicated Nea2_replicated Nea3_replicated STAR1_replicated VC3-2_replicated VLASA32_replicated VLASA7_replicated )
group=( HG HG NE NE NE NE NE NE NE HG HG HG HG HG HG HG HG HG NE NE NE NE NE NE NE NE NE NE NE NE NE NE HG HG )
tissue=( Bone Tooth Tooth Tooth Tooth Tooth Bone Bone Bone Bone Bone Tooth Bone Bone Tooth Tooth Tooth Tooth Tooth Bone Bone Bone Bone Bone Bone Bone Bone Bone Bone Bone Bone Bone Bone Bone )
genetic_sex=( XY XY XX XY XY XY XX XX XY XY XY XY XY XX XX XY XY XY XX XX XY XY XY XY XX XX XY XY XX XX XX XY XY XY )
lab=( Fu Lazaridis Quinto Quinto Quinto Quinto Antonio Antonio Antonio Antonio Antonio Lazaridis Seguin Gunther Kilinc Kilinc Kilinc Kilinc Lazaridis Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi )

for i in "${!individuals[@]}"
do
	python3 ../../../../replicate.py "${individuals[i]}" "${output_names[i]}" "${individuals[i]}" "${group[i]}" "${tissue[i]}" "${genetic_sex[i]}" "${lab[i]}"
done

# join all the replicated files
cat *_replicated > all_replicated_ds_14_a