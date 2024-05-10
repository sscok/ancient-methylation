#!/bin/bash -l

intersectBed=/usr/local/sw/bedtools2-2.25.0/bin/intersectBed							# /path/to/bedtools/intersect
epipaleo=$1	# /path/to/epiPALEOMIX/outputs/
output=$2	# /output/directory/name/
refanno=$3	# /path/to/reference/annotation/file		

cd ${epipaleo}

mkdir ${output}
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
	cd ${output}
	# replicate.sh; replicate the positions according to MS to produce reads (deaminated: 1, not deaminated 0)
	individuals=( UstIshim Loschbour prs002 prs009 prs013 prs016 R2 R3 R9 R7 R15 Motala12 K14 Sf12 cta016 irk025 irk061 kra001 LBK AKT16 Asp6 BAR25 Dil16 Ess7 Herx Klein7 LEPE48 LEPE52 Nea2 Nea3 STAR1 VC3-2 VLASA32 VLASA7 )
	output_names=( UstIshim_replicated Loschbour_replicated prs002_replicated prs009_replicated prs013_replicated prs016_replicated R2_replicated R3_replicated R9_replicated R7_replicated R15_replicated Motala12_replicated K14_replicated Sf12_replicated cta016_replicated irk025_replicated irk061_replicated kra001_replicated LBK_replicated AKT16_replicated Asp6_replicated BAR25_replicated Dil16_replicated Ess7_replicated Herx_replicated Klein7_replicated LEPE48_replicated LEPE52_replicated Nea2_replicated Nea3_replicated STAR1_replicated VC3-2_replicated VLASA32_replicated VLASA7_replicated )
	group=( HG HG NE NE NE NE NE NE NE HG HG HG HG HG HG HG HG HG NE NE NE NE NE NE NE NE NE NE NE NE NE NE HG HG )
	tissue=( Bone Tooth Tooth Tooth Tooth Tooth Bone Bone Bone Bone Bone Tooth Bone Bone Tooth Tooth Tooth Tooth Tooth Bone Bone Bone Bone Bone Bone Bone Bone Bone Bone Bone Bone Bone Bone Bone )
	genetic_sex=( XY XY XX XY XY XY XX XX XY XY XY XY XY XX XX XY XY XY XX XX XY XY XY XY XX XX XY XY XX XX XX XY XY XY )
	lab=( Fu Lazaridis Quinto Quinto Quinto Quinto Antonio Antonio Antonio Antonio Antonio Lazaridis Seguin Gunther Kilinc Kilinc Kilinc Kilinc Lazaridis Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi Marchi )

	for i in "${!individuals[@]}"
	do
		python3 ../../replicate.py "${individuals[i]}" "${output_names[i]}" "${individuals[i]}" "${group[i]}" "${tissue[i]}" "${genetic_sex[i]}" "${lab[i]}"
	done