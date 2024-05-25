
# Investigating food production-associated DNA methylation changes in paleogenomes: lack of consistent signals beyond technical noise

## Pipeline usage

The file that operates all the steps is called pipeline.sh. One can use it as follows:

```
bash pipeline.sh output_folder_name processed_outputs_folder_name /path/to/sorted/annotation/file.bed /path/to/sorted/reference/cpgsite.bed /path/to/CGIs/dataset /path/to/scripts/
```

## Descriptions related to all steps

### Generating replicated (recoded) datasets

First, one should generate the annotated methylation score (MS) datasets (i.e. including gene names) using the script below using the epiPALEOMIX software (https://github.com/KHanghoj/epiPALEOMIX). 

```
Rscript --verbose --vanilla process_and_filter_epiPALEOMIX.R
```

The output of these files is used in the script:

```
bash replicated_dataset_generation.sh output_folder_name processed_outputs_folder_name /path/to/sorted/annotation/file.bed /path/to/scripts/
```

This uses the replicate.py script to recode the data so that each row describes the state of one observation (deaminated or not) per CpG position per library. 

In addition to deamination, the columns include information about the library (tissue, sex, subsistence, lab, and gene). 

This dataset is referred to as "Full data for linear mixed models" in the article. See the example file (replicated-dataset-preview).


### Generating downsampled (normalized) datasets and merging

The script downsampled_dataset_generation.R generates 20 downsampled / normalised datasets for each library. The goal here is to ensure all libraries processed have the same average deamination rate. 

The script is run by navigating to the directory which includes generated replicated datasets per individual.

Note: the names *_replicated should be replaced by the actual file names generated when reproducing the results.
  
```
Rscript --verbose --vanilla downsampled_dataset_generation.R output_folder_name/processed_outputs_folder_name
```

After generating these files, one can merge all the individual files using the cat command.

One can then generate a file with the mean MS per individual per gene using the scripts below:

```
Rscript --vanilla reshape_replicated_dataset.R
Rscript --vanilla downsampled_means_per_gene.R
```
This resulting dataset is referred in the article as "Gene-averaged data".
### Analysis of Variance (ANOVA) on the replicated (recoded) dataset

One can use the ANOVA.R file to carry out ANOVA on the downsampled (normalized) datasets using the command below:

```
Rscript --vanilla ANOVA.R /path/to/replicated/downsampled/dataset /path/to/anova/summary/output/file /path/to/resulting/genes
```

This script is called for all 20 downsampled (normalized) datasets separately.


### Matrix generation

The second type of MS dataset (in addition to the replicated MS dataset) is a matrix with MS values per CpG positions in the rows and individuals in the columns.

The CpG site reference dataset (must be sorted by position) is merged with the processed epiPALEOMIX output of every individual and is annotated. Missing MS values are represented by NA. The matrix is generated using the script:

```
Rscript --vanilla join_matrix.R /path/to/sorted/reference/cpgsite.bed
```


### Plotting figures  

The below script plots the violin plots shown in Figure 1. It uses the generated matrix described in the previous section as input.

  
```
Rscript --vanilla violin_plots_figure1B.R
```

The scatter plots shown in Figure 1 of our article (Çokoğlu et al.) can be replicated via running the next script. Genomic annotation files for CpG Islands (CGIs) is needed as input for this script. A preview of the dataset is included inside the example_previews/ folder.

```
Rscript --vanilla mean_MS_per_individual_on_genomic_areas.R /path/to/CGIs/dataset
```

 
In Figure 2 of our article (Çokoğlu et al.), we see the plots of candidate genes for each category (subsistence type, tissue type and genetic sex) and for each model described in the Methods section. The "/path/to/genes" is located here: data/genes/. The figure can be plotted using the following script:

```
Rscript --vanilla scatterplots_of_candidate_genes.R /path/to/genes
```

The scripts inside this R file require the input files with the information related to the selected genes. These are the grepped positions from the generated replicated datasets.

Multi-dimensional Scaling (MDS) analysis for Figure 3 was carried out using the following script. It requires the means of all 20 downsampled datasets with respect to observed genes and information file related to the individuals. The info file is located here: data/Shotgun_inds.tsv. One can run the script using the command below. 

```
Rscript --vanilla multi_dimensional_scaling.R info_file.tsv
```

The X-chromosome analyses in Figure 4 uses normalized MSs.

### Other analyses

All the other analyses described in the article were carried out via using the scripts inside the general_analysis.R file.

For further information, please refer to the article.

A preview of all the datasets mentioned in the codes are provided under the "example_previews/" directory.
