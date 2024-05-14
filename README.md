
# Investigating food production-associated DNA methylation changes in paleogenomes: lack of consistent signals beyond technical noise

## Pipeline Usage
A file that operates all the steps is called pipeline.sh
One can use it as follows:

```
bash pipeline.sh output_folder_name processed_outputs_folder_name /path/to/sorted/annotation/file.bed /path/to/sorted/reference/cpgsite.bed /path/to/CGIs/dataset /path/to/scripts/
```

## Descriptions Related to all Steps

### Generation of replicated datasets

First, one should generate the annotated methylation score (MS) datasets (i.e. including gene names) using the script below.

```
Rscript --verbose --vanilla process_and_filter_epiPALEOMIX.R
```

After that the output of these files are used in the script:

```
bash replicated_dataset_generation.sh output_folder_name processed_outputs_folder_name /path/to/sorted/annotation/file.bed /path/to/scripts/
```

This uses replicate.py script to generate the desired output for replicated dataset generation.

### Downsampled Datasets Generation

20 Downsampled datasets are generated using downsampled_dataset_generation.R

by navigating to the directory which includes generated replicated datasets per individual.

Note: the names *_replicated should be replaced by the actual file names generated when reproducing the results.

  
```
Rscript --verbose --vanilla downsampled_dataset_generation.R output_folder_name/processed_outputs_folder_name
```
After generating these files, one can merge all the individual files by using cat command.

Then, generate a file that has the mean Methylation Scores per individual per gene using the script below.

```
Rscript --vanilla reshape_replicated_dataset.R
Rscript --vanilla downsampled_means_per_gene.R
```

### Analysis of Variance on Replicated Dataset

  

One can use the ANOVA.R file to carry out the Analysis of Variance on our 20 downsampled datasets by simply using the command below.

  
```
Rscript --vanilla ANOVA.R /path/to/replicated/downsampled/dataset /path/to/anova/summary/output/file /path/to/resulting/genes
```

This script is called for all 20 downsampled replicated files.


### Matrix Generation

The second type of MS dataset aside from the replicated MS dataset is an annotated matrix with CpG positions and individuals is used for different analyses. The CpG site reference dataset (it must be sorted by position) is merged with the processed epiPALEOMIX output of every single individual and annotated. Missing MS values are represented by NA. The matrix is generated by using the script:

  
```
Rscript --vanilla join_matrix.R /path/to/sorted/reference/cpgsite.bed
```

### Plotting the Figures

  

Below script plots the Violin plots shown in Figure 1. It uses the generated matrix which is described in the previous section as input.

  
```
Rscript --vanilla violin_plots_figure1B.R
```

The scatter plots shown in Figure 1 can be replicated via running the next script. If one investigates the script, they will see a special type of input file is required for this script (CGIs). This file includes the genomic areas categorized merged with our MS dataset. A preview of the dataset is included inside the example_previews/ folder.

```
Rscript --vanilla mean_MS_per_individual_on_genomic_areas.R /path/to/CGIs/dataset
```

 
In Figure 2, we see the plots of candidate genes for each category (subsistence type, tissue type and genetic sex) and for each model described in the Methods section of the article. "/path/to/genes" is located here: data/genes/. The figure can be plotted using the following script:

```
Rscript --vanilla scatterplots_of_candidate_genes.R /path/to/genes
```

The scripts inside this R file require the input files with the information related to the selected genes. They are basically the grepped positions from the generated replicated datasets.

Multi-dimensional Scaling (Figure 3) analysis was carried out using the following file. It requires the means of all 20 downsampled datasets with respect to observed genes and information file related to the individuals. The info file is located here: data/Shotgun_inds.tsv. One can run the script using the command below. 

```
Rscript --vanilla multi_dimensional_scaling.R info_file.tsv
```

X-chromosome analysis just uses the normalized MSs and basic violin plots are realized after (Figure 4).

### Other analyses mentioned

All the other analyses mentioned in the article are carried out via using the scripts inside the general_analysis.R file.


For further information, please refer to the article.

A preview of all the datasets mentioned in the codes are provided under the "example_previews/" directory.
