# STFactor



## Introduction

STFactor is a method to identify underlying factors of spatially resolved transcriptomics

For more details, please read our paper: **Deciphering the spatial patterns and associated genes of factors underlying spatially resolved transcriptomics**


## Prerequisite
**R 4.0.3**

**Required libraries**: Seurat, stringr, ggplot2

**Platform**: Linux


## Data

A file of the count matrix:

  - a matrix of spots (rows) * genes (columns)
  - rownames: string of coordinates of the spots, x and y position are connected with a separator (e.g., '10x10', '12_16')
  - colnames: genes


Demonstration data is provided in 'data/':

  - 'Layer2_BC_count_matrix-1.tsv': count matrix of human breast cancer
  - 'Rep11_MOB_count_matrix-1.tsv': count matrix of mouse olfactory bulb

The count matrices are downloaded from https://www.spatialresearch.org/resources-published-datasets/


## Usage

```
> source('run_STFactor.r')
> STFactor(file_name, num_components, output_flag, file_sep, top_hvg, gene_filtering, loc_sep)
```

Parameters:
  - file_name: path of the count matrix file
  - num_components: number of underlying factors (maximum limit)
  - output_flag: a string used in the names of the output files to indicate the input data
  - file_sep: the separator used to read the count matrix file, default '\t
  - top_hvg: number of highly variable genes to select, default 2000
  - loc_sep: the separator between x and y position of the spot names, default 'x'
  - gene_filtering: used to remove genes expressed in less than gene_filtering locations, default 0.1


### Examples

```
> STFactor('data/Layer2_BC_count_matrix-1.tsv', 10, 'bc2')
```

```
> STFactor('data/Rep11_MOB_count_matrix-1.tsv', 10, 'mob11')
```

## Output

  - Spatial patterns of the underlying factors
    - Activities of the spatial patterns at each spot, saved in the text file: 'results/output_flag__factor_spatial_patterns.txt'. Rows are spots. Columns are underlying factors.
    - Visualization of the spatial patterns, plotted to 'results/output_flag_factor_spatial_patterns_visualization.pdf'.

  - Associated genes of the underlying factors

      Activities of the associated genes in each underlying factor, saved in 'results/output_flag_factor_associated_genes.txt'.


### Examples

Example output of the breast cancer data:

  - Spatial patterns of the underlying factors
    - Saved activities
    
    &nbsp;
    <img src = "https://user-images.githubusercontent.com/57746198/176135276-a6ede201-e4bb-4322-9978-c7323a349e98.png" width = 900>
    &nbsp;
    
    
    - Visualization
    
    &nbsp;
    <img src = "https://user-images.githubusercontent.com/57746198/175618820-8a530f42-cb84-4e2d-bde3-0063508ac721.png" width = 150>
    <img src = "https://user-images.githubusercontent.com/57746198/175618859-7b65f740-7bbb-4008-8d7f-6c2f797a07c2.png" width = 150>
    <img src = "https://user-images.githubusercontent.com/57746198/175618885-91c96c89-426f-496d-8459-bfe6b09fe44c.png" width = 150>
    <img src = "https://user-images.githubusercontent.com/57746198/175618910-5f5ad6fc-0cb3-4477-957a-1296b7f5a425.png" width = 150>
    &nbsp;
    

  - Associated genes of the underlying factors

    &nbsp;
    <img src = "https://user-images.githubusercontent.com/57746198/176133810-074c0469-968a-45d6-86ee-33a5efc02f3c.png" width = 250>
    &nbsp;
