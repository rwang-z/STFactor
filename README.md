# STFactor

&nbsp;

## Introduction

STFactor is a method to identify underlying factors of spatially resolved transcriptomics

For more details, please read our paper: **Bayesian factorization of spatially resolved transcriptomics reveals underlying factors of spatial gene expression**.

&nbsp;

## Prerequisite
**R 4.0.3**

**Required libraries**: Seurat, stringr, ggplot2

**Platform**: Linux

&nbsp;

## Data

A file of the count matrix:

  - a matrix of spots (rows) * genes (columns)
  - rownames: string of coordinates of the spots, x and y position are connected with a separator (e.g., '10x10', '12_16')
  - colnames: genes

&nbsp;

Demonstration data is provided in 'data/':

  - 'Layer2_BC_count_matrix-1.tsv': count matrix of human breast cancer
  - 'Rep11_MOB_count_matrix-1.tsv': count matrix of mouse olfactory bulb

The count matrices are downloaded from https://www.spatialresearch.org/resources-published-datasets/

&nbsp;


## Usage

```
> source('run_STFactor.r')
> STFactor(file_name, num_components, output_flag, file_sep, top_hvg, gene_filtering, loc_sep)
```

Parameters:
  - file_name: path of the count matrix file
  - num_components: number of underlying factors (maximum limit)
  - output_flag: a string used in the names of the output files to indicate the input data
  - file_sep: the separator used to read the count matrix file, default '\t'
  - top_hvg: number of highly variable genes to select, default 2000
  - loc_sep: the separator between x and y position of the spot names, default 'x'
  - gene_filtering: used to remove genes expressed in less than gene_filtering locations, default 0.1

&nbsp;

### Examples

```
> STFactor('data/Layer2_BC_count_matrix-1.tsv', 10, 'bc2')
```

```
> STFactor('data/Rep11_MOB_count_matrix-1.tsv', 10, 'mob11')
```

&nbsp;

## Output

  - Spatial patterns of the underlying factors
    - Activities of the spatial patterns at each spot, saved in the text file: 'results/output_flag__factor_spatial_patterns.txt'. Rows are spots. Columns are underlying factors.
    - Visualization of the spatial patterns, plotted to 'results/output_flag_factor_spatial_patterns_visualization.pdf'.

  - Associated genes of the underlying factors

      Activities of the associated genes in each underlying factor, saved in 'results/output_flag_factor_associated_genes.txt'.

&nbsp;

### Examples

Example output of the breast cancer data:

  - Spatial patterns of the underlying factors
    - Saved activities
    
    &nbsp;
    <img src = "https://user-images.githubusercontent.com/57746198/176135276-a6ede201-e4bb-4322-9978-c7323a349e98.png" width = 900>
    &nbsp;
    
    
    - Visualization
    
    &nbsp;
    <img src = "https://user-images.githubusercontent.com/57746198/176137397-68d1be15-0eb0-465d-a6e6-99347986743c.png" width = 150>
    <img src = "https://user-images.githubusercontent.com/57746198/176137563-b7d2c263-9639-46c0-bec3-92ed9915213d.png" width = 150>
    <img src = "https://user-images.githubusercontent.com/57746198/176137635-8e245d8c-685d-4aa2-b9bf-a4656e06ae65.png" width = 150>
    <img src = "https://user-images.githubusercontent.com/57746198/176137769-daebfe80-6d7d-46f2-8836-877f912f7117.png" width = 150>
    <img src = "https://user-images.githubusercontent.com/57746198/176137836-c0ed4af2-5ab4-42c4-9648-118f72d12847.png" width = 150>
    
    <img src = "https://user-images.githubusercontent.com/57746198/176137885-1c6a8f13-bca1-4a51-9539-f12d54818544.png" width = 150>
    <img src = "https://user-images.githubusercontent.com/57746198/176137952-a4d22bfe-edbe-48f9-b959-29f57d2306b0.png" width = 150>
    <img src = "https://user-images.githubusercontent.com/57746198/176138014-91a5f259-2e07-4bdc-8c3a-a7da82b0133d.png" width = 150>
    <img src = "https://user-images.githubusercontent.com/57746198/176138079-f6c31f6c-ccaa-4d1f-80f3-62fed7885284.png" width = 150>
    <img src = "https://user-images.githubusercontent.com/57746198/176138145-d36b63d0-74a9-4ddb-b5dd-653f13cff00b.png" width = 150>
    &nbsp;
    

  - Associated genes of the underlying factors

    &nbsp;
    <img src = "https://user-images.githubusercontent.com/57746198/176133810-074c0469-968a-45d6-86ee-33a5efc02f3c.png" width = 250>
    &nbsp;
