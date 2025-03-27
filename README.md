# Molecular subtypes of high grade serous ovarian cancer across racial groups and gene expression platforms

## Overview
This repository and the subtype clustering repository ([link](https://github.com/greenelab/hgsc_rnaseq_clustering_pipeline)) contains all the code needed to recreate the analyses, tables, and figures in the paper "Molecular subtypes of high grade serous ovarian cancer across racial groups and gene expression platforms"
Black individuals with ovarian cancer experience poorer survival compared to non-Hispanic White individuals.
High-grade serous carcinoma (HGSC) is the most common ovarian cancer histotype, with gene expression subtypes that are associated with differential survival.
We characterized HGSC gene expression in Black individuals and considered whether gene expression differences by race may contribute to disparities in survival. 

We performed gene expression clustering using RNA-Seq data from Black and White individuals with HGSC, as well as array-based genotyping data from four existing studies of HGSC.
Our main analysis assigned subtypes by identifying dataset-specific clusters using K-means clustering for K=2-4.
The cluster- and dataset-specific gene expression patterns were summarized by moderated t-scores that differentiate an individual cluster from all other clusters within each dataset.
We compared the calculated gene expression patterns for each cluster across datasets by calculating the Pearson correlation between the two summarized vectors of moderated t-scores.
Following K=4 subtype assignment and mapping to The Cancer Genome Atlas (TCGA)-derived HGSC subtypes, we used multivariable-adjusted Cox proportional hazards models to estimate subtype-specific survival separately for each dataset. 

This repository contains all of the code used to quantify and QC the RNA-Seq samples before clustering, all downstream analyses after clustering, and all code to generate the figures and tables in the manuscript.
All clustering analyses can be re-run using Docker container for the clustering part of this project. Instructions to build the docker container are here: [link](https://github.com/greenelab/hgsc_rnaseq_clustering_pipeline).


## Data Availability
The processed data used in this analysis is included in the publication. Raw data will be made available through dbGaP under study ID: phs002262.v3.p1. 

## Code overview
- `data`: This is an empty folder, once the paper is published, we will release the processed data which can be put here.
- `reference_data`: This contains all the supplemental data needed to perform analyses. This also contains the main metadata files `reference_data/main_AA_metadata_table.tsv` and `reference_data/main_W_metadata_table.tsv`, which contain the cluster labels, as well as sample quality features.
- `figure_notebooks`: This contains code for generating most figures in the manuscipt
    - `figure_notebooks/compare_centroids.Rmd`: Figure 1
    - `figure_notebooks/all_method_comparison.Rmd`: Figure 2, Supp. Figure 9, 10
    - `figure_notebooks/make_qc_plots.Rmd`: Supp. Fig. 3-6
    - `figure_notebooks/K3_kmeans_vs_nmf.Rmd`: Supp. Fig 7
    - `figure_notebooks/rerun_clustering.Rmd`: Supp. Fig. 8
- `analyze_rnaseq`: This contains code for QC-ing the RNA-seq quantification. Below I delineate the main scripts used in this folder.
    - `analyze_rnaseq/analyze_rnaseq.Rmd`: This QC's the SchildkrautB samples.
    - `analyze_rnaseq/analyze_rnaseq_white.Rmd`: This QC's the SchildkrautW samples.
    - `analyze_rnaseq/make_log10_counts_for_way_pipeline.Rmd`: This makes the SchildkrautB counts to be used by the downstream clustering pipeline.
    - `analyze_rnaseq/make_log10_counts_for_way_pipeline_whites.Rmd`: This makes the SchildkrautW counts to be used by the downstream clustering pipeline.
    - `analyze_rnaseq/make_metadata_table.Rmd`: This generates the main metadata tables that are used in all downstream analyses.
    - `analyze_rnaseq/plot_utils.R`: Plotting methods
    - `analyze_rnaseq/test_SAM_performance.Rmd`: QC-ing the performance of SAM on RNA-Seq data.
- `quantify_rnaseq`: This contains code for RNA-seq quantification. 
- `survival`: This folder contains code for creating the tables and figures related to the survival analyses 
- `way_pipeline_downstream_analysis`: This folder contains code for comparing the labels of the same samples in the previous and updated Way pipeline code.

## Installation

This code uses renv ([link](https://rstudio.github.io/renv/articles/renv.html)) to manage the environment needed to run all the code within this repo. 
After installing R (tested with R version 4.1.2) and renv, and cloning this repository, type `renv::restore()` in order to restore the environment. 
