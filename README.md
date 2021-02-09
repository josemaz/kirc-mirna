# Introduction

This repository contains code and supplementary materials for paper: *miRNA and gene expression and co-expression networks are strongly altered through stages in clear cell renal carcinoma*. 

Jose Maria Zamora-Fuentes, Jesus Espinal-Enriquez, Enrique Hernandez-Lemus



## Pre-requisites

Considerations:

- R (3.6.3)
- Python 3

Pre-requisites to run scripts in theses paheses are obtained with:

`$ Rscript pkgs-requiremnts.R`

Python Pre-requisites (optional):

`$ bash install-miniconda.sh`



## Directory structure

<!-- - */pipeline* : Files to fetch data of GDC. Make preQC, postQC and normalization
- */Results/DEG* : Output of Differential Expression of genes (DEG)
- */Results/Expression* : Clean Expression Matrix (Genes x Samples)
- */Results/MI* : 10000 biggest Mutual Information pairs of genes
- */Extras* : Data, Plots and Utils for this paper. -->

# RNA-SEQ

## 01 - Data Aquisition, Quality Control  and Normalization

These three phases of the process are performed by next R script:

`$ Rscript R/01-rnaseq.R`

The expression data for control and each stage are saved on:

`$ ls pipeline/TCGA/Outs/`

Plots for pre Normalization Quality Control are saved on:

`$ ls pipeline/TCGA/Plots/QC_PRE/`

Plots for post Normalization Quality Control are saved on:

`$ ls pipeline/TCGA/Plots/QC_POST/`



<!-- ## 02 - Mutual Information

You can use repository of parallel mulitcore [ARACNe](https://github.com/CSB-IG/ARACNE-multicore)

From ARACNe2 output writing out MI network cuts of 10k interactions into `Results`  directory in this repository:

`$ cp kirc-ctrl-10k.tsv kirc-stagei-10k.tsv kirc-stageii-10k.tsv kirc-stageiii-10k.tsv kirc-stageiv-10k.tsv Results/MI/`

YOu need



## 03 - Differential Expression (DEG)

Script to write out DEG in a *tsv* file:

`$ Rscript DEG.R`

Output for each stage is saved on:

`$ ls Results/DEG/`

This output contains volcano plot in html format.



## 04 - Networks Analysis by MI cut off interaction

We realized a network analysis and biological process cutting MI interactions in 100, 1K, 10K, 100K, 1M.

Pipeline to this analysis is done by:

`$ bash scripts/make-cuts.sh`

Steps are executed in the following order:

1. Intersections and diferences networks are saved on:

`$ ls Results/cuts-mi/100/intersections/*`

i.e., for cut off of 100 interactions in MI.

2. Heatmaps of intersections and diferences are saved on:

`$ ls Plot/heat-interacciones-100.png`

i.e., for cut off of 100 interactions in MI.

3. Venn diagrams of interactions between experimental groups are saved on:

`$ ls Plot/venn-100.png`

i.e., for cut off of 100 interactions in MI.

4. Calculation of enrichment for all groups are saved on:

`$ Results/cuts-mi/100/inter-all-groups-100.go.txt`

where last column is the component id asociated to its rows of genes enrichment.

And the file with genes asociated to its component is saved on:

`$ Results/cuts-mi/100/inter-all-groups-100.comp.txt`

5. The same notation was used to enrichment of only stages. We have those files saved on:

`$ Results/cuts-mi/100/inter-only-stages-100.go.txt`
`$ Results/cuts-mi/100/inter-only-stages-100.comp.txt`



## 05 - DEG Contrast analysis

To get genes (PLG and SLC) underexpressed in all DEG contrast:

`$ Rscript  R/contrasts-deseq.R`


## 06 - Plot Genes Under and Over expressed

In Rstudio you can use 

`$ Rscript  R/gene-boxplot.R`

## 07 - Enrichment by communities

We only enriched 1M networks with:

python Py/community-GO.py Results/cuts-mi/1M/inter-all-groups-1M.tsv

 -->

