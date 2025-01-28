# PoxHost
Predicting the orthopoxvirus-mammal network with viral genetic features

## Description
This repository contains data and code to reproduce all analyses from "Viral genomic features predict orthopoxvirus reservoir hosts" by Tseng et al (2025).

## Instructions 
To run the analysis, you can fork and clone this repository or download the "Tseng2020" folder directly to your local desktop (see *File organization (recommended)*). The code is separated into two markdown files, each contained within their corresponding folders: one for the host prediction model – *HostPrediction_Code.Rmd* – and the other for the link prediction model – *LinkPrediction_Code.Rmd*. Both markdown files are organized similarly into five parts (see *Code organization*) and draw data from the same source file, *Data_raw.RData*. 

### Table of Contents 
1. PoxHost > data > ... 
      All source data can be found here with the exception of MAMMALS.shp, a large shape file (>1GB) of mammal geographical range required for mapping host distributions in the section titled "3. Mapping Host Distribution" (PoxHost > scripts > 3_Results&Figures.Rmd). We recommend you download this file separatel and save it directly to the the working directory of your local desktop. The file can be obtained from the IUCN Red List Spatial Database <https://www.iucnredlist.org/resources/spatial-data-download>.
   
3. PoxHost > figures > ...
      All figures, tables, and other outputs are saved to this file, which contain the following additional subfolders:
      > other > hosttrait
      > other > linkpred > pca > ...
      > supplementary > ...

4. PoxHost > renv > ...
      This folder contains the reproducible environment (package versions) for running the scripts in R.
6. PoxHost > scripts > ...
      All code for data cleaning, analysis, and visualizations can be found in the markdown files below and progress in the order as listed:
      > 1_HostTraitModel.Rmd
      > 2_LinkPredictionModel.Rmd
      > 3_Results&Figures.Rmd   

### Code organization (input/output data files)
For *HostTraitModel_Code.Rmd*:
1. Data Preparation
     - Input: *hosttrait_rawdata.RData*
     - Output: *hosttrait_cleandata.RData*
2. Phylogenetic analysis
3. BRT Model
     - Input: *hosttrait_cleandata.RData*
     - Output: *pcr_brts.rds*, *comp_brts.rds*, *pm_brts.rds*

For *LinkPredictionModel_Code.Rmd*:
1. Dimensionality Reduction
     - Input: *OPVnew_nowwithVirus.xlsx*
     - Output: *wgs_PCs.RData*
2. Data Preparation
     - Input: *linkpred_rawdata*, *wgs_PCs.RData*
     - Output: *linkpred_cleandata.RData*
3. BRT Model
     - Input: *LinkData_clean.RData*
     
### Running the BRT models on an HPC cluster
To run the BRT models on a high performance computing (HPC) cluster (highly recommended), sample scripts and examples are available in the folder ~/PoxHost/[hpc](https://github.com/viralemergence/PoxHost/tree/main/hpc).

### Please direct questions regarding model code to katie.tseng@wsu.edu ###
