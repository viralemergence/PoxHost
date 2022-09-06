# PoxHost
Predicting the orthopoxvirus-mammal network with viral genetic features

## Description
This repository contains data and code to reproduce all analyses from <...>.

## Instructions 
To run the analysis, you can fork and clone this repository or download the "Tseng2020" folder to your local desktop. The code is separated into two markdown files each contained within their corresponding model folders: one for the host prediction model – *HostPrediction_Code.Rmd* – and the other for the link prediction model – *LinkPrediction_Code.Rmd*. Both markdown files are organized similarly into five parts (see *File organization*) and draw data from the same source file, *Data_raw.RData*. All output (e.g., cleaned datasets, model output, figures and tables) are saved in the "Output" folder within each corresponding model folder. 

### File organization
1. Data Preparation
     - Input: *Data_raw.RData*
     - Output: *HostData_clean.RData*
2. Phylogenetic analysis
     - Input: *HostData_clean.RData*
     - Output: Figure 1
3. Boosted regression trees (BRT)
     - Input: *HostData_clean.RData*
     - Output: Figure S1, TableS1, par_tuning_data_summary.csv, Figure S2, *HostData_results.RData*
4. BRT figures 
     - Input: *HostData_results.RData*
     - Output: 
5. Maps
     - Input: *HostData_results.RData*
     - Output: 

### Running the analysis on an HPC (high performance computing) cluster

