# PoxHost
Predicting the orthopoxvirus-mammal network with viral genetic features

## Description
This repository contains data and code to reproduce all analyses from <...>.

## Instructions 
To run the analysis, you can fork and clone this repository or download the "Tseng2020" folder directly to your local desktop (see *File organization (recommended)*). The code is separated into two markdown files, each contained within their corresponding folders: one for the host prediction model – *HostPrediction_Code.Rmd* – and the other for the link prediction model – *LinkPrediction_Code.Rmd*. Both markdown files are organized similarly into five parts (see *Code organization*) and draw data from the same source file, *Data_raw.RData*. 

### File organization (recommended)
1. Tseng2020 > Host Prediction Model > ... 
      - HostPrediction_Code.Rmd
      - HostPrediction_Code.pdf
      - Data_raw.RData
      - MAMMALS.shp *see NOTES below*
      - Output/ *see NOTES below*
2. Tseng 2020 > Link Prediction Model > ...
      - LinkPrediction_Code.Rmd
      - LinkPrediction_Code.pdf
      - Data_raw.RData
      - MAMMALS.shp *see NOTES below*
      - Output/ *see NOTES below*

*NOTE: 
- MAMMALS.shp is a large shape file (>1GB) of mammal geographical range required in the last section of the code ("5. Mapping host distributions"). We recommend you download and save it directly to the the working directory of your local desktop. The file can be obtained from the IUCN Red List Spatial Database <https://www.iucnredlist.org/resources/spatial-data-download>.
- Before proceeding to run the code, we recommend you create an "Output" sub-folder (e.g., ~/Tseng2022/#### Prediction Model/Output/) contained within each model folder. The code for both models will save all output (e.g., cleaned datasets, model output, figures and tables) to the corresponding "Output" folders. 

### Code organization
1. Data Preparation
     - Input: *Data_raw.RData*
     - Output: *HostData_clean.RData*
2. Phylogenetic analysis
     - Input: *HostData_clean.RData*
     - Output: Figure1
3. Boosted regression trees (BRT)
     - Input: *HostData_clean.RData*
     - Output: FigureS1, TableS1, *par_tuning_data_summary.csv*, FigureS2, *HostData_results.RData*
4. BRT figures 
     - Input: *HostData_results.RData*, *HostData_clean.RData*
     - Output: FigureS3, Figure2, TableS5, FigureS4, *PoxHost_predictions.csv*, Figure3, TableS6.csv
5. Mapping host distributions
     - Input: *PoxHost_predictions.csv*, *MAMMALS.shp*
     - Output: Figure4

### Running the BRT models on an HPC (high performance computing) cluster

