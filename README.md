# PoxHost
Predicting the orthopoxvirus-mammal network with viral genetic features

## Description
This repository contains data and code to reproduce all analyses from <...>.

## Instructions 
To run the analysis, you can fork and clone this repository or download the "Tseng2020" folder directly to your local desktop (see *File organization (recommended)*). The code is separated into two markdown files, each contained within their corresponding folders: one for the host prediction model – *HostPrediction_Code.Rmd* – and the other for the link prediction model – *LinkPrediction_Code.Rmd*. Both markdown files are organized similarly into five parts (see *Code organization*) and draw data from the same source file, *Data_raw.RData*. 

### File organization (recommended)
1. Tseng2020 > Host Trait Model > ... 
      - HostTraitModel_Code.Rmd
      - HostTraitModel_Code.pdf
      - PoxHost_RawData.RData
      - MAMMALS.shp (*see NOTES*)
      - Output/ (*see NOTES*)
2. Tseng 2020 > Link Prediction Model > ...
      - LinkPredictionModel_Code.Rmd
      - LinkPredictionModel_Code.pdf
      - PoxHost_RawData.RData
      - MAMMALS.shp (*see NOTES*)
      - Output/ (*see NOTES*)

*NOTE: 
- MAMMALS.shp is a large shape file (>1GB) of mammal geographical range required in the last section of the code ("5. Mapping host distributions"). We recommend you download and save it directly to the the working directory of your local desktop. The file can be obtained from the IUCN Red List Spatial Database <https://www.iucnredlist.org/resources/spatial-data-download>.
- Before proceeding to run the code, we recommend you create an "Output" sub-folder (e.g., ~/Tseng2022/#### Prediction Model/Output/) contained within each model folder. The code for both models will save all output (e.g., cleaned datasets, model output, figures and tables) to the corresponding "Output" folders. 

### Code organization
For *HostTraitModel_Code.Rmd*:
1. Data Preparation
     - Input: *PoxHost_RawData.RData*
     - Output: *HostTraitModel_CleanData.RData*
2. Phylogenetic analysis
     - Input: *HostTraitModel_CleanData.RData*
     - Output: Figure1
3. Boosted regression trees (BRT)
     - Input: *HostTraitModel_CleanData.RData*
     - Output: FigureS1, TableS1, *par_tuning_data_summary.csv*, FigureS2, *HostData_results.RData*
4. BRT figures 
     - Input: *pcr_brts.rds*, *comp_brts.rds*, *pm_brts.rds*, *HostTraitModel_CleanData.RData*
     - Output: FigureS3, Figure2, TableS5, FigureS4, *PoxHost_predictions.csv*, Figure3, TableS6.csv
5. Mapping host distributions
     - Input: *PoxHost_predictions.csv*, *MAMMALS.shp*
     - Output: Figure4

For *LinkPredictionModel_Code.Rmd*:
1. Dimensionality Reduction
     - Input: *PresenceAbsence_OPV.xlsx*
     - Output: *pc_genes.RData*
2. Data Preparation
     - Input: *PoxHost_RawData.RData* & *pc_genes.RData*
     - Output: *LinkData_clean.RData*
3. Boosted regression trees (BRT)
     - Input: *LinkData_clean.RData*
     - Output: FigureS1, TableS1, *par_tuning_data_summary.csv*, FigureS2, *LinkData_results.RData*
4. BRT figures 
     - Input: *LinkData_results.RData*, *LinkData_clean.RData*
     - Output: FigureS3, Figure2, TableS5, FigureS4, *PoxHost_predictions.csv*, Figure3, TableS6.csv
5. Mapping host distributions
     - Input: *PoxHost_predictions.csv*, *MAMMALS.shp*
     - Output: Figure4
     
### Running the BRT models on an HPC cluster
To run the BRT models on a high performance computing (HPC) cluster (highly recommended), sample scripts and examples are available in the folder ~/Tseng2022/[HPC Example](https://github.com/viralemergence/PoxHost/tree/main/Tseng2022/HPC%20Example).
