# Running a job using R script on an HPC cluster
The following Unix command lines (Terminal in Mac) provides a walk-through of how to submit an R job on Kamiak, Washington State University's high performance computing cluster. You can adapt the script to your HPC cluster.

## Description
To follow along with this example, download the folder ~/Tseng2022/[HPC Example](https://github.com/viralemergence/PoxHost/tree/main/Tseng2022/HPC%20Example) from the [PoxHost repository](https://github.com/viralemergence/PoxHost).  In this example, we execute part three (“3. Boosted regression trees (BRT)”) of the markdown file [HostPrediction_Code.Rmd](https://github.com/viralemergence/PoxHost/blob/main/Tseng2022/Host%20Prediction%20Model/HostPrediction_Code.Rmd). This selection of code builds boosted regression tree models trained on infection and competence data to predict mammal hosts of Orthopoxviruses. The following files required for this job are listed below: 

| Files                     | Description                                                                          |
| ------------------------- |------------------------------------------------------------------------------------- |
| Kamiak_BRT_07Sep2022.R    | Model code in R script |
| HostData_clean.RData      | Input data obtained from part one "1. Data preparation" of [HostPrediction_Code.Rmd](https://github.com/viralemergence/PoxHost/blob/main/Tseng2022/Host%20Prediction%20Model/HostPrediction_Code.Rmd) |
| "Output" folder           | Empty folder location for saving model output                                        |
| KatieJob_07Sep2022.sh     | Bash shell script (aka Slurm script) for submitting job                              |

## Instructions 
To run the analysis, you can fork and clone this repository or download the "Tseng2020" folder directly to your local desktop (see *File organization (recommended)*). The code is separated into two markdown files, each contained within their corresponding folders: one for the host prediction model – *HostPrediction_Code.Rmd* – and the other for the link prediction model – *LinkPrediction_Code.Rmd*. Both markdown files are organized similarly into five parts (see *Code organization*) and draw data from the same source file, *Data_raw.RData*. 

### File organization (recommended)
1. Tseng2020 > Host Prediction Model > ... 
      - HostPrediction_Code.Rmd
      - HostPrediction_Code.pdf
      - Data_raw.RData
      - MAMMALS.shp (*see NOTES*)
      - Output/ (*see NOTES*)
2. Tseng 2020 > Link Prediction Model > ...
      - LinkPrediction_Code.Rmd
      - LinkPrediction_Code.pdf
      - Data_raw.RData
      - MAMMALS.shp (*see NOTES*)
      - Output/ (*see NOTES*)

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

### Running the BRT models on an HPC cluster
To run the BRT models on a high performance computing (HPC) cluster (highly recommended), sample scripts and examples are available in the folder ~/Tseng2022/[HPC Example](https://github.com/viralemergence/PoxHost/tree/main/Tseng2022/HPC%20Example).

