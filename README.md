# PoxHost
Predicting the poxvirus-mammal network with viral genetic features

The folder "Tseng2022" contains data and code to reproduce all analyses from ......

# File organization
1. *Host Prediction Model*
  a) HostPrediction_Code.pdf
  b) HostPrediction_Code.Rmd: 
     I. Data Preparation
     II. Phylogenetic analysis
     III. Boosted regression trees (BRT)
     IV. BRT figures 
     V. Maps
  c) Data: 
     - Data_raw.RData: contains the raw datasets for cleaning in part I (Data Preparation)
     - HostData_clean.RData: contains the cleaned datasets produced at the end of part I 
  d) Output
     - All BRT model output, figures and tables are output into this folder. This folder contains a Git ignore pattern to excludes files output to this folder from the Git history due to file size constraints.
2. *Link Prediction Model*
3. *HPC Instructions*
