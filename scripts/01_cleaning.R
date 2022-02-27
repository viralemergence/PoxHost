#PoxHost 01: Mammal orthopoxvirus (OPV) cleaning
#katie.tseng@wsu.edu
#2022-02-15

## clean environment & plots
rm(list=ls()) 
graphics.off()

## libraries
library(ape)
library(dplyr)
library(vroom)
library(tidyverse)

## set working directory
setwd("/Users/CDDEP/Downloads/OPV Host Prediction/PoxHost")

######################################### POXVIRUS DATA #########################################
## load in orthopoxvirus data
virion <- vroom("data/raw/Virion.csv.gz")

## view detection methods
virion %>%
  filter(VirusGenus == "orthopoxvirus") %>%
  pull(DetectionMethod) %>%
  unique()

## extract associations detected by PCR & Isolation only
poxdata <- virion %>% 
  filter(VirusGenus == "orthopoxvirus" & (DetectionMethod %in% c("PCR/Sequencing","Isolation/Observation","Antibodies"))) %>% 
  select(Host,DetectionMethod) %>%
  unique()
poxdata <- poxdata[!is.na(poxdata$Host),]

## subset detection by Antibodies
antibodies <- poxdata[which(poxdata$DetectionMethod=="Antibodies"),]
antibodies$Antibodies <- 1

## subset detection by PCR
pcr <- poxdata[which(poxdata$DetectionMethod=="PCR/Sequencing"),]
pcr$PCR <- 1

## subset detection by Isolation
competence <- poxdata[which(poxdata$DetectionMethod=="Isolation/Observation"),]
competence$Competence <- 1

## merge PCR & Isolation data
poxdata <- merge(pcr, competence, by="Host", all=TRUE)
poxdata <- merge(poxdata, antibodies, by="Host", all=TRUE)
poxdata <- subset(poxdata, select=-c(DetectionMethod,DetectionMethod.x,DetectionMethod.y))

## recode as binary (replace NA with zeroes)
poxdata$Antibodies=ifelse(is.na(poxdata$Antibodies),0,poxdata$Antibodies)
poxdata$PCR=ifelse(is.na(poxdata$PCR),0,poxdata$PCR)
poxdata$Competence=ifelse(is.na(poxdata$Competence),0,poxdata$Competence)

## create genus variable
poxdata$Host=str_to_title(poxdata$Host)
poxdata$gen=sapply(strsplit(as.character(str_to_title(poxdata$Host))," "),function(x) x[1])

  # # save file and clean environment
  # write.csv(poxdata,"data/cleaned/poxdata")  
  # rm(antibodies,pcr,competence,virion)

  ###### ALTERNATIVE FORMATTING #####
  ## Remove duplicate host rows giving preference to Isolation (over PCR)
  #   poxdata <- poxdata[order(poxdata$Host, poxdata$DetectionMethod),]
  #   poxdata_nodup <- poxdata[!duplicated(poxdata$Host),]
  ## Extract potential hosts associated with OPV detected by Isolation only
  #   poxdata_isol <- poxdata_nodup %>% 
  #   filter(DetectionMethod=="Isolation/Observation") %>% 
  #   select(Host,DetectionMethod) %>%
  #   unique()


## clean environment
rm(virion,antibodies,pcr,competence)

######################################### TAXONOMY DATA #########################################

## load in taxonomy (source: https://data.vertlife.org/)
taxa=read.csv('phylo/taxonomy_mamPhy_5911species.csv',header=T) 

## simplify to genus
gtaxa=taxa[!duplicated(taxa$gen),] #note: drops duplicates at the genus level
gtaxa=gtaxa[c('gen','fam','ord','clade','higher')]

## check that poxdata genus names match those in gdata
poxdata$gen[!poxdata$gen %in% gtaxa$gen]

## rename pox data taxonomy for "Piliocolobus badius" to "Procolobus badius" per IUCN taxonomy (https://www.iucnredlist.org/species/161247840/161259430)
poxdata$Host[data$gen=="Piliocolobus"] <- "procolobus badius"
poxdata$gen[data$gen=="Piliocolobus"] <- "Procolobus"

## combine taxonomy data with pox data (left-join)
data=merge(poxdata,gtaxa,by='gen',all.x=TRUE)

## clean environment
rm(poxdata,taxa,gtaxa)

######################################### PHYLOGENY DATA #########################################

## load Upham phylogeny 
# Source: https://doi.org/10.5061/dryad.tb03d03 (in subfolder "Data_S8_finalFigureFiles" > "_DATA")
tree=read.nexus('phylo/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')









## load in COMBINE trait data: https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.3344#support-information-section
combine=read.csv("data/raw/trait_data_imputed.csv",header=T)

## load in PanTHERIA trait data: https://esapubs.org/archive/ecol/E090/184/metadata.htm
pantheria =read.delim("data/raw/PanTHERIA_1-0_WR05_Aug2008.txt",header=T)

## redefine "999.00" (unknown info) records with "NA"
pantheria[pantheria==-999.00]<-NA

## keep PanTHERIA variables not captured by COMBINE database (see Trait Data Dictionary)
## Note: all variables in EltonTraits (see Trait Data Dictionary) were already captured in COMBINE database (source: https//figshare.com/collections/EltonTraits_1_0_Species-level_foraging_attributes_of_the_world_s_birds_and_mammals/3306933)
keeps <- c("MSW05_Order","MSW05_Family","MSW05_Genus","MSW05_Species","MSW05_Binomial","X2.1_AgeatEyeOpening_d","X18.1_BasalMetRate_mLO2hr","X5.2_BasalMetRateMass_g","X7.1_DispersalAge_d","X22.2_HomeRange_Indiv_km2","X13.2_NeonateHeadBodyLen_mm","X10.1_PopulationGrpSize","X13.3_WeaningHeadBodyLen_mm")
pantheria = pantheria[keeps]

## rename column names on which to match for merging with COMBINE database
names(pantheria)[1:5] <- c("order","family","genus","species","pan_binomial")  

## merge COMBINE and PanTHERIA databases
  # traits <- full_join(combine, pantheria, by = c("order","family","genus","species"),keep=TRUE)
  # ~1300 names do not match by O/F/G/S between pantheria and combine
traits <- full_join(combine, pantheria, by = c("genus","species"),keep=TRUE)
  # ~500 names do not match by G/S between pantheria and combine

## export PanTHERIA names that do not match with COMBINE names
pan_mismatch <- traits[is.na(traits$order.x),]
pan_mismatch <- pan_mismatch[order(pan_mismatch$order.y, pan_mismatch$family.y, pan_mismatch$genus.y, pan_mismatch$species.y),]
write.csv(pan_mismatch, 'data/names/pan_mismatch.csv')

## For now, let's work with only observations that matched in COMBINE
traits <- subset(traits, (!is.na(traits$order.x)), select=-c(order.y,genus.y,family.y,species.y,pan_binomial))
rm(combine, pantheria, pan_mismatch, keeps)

