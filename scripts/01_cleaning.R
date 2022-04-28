## PoxHost 01: mammal orthopoxvirus (OPV) cleaning
## katie.tseng@wsu.edu
## updated 2/15/2022

## clean environment & plots
rm(list=ls()) 
graphics.off()

## libraries
library(ape)
library(dplyr)
library(nlme)
library(tidyverse)
library(vroom)
library(treespace) #dependency XQuartz v2.7.11 (https://www.xquartz.org/releases/XQuartz-2.7.11.html); for help installing, see https://stackoverflow.com/questions/66011929/package-rgl-in-r-not-loading-in-mac-os/66127391#66127391

## set working directory
setwd("~/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Fernandez Lab/Projects (Active)/OPV Host Prediction/GitHub/PoxHost/")


######################################### POXVIRUS DATA #########################################

## load in orthopoxvirus data (source: https://github.com/viralemergence/virion/tree/main/Virion)
virion <- vroom("/Users/katietseng/Fernandez Lab Dropbox/Katie Tseng/Mac/Desktop/PoxHost(copy)/data/raw/Virion.csv.gz")

## view detection methods
virion %>%
  filter(VirusGenus == "orthopoxvirus") %>%
  pull(DetectionMethod) %>%
  unique()

## extract associations detected by PCR & Isolation only
poxdata <- virion %>% 
  filter(VirusGenus == "orthopoxvirus" & (DetectionMethod %in% c("PCR/Sequencing","Isolation/Observation","Antibodies"))) %>% 
  dplyr::select(Host,Virus,DetectionMethod) %>%
  unique()

## exclude NAs and variola virus
poxdata <- poxdata[!is.na(poxdata$Host),]
poxdata <- poxdata[!is.na(poxdata$Virus),]
poxdata <- poxdata[!(poxdata$Virus=="variola virus"),]

## clean
#poxdata <- subset(poxdata, select=-c(Virus))

## subset detection by Antibodies
antibodies <- poxdata[which(poxdata$DetectionMethod=="Antibodies"),]
antibodies$antibodies <- 1

## subset detection by PCR
pcr <- poxdata[which(poxdata$DetectionMethod=="PCR/Sequencing"),]
pcr$pcr <- 1

## subset detection by Isolation
competence <- poxdata[which(poxdata$DetectionMethod=="Isolation/Observation"),]
competence$competence <- 1

## merge PCR & Isolation data
poxdata <- merge(pcr, competence, by=c("Host","Virus"), all=TRUE)
poxdata <- merge(poxdata, antibodies, by=c("Host","Virus"), all=TRUE)

## clean dataset
poxdata <- subset(poxdata, select=-c(DetectionMethod,DetectionMethod.x,DetectionMethod.y))

## recode as binary (replace NA with zeroes)
poxdata$antibodies=ifelse(is.na(poxdata$antibodies),0,poxdata$antibodies)
poxdata$pcr=ifelse(is.na(poxdata$pcr),0,poxdata$pcr)
poxdata$competence=ifelse(is.na(poxdata$competence),0,poxdata$competence)

## create genus variable
poxdata$Host=(str_replace(poxdata$Host," ","_"))
poxdata$Host=str_to_title(poxdata$Host)
poxdata$gen=sapply(strsplit(as.character(str_to_title(poxdata$Host)),"_"),function(x) x[1])

# # examine pcr- or competence-only positive associations by genus and virus
 temp <- subset(poxdata, pcr==1 | competence==1)
 temp %>% count(gen, sort = TRUE)
 temp %>% count(Virus, sort = TRUE)

## clean environment
rm(virion,antibodies,pcr,competence,temp)

######################################### MERGE POXVIRUS & TAXONOMY DATA######################################### 

## load in taxonomy (source: https://data.vertlife.org/)
taxa=read.csv('phylo/taxonomy_mamPhy_5911species.csv',header=T) 

## check that poxdata genus names match those in gdata
poxdata$gen[!poxdata$gen %in% taxa$gen]

## rename VertLife taxonomy for "Procolobus badius" to "Piliocolobus badius" per IUCN taxonomy (https://www.iucnredlist.org/species/161247840/161259430)
taxa$gen[taxa$gen=="Procolobus"] <- "Piliocolobus"

## simplify to genus
gtaxa=taxa[!duplicated(taxa$gen),] #note: drops duplicates at the genus level
gtaxa=gtaxa[c('gen','fam','ord','clade','higher')]

## merge pox data & taxonomy data (left-join)
data=merge(gtaxa,poxdata,by='gen',all.x=TRUE)

## sort by taxonomy level (high to low)
data <- data[order(data$higher,data$clade,data$ord,data$fam,data$gen),]

## note genera that are unsampled/pseudoabsences
data$studied=ifelse(is.na(data$pcr) & is.na(data$competence),0,1)

## create pseudo-absences for viral detection
data$pcr=ifelse(is.na(data$pcr),0,data$pcr)
data$competence=ifelse(is.na(data$competence),0,data$competence)
data$antibodies=ifelse(is.na(data$antibodies),0,data$antibodies)

## subset data of pcr/isolation positive associations only
keep=subset(data, pcr==1 | competence==1)

## create variable denoting whether the order of genera with positive associations exists in merged data 
data$keep=ifelse(data$ord %in% keep$ord,TRUE,FALSE)

## subset true values
data=subset(data,keep==TRUE)

## clean dataset
data=subset(data, select=-c(antibodies,keep))

## clean environment
rm(poxdata,taxa,gtaxa,keep)


######################################### LOAD IN TRAIT DATA #########################################

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

## clean environment
rm(combine, pantheria, pan_mismatch, keeps)


######################################### LOAD IN PHYLOGENY DATA #########################################

## load Upham phylogeny 
# Source: https://doi.org/10.5061/dryad.tb03d03 (in subfolder "Data_S8_finalFigureFiles" > "_DATA")
tree=read.nexus('phylo/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')

## simplify species names
tree$tip.label[tree$tip.label=="_Anolis_carolinensis"] <- "Anolis_carolinensis"
tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep='_'))
tree$gtip=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],sep='_'))


######################################### CHECK FOR MISMATCHED NAMES #########################################

## are all poxdata in tree?
data$gtip=data$gen
data$intree=ifelse(data$gtip%in%setdiff(data$gtip,tree$gtip),'missing','upham')

## are all poxdata in traits?
traits$gtip=traits$genus.x
data$intraits=ifelse(data$gtip%in%setdiff(data$gtip,traits$gtip),'missing','traits')

## create df of just obs w/ missing names
miss=data[c('Host','gtip','intree','intraits')]
miss=miss[miss$intree=='missing'|miss$intraits=='missing',]
miss=miss[order(miss$intree,miss$intraits),]

## export missing names: 
write.csv(miss,'data/names/poxhost_name_mismatch.csv')
rm(miss)

## to species name 
  # tree$tip.label[which(grepl('Colobus',tree$tip.label))]
  # traits$iucn2020_binomial[which(grepl('Piliocolobus',traits$iucn2020_binomial))]
  # to identify related species see IUCN (https://www.iucnredlist.org/) and NCBI (https://www.ncbi.nlm.nih.gov/taxonomy)

## load in revised treename for "Piliocolobus badius": see NCBI Taxonomy Browser for homotypic synonym
fix=read.csv('data/names/poxhost_name_mismatch_KTedits.csv',header=TRUE)
fix=fix[c('gtip','treename','traitname','proxy')]

## merge in revised names with poxdata (left-join) 
data=merge(data,fix,by='gtip',all.x=T)

## if blank, NA
data$traitname=ifelse(data$traitname=='',NA,as.character(data$traitname))
data$treename=ifelse(data$treename=='',NA,as.character(data$treename))

## if NA, tip
data$treename=ifelse(is.na(data$treename),as.character(data$gtip),as.character(data$treename))
data$traitname=ifelse(is.na(data$traitname),as.character(data$gtip),as.character(data$traitname))

## trait name includes proxy
data$traitname=ifelse(data$intraits=='missing' & is.na(data$traitname),as.character(data$proxy),
                      ifelse(data$intraits=='missing' & !is.na(data$traitname),as.character(data$traitname),
                             as.character(data$gtip)))

## simplify
data=data[c('gtip','gen','fam','ord','clade','treename','traitname','pcr','competence','studied')]
rm(fix)


######################################### AGGREGATE DATA AT GENUS LEVEL #########################################

#### POX DATA ####

## check/fix duplicates in phylogeny
set=data[c('pcr','competence','treename')]

## which treenames occur multiple times
unique(set$treename)[table(set$treename)>1]

## aggregate
set=aggregate(.~treename,set,sum)

## merge with meta
set2=merge(set,data[!duplicated(data$treename),c('studied','gtip','gen','fam','ord','clade','treename','traitname')],all.x=T,by='treename')

## simplify
gdata=set2
rm(data,set,set2)

## fix binomial
gdata$pcr=ifelse(gdata$pcr>0,1,0)
gdata$competence=ifelse(gdata$competence>0,1,0)

#### TRAIT DATA ####

## to aggregate continuous/integer variables, use the median as the summary measure
colnames(traits)
gtraits_continuous=aggregate(cbind(adult_mass_g,brain_mass_g,adult_body_length_mm,adult_forearm_length_mm,
                                   max_longevity_d,maturity_d,female_maturity_d,male_maturity_d,
                                   age_first_reproduction_d,gestation_length_d,teat_number_n,
                                   litter_size_n,litters_per_year_n,interbirth_interval_d,
                                   neonate_mass_g,weaning_age_d,weaning_mass_g,generation_length_d,
                                   dispersal_km,density_n_km2,home_range_km2,social_group_n,
                                   dphy_invertebrate,dphy_vertebrate,dphy_plant,
                                   det_inv,det_vend,det_vect,det_vfish,det_vunk,det_scav,det_fruit,det_nect,det_seed,det_plantother,det_diet_breadth_n,
                                   upper_elevation_m,lower_elevation_m,altitude_breadth_m,habitat_breadth_n,
                                   X2.1_AgeatEyeOpening_d,X18.1_BasalMetRate_mLO2hr,X5.2_BasalMetRateMass_g,
                                   X7.1_DispersalAge_d,X22.2_HomeRange_Indiv_km2,
                                   X13.2_NeonateHeadBodyLen_mm,X10.1_PopulationGrpSize,X13.3_WeaningHeadBodyLen_mm) 
                             ~ order.x+family.x+genus.x+gtip, data=traits, FUN=median, na.action=na.pass, na.rm=TRUE)
  # 'na.action=na.pass, na.rm=TRUE' is specified such that if species w/in a genus has a combination of real values & NAs, the median of real values will be returned (as opposed to omitting the genus or returning NA)

## to aggregate binary variables, use the mean as the summary measure
traits$fossoriality[traits$fossoriality==2]<-0  #recode 0/1
gtraits_binary=aggregate(cbind(hibernation_torpor,fossoriality,freshwater,marine,terrestrial_non.volant,terrestrial_volant,
                               island_dwelling,disected_by_mountains,glaciation) ~ order.x+family.x+genus.x+gtip, data=traits, FUN=mean, na.action=na.pass, na.rm=TRUE)

## to aggregate categorical variables transform into binary
gtraits_cat <- traits
gtraits_cat$trophic_herbivores <- ifelse(gtraits_cat$trophic_level==1,1,0)
gtraits_cat$trophic_omnivores <- ifelse(gtraits_cat$trophic_level==2,1,0)
gtraits_cat$trophic_carnivores <- ifelse(gtraits_cat$trophic_level==3,1,0)
gtraits_cat$activity_nocturnal <- ifelse(gtraits_cat$activity_cycle==1,1,0)
gtraits_cat$activity_crepuscular <- ifelse(gtraits_cat$activity_cycle==2,1,0) #nocturnal/crepuscular, cathemeral, crepuscular or diurnal/crepuscular
gtraits_cat$activity_diurnal <- ifelse(gtraits_cat$activity_cycle==3,1,0)
gtraits_cat$forager_marine <- ifelse(gtraits_cat$foraging_stratum=="M",1,0)
gtraits_cat$forager_ground <- ifelse(gtraits_cat$foraging_stratum=="G",1,0) 
gtraits_cat$forager_scansorial <- ifelse(gtraits_cat$foraging_stratum=="S",1,0)
gtraits_cat$forager_arboreal <- ifelse(gtraits_cat$foraging_stratum=="Ar",1,0)
gtraits_cat$forager_aerial <- ifelse(gtraits_cat$foraging_stratum=="A",1,0)
gtraits_cat$island_end_marine <- ifelse(gtraits_cat$island_endemicity=="Exclusively marine",1,0)
gtraits_cat$island_end_mainland <- ifelse(gtraits_cat$island_endemicity=="Occurs on mainland",1,0)
gtraits_cat$island_end_lgbridge <- ifelse(gtraits_cat$island_endemicity=="Occurs on large land bridge islands",1,0)
#gtraits_cat$island_end_smbridge <- ifelse(gtraits_cat$island_endemicity=="Occurs on small land bridge islands",1,0)
gtraits_cat$island_end_isolated <- ifelse(gtraits_cat$island_endemicity=="Occurs only on isolated islands",1,0)
gtraits_cat$biogeo_afrotropical <- ifelse(grepl("Afrotropical",gtraits_cat$biogeographical_realm),1,0)
gtraits_cat$biogeo_antarctic <- ifelse(grepl("Antarctic",gtraits_cat$biogeographical_realm),1,0)
gtraits_cat$biogeo_australasian <- ifelse(grepl("Australasian",gtraits_cat$biogeographical_realm),1,0)
gtraits_cat$biogeo_indomalayan <- ifelse(grepl("Indomalayan",gtraits_cat$biogeographical_realm),1,0)
gtraits_cat$biogeo_nearctic <- ifelse(grepl("Nearctic",gtraits_cat$biogeographical_realm),1,0)
gtraits_cat$biogeo_neotropical <- ifelse(grepl("Neotropical",gtraits_cat$biogeographical_realm),1,0)
gtraits_cat$biogeo_oceanian <- ifelse(grepl("Oceanian",gtraits_cat$biogeographical_realm),1,0)
gtraits_cat$biogeo_palearctic <- ifelse(grepl("Palearctic",gtraits_cat$biogeographical_realm),1,0)

## to aggregate transformed categorical (binary) variables, use the mean as the summary measure
gtraits_cat=aggregate(cbind(trophic_herbivores,trophic_omnivores,trophic_carnivores,
                            activity_nocturnal,activity_crepuscular,activity_diurnal,
                            forager_marine,forager_ground,forager_scansorial,forager_arboreal,forager_aerial,
                            island_end_marine,island_end_mainland,island_end_lgbridge,island_end_isolated,
                            biogeo_afrotropical,biogeo_antarctic,biogeo_australasian,biogeo_indomalayan,biogeo_nearctic,biogeo_neotropical,biogeo_oceanian,biogeo_palearctic)
                      ~ order.x+family.x+genus.x+gtip, data=gtraits_cat, FUN=mean, na.action=na.pass, na.rm=TRUE)

## merge continuous variables with binary variables and clean environment
gtraits <- full_join(gtraits_continuous, gtraits_binary, by = c("order.x","family.x","genus.x","gtip"),keep=TRUE)
colnames(gtraits)[1:4] <- c("order.x","family.x","genus.x","gtip")
gtraits=subset(gtraits, select=-c(order.x.y,family.x.y,genus.x.y,gtip.y))

## merge transformed categorical variables and clean environment
gtraits <- full_join(gtraits, gtraits_cat, by = c("order.x","family.x","genus.x","gtip"),keep=TRUE)
colnames(gtraits)[1:4] <- c("ord","fam","gen","gtip")
gtraits=subset(gtraits, select=-c(order.x.y,family.x.y,genus.x.y,gtip.y))

## clean 
rm(traits,gtraits_binary,gtraits_cat,gtraits_continuous)

#### TREE DATA ####

## create df linking tip labels with their corresponding categories (genus and species)
tdata <- data.frame(matrix(NA,nrow=length(tree$tip.label),ncol=0))
tdata$genus=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],sep='_'))
tdata$species=tree$tip.label

## collapse tree to genus level and clean environment
gtree=makeCollapsedTree(tree=tree,df=tdata[c('genus','species')])
rm(tdata,tree)


######################################### MERGE POX DATA, TRAITS, & TREE #########################################

## merge traits
gtraits$traitname=gtraits$gtip
gtraits$trait=1
gdata=merge(gdata,gtraits,by=c('traitname'),all.x=T)
rm(gtraits)

## simplify
gdata <- subset(gdata, select=-c(gen.y,ord.y,fam.y,gtip.y,trait))
gdata <- gdata %>% rename(gtip=gtip.x, gen=gen.x, fam=fam.x, ord=ord.x)

## trim tree
gtree=keep.tip(gtree,gtree$tip.label[gtree$tip.label%in%gdata$treename])

## fix
gtree$gtip=NULL
gtree=makeLabel(gtree)


######################################### ADD PUBMED CITATIONS MEASURE #########################################

## pubmed citations
library(easyPubMed)

## function
counter=function(name){
  as.numeric(as.character(get_pubmed_ids(gsub('_','-',name))$Count))
}
citations=c()

## loop through
for(i in 1:length(gdata$treename)) {
  citations[i]=counter(gdata$treename[i])
  print(i)
}

## compile
cites=data.frame(treename=gdata$treename,
                 cites=citations)

## merge
gdata=merge(gdata,cites,by='treename')

## clean
rm(cites,citations,i,counter)


######################################### ADD EVOLUTIONARY DISTINCTIVENESS MEASURE #########################################

## get ed
library(picante) #before loading picante, make sure latest version of nlme package is loaded
ed=evol.distinct(gtree,type='equal.splits')
  # note: calculates evolutionary distinctiveness measures for a suite of species by equal splits and fair proportions; returns species score

## treename
ed$treename=ed$Species
ed$Species=NULL

## rename
ed$ed_equal=ed$w
ed$w=NULL

## merge into data
gdata=merge(gdata,ed,by='treename',all.x=T)
rm(ed)

## export files
write.csv(gdata,'data/cleaned/pox cleaned response and traits.csv')
saveRDS(gtree,'data/cleaned/mammal phylo trim.rds')
