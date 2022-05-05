## PoxHost 05: maps of predicted probabilities
## katie.tseng@wsu.edu
## updated 04/27/2022


## clean environment & plots
rm(list=ls()) 
graphics.off()

## libraries
library(classInt)
library(tidyverse)
library(raster)
library(rgdal)
library(dismo)
library(XML)
library(maps)
library(sp)
library(dplyr)

library(devtools)
install_github("hunzikp/velox")
library(velox)

## load files
setwd("/Users/katietseng/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Fernandez Lab/Projects (Active)/OPV Host Prediction/GitHub/PoxHost")
set.seed(12345)
pred <- read.csv("data/cleaned/PoxHost predictions.csv")

# 1. Threshold the results
library(PresenceAbsence)

###################### INTERMISSION: THRESHOLD IMPACTS

ts.p <- optimal.thresholds(data.frame(pred[,c('treename','PCR','pred_pcr')]),
                           threshold = 10001,
                           opt.methods = c(2,4,5,10),
                           req.sens = 0.90,
                           na.rm = TRUE)

cut.p <- function(x) {sum(pred$pred_pcr[pred$PCR==0] > x)}

sapply(unlist(ts.p[2]), cut.p)


ts.c <- optimal.thresholds(data.frame(pred[,c('treename','competence','pred_comp')]),
                           threshold = 10001,
                           opt.methods = c(2,4,5,10),
                           req.sens = 0.90,
                           na.rm = TRUE)

cut.c <- function(x) {sum(pred$pred_comp[pred$competence==0] > x)}

sapply(unlist(ts.c[2]), cut.c)




###################### MOVE FORWARD

t.pcr <- optimal.thresholds(data.frame(pred[,c('treename','PCR','pred_pcr')]),
                            threshold = 10001,
                            opt.methods = 10,
                            req.sens = 0.95,
                            na.rm = TRUE)

t.comp <- optimal.thresholds(data.frame(pred[,c('treename','competence','pred_comp')]),
                             threshold = 10001,
                             opt.methods = 10,
                             req.sens = 0.95,
                             na.rm = TRUE)

# Threshold the results to binary outputs

pred %>%
  mutate(bin_comp = (pred_comp > t.comp$pred_comp),
         bin_pcr = (pred_pcr) > t.pcr$pred_pcr) -> pred

# How many predicted undiscovered hosts by PCR?

table(pred$pred_pcr[pred$PCR==0] > t.pcr$pred_pcr)

# How many predicted undiscovered hosts by competence

table(pred$pred_comp[pred$competence==0] > t.comp$pred_comp)

# Do the predicted competence hosts overlap with PCR

pred %>% filter(competence==0, pred_comp > t.comp$pred_comp)

# Pull out the relevant lists 

pred %>% filter(competence==1) %>% pull(treename) %>% gsub("_"," ",.) -> known.comp

pred %>% filter(PCR==1) %>% pull(treename) %>% gsub("_"," ",.) -> known.pcr

pred %>% filter(bin_comp==1) %>% pull(treename) %>% gsub("_"," ",.) -> pred.comp

pred %>% filter(bin_pcr==1) %>% pull(treename) %>% gsub("_"," ",.) -> pred.pcr

sort(pred.pcr[!(pred.pcr %in% known.pcr)])

sort(pred.comp[!(pred.comp %in% known.comp)])

# 2. Let's make some maps?

library(fasterize)
library(rgdal)   # switch to sf in 2023
library(raster)
library(sf)

#iucn <- st_read(dsn = "C:/Users/cjcar/Dropbox/HowManyHelminths2019", layer='TERRESTRIAL_MAMMALS')
iucn <- st_read(dsn = "/Users/katietseng/Fernandez Lab Dropbox/Katie Tseng/Mac/Desktop/PoxHost(copy)/data/raw/MAMMALS/MAMMALS.shp", layer='MAMMALS')

r <- disaggregate(getData("worldclim",var="alt",res=2.5)*0,2) # Make a blank raster

# Create four layers
iucn$treename=sapply(strsplit(iucn$binomial,' '),function(x) paste(x[1],sep=' '))

iucn.1 <- iucn[iucn$treename %in% known.comp,] 
iucn.2 <- iucn[iucn$treename %in% known.pcr,] 
iucn.3 <- iucn[iucn$treename %in% pred.comp,] 
iucn.4 <- iucn[iucn$treename %in% pred.pcr,] 


map.knc <- (fasterize(iucn.1, r, fun="sum"))
map.knp <- (fasterize(iucn.2, r, fun="sum"))
map.prc <- (fasterize(iucn.3, r, fun="sum"))
map.prp <- (fasterize(iucn.4, r, fun="sum"))

fix <- function(x) {sum(x,r,na.rm=TRUE)+r} # This adds zeros for the continental area

map.knc <- fix(map.knc)
map.knp <- fix(map.knp)
map.prc <- fix(map.prc)
map.prp <- fix(map.prp)

raster::stack(map.knp, map.knc, map.prp, map.prc) %>%    #tera package
  crop(c(-170,-25,-90,90)) %>% # sf::intersection
  raster::trim() -> maps

names(maps) <- c('KnownPCR', 'KnownComp', 'PredPCR', 'PredComp')

# Generate the actual visualization

library(rasterVis)
library(RColorBrewer)

mycolors <- colorRampPalette(rev(brewer.pal(10,"Spectral")))(21)
mycolors[1] <- "#C0C0C0"


setwd("/Users/katietseng/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Fernandez Lab/Projects (Active)/OPV Host Prediction/GitHub/PoxHost/figs")
png("Figure 4.png",width=10,height=10,units="in",res=300)
rasterVis::levelplot(maps,  
                     col.regions = mycolors,
                     #at = seq(0, 15, 1),
                     alpha = 0.5, 
                     scales=list(alternating=FALSE),
                     par.strip.text=list(cex=0),
                     xlab = NULL, ylab = NULL,
                     #labels = labels,
                     maxpixels = 5e6)
# 
# ggarrange(maps.all,ncol=1,nrow=1,
#           labels=c("(A)",paste("n=",length(known.pcr)),
#                    "(B)",paste("n=",length(known.comp)),
#                    "(C)",paste("n=",length(pred.pcr)),
#                    "(D)",paste("n=",length(pred.comp))),
#           label.x=c(0.05,0.05,0.55,0.55,0.05,0.55,0.05,0.55),    #defaults to 0 - left
#           label.y=c(0.95,0.55,0.95,0.55,0.45,0.05,0.45,0.05),       #defaults to 1 - top
#           font.label=list(face="plain",size=12))

dev.off()
