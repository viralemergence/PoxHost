## PoxHost 03: mammal orthopoxvirus brt
## katie.tseng@wsu.edu
## updated 03/21/2022

## clean environment & plots
rm(list=ls()) 
graphics.off()

## libraries
library(gbm)
library(fastDummies)
library(rsample)
library(ROCR)
library(sciplot)
library(ggplot2)
library(pdp)
library(PresenceAbsence)
library(tidyr)
library(viridis)
library(caper)
library(phylofactor)
library(ggtree)
library(treeio)
library(caret) 
library(InformationValue)
library(mgcv)

## set working directory and load files
setwd("/Users/CDDEP/Downloads/OPV Host Prediction/PoxHost")
data=read.csv('data/cleaned/pox cleaned response and traits.csv')

## classify true negatives
data$type=ifelse(data$pcr==0 & data$competence==0,"true negative","other")

## which species is competent but no PCR record?
set=data
set$treename[set$pcr==0 & set$competence==1]

## tabulate PCR/infection and isolation
set$inf=ifelse(set$pcr==0,"PCR negative","PCR positive")
set$iso=ifelse(set$competence==0,"no isolation","isolation")
table(set$inf,set$iso)

## make binary columns for genus
dums=dummy_cols(data["gen"])

## unique
dums=dums[!duplicated(dums$gen),]

## ensure all factor
for(i in 1:ncol(dums)){
  
  ## column as factor
  dums[,i]=factor(dums[,i])
}

## merge
data=merge(data,dums,by="gen",all.x=T)
rm(dums)

## drop unnecessary columns
data$X=NULL
data$traitname=NULL

## mode function
mode.prop <- function(x) {                
  ux <- unique(x[is.na(x)==FALSE])        # creates array of unique values
  tab <- tabulate(match(na.omit(x), ux))  # creates array of the frequency (number of times) a unique value appears in a column 
  max(tab)/length(x[is.na(x)==FALSE])     # max-frequency / number of elements in each column that are not NA
}

## assess variation across columns (2 indicates columns)
vars=data.frame(apply(data,2,function(x) mode.prop(x)),
                apply(data,2,function(x) length(unique(x))))    # number of unique elements in each column

## get names
vars$variables=rownames(vars)
names(vars)=c("var","uniq","column")

## round values
vars$var=round(vars$var,2)

## label variables "cut" if homogeneous (100%)
vars$keep=ifelse(vars$var<1,"keep","cut")
vars$keep=ifelse(vars$column%in%c('pcr','competence','antibodies','fam'),'keep',vars$keep) # ensures we keep these columns
vars=vars[order(vars$keep),]

## trim (creates array of column names to cut and removes from df)
keeps=vars[-which(vars$keep=="cut"),]$column

## drop if no variation
data=data[keeps]
rm(keeps,vars)

## assess missing values
mval=data.frame(apply(data,2,function(x) length(x[!is.na(x)])/nrow(data))) # proportion of values that are not NA

## get names
mval$variables=rownames(mval)
names(mval)=c("comp","column")

## visualize distribution of NA
png("figs/Figure S1.png", width=4,height=4,units="in",res=600)
ggplot(mval[!mval$column%in%c("gen","treename","pcr","competence","tip.label","fam"),],
       aes(comp))+
  geom_histogram(bins=50)+
  geom_vline(xintercept=0.70,linetype=2,size=0.5)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  labs(y="frequency",
       x="trait coverage across mammal species (genus)")+
  scale_x_continuous(labels = scales::percent)
dev.off()

## round values
mval$comp=round(mval$comp,2)

## label variables "cut" if >30% values are NA
mval$keep=ifelse(mval$comp>=0.70,"keep","cut")
table(mval$keep)
mval=mval[order(mval$keep),]

## trim (creates array of column names to cut and removes from df)
keeps=mval[-which(mval$keep=="cut"),]$column

## order
mval=mval[order(mval$comp),]

## drop if not well represented
data=data[keeps]
rm(keeps,mval)

## simplify
set=data
set$treename=NULL
set$tip.label=NULL
set$type=NULL
set$fam=NULL
set$ord=NULL
set$clade=NULL

## covariates
ncol(set)-3

## coverage table s1
ts1=data.frame(apply(set,2,function(x) length(x[!is.na(x)])/nrow(set)))

## get names
ts1$variables=rownames(ts1)
names(ts1)=c("comp","column")

## trim disease data
ts1=ts1[!ts1$column%in%c("pcr","competence","antibodies"),]
ts1$feature=ts1$column
ts1$column=NULL
ts1$coverage=ts1$comp
ts1$comp=NULL

## write file
write.csv(ts1, "figs/TableS1.csv")

## hyperparameter tuning ifelse
hok="ok"
if(hok!="ok"){
  
  ## hyperparameter grid (creates df from all combos of supplied vectors or factors)
  hgrid=expand.grid(n.trees=5000,                   # number of trees or iterations (increasing trees reduces error on training set, but set too high may lead to over-fitting)
                    interaction.depth=c(2,3,4),     # number of splits it has to perform on a tree (N; max nodes per tree), where each split increases total nodes by 3 (3*N+1) and terminal nodes by 2 (2*N+1)
                    shrinkage=c(0.01,0.001,0.0005), # reduces/shrinks the impact of each additional fitted tree (aka learning rate); generally gives better result, but at the expense of more trees/iterations required
                    n.minobsinnode=4,               # min number of obs in trees' terminal nodes (splitting of nodes ceases when 4 obs are in each terminal node); when training samples are small, may need to lower this setting
                    seed=seq(1,10,by=1))
              
  ## fix trees
  hgrid$n.trees=ifelse(hgrid$shrinkage<0.001,hgrid$n.trees*3,hgrid$n.trees)
  
}

#Katie: before continuing, make sure binary variables are numeric and not factor