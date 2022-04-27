## PoxHost 02: mammal orthopoxvirus phylofactor
## katie.tseng@wsu.edu
## updated 04/08/2022

## clean environment & plots
rm(list=ls()) 
graphics.off()

## libraries (for help installing phylofactor: https://reptalex.github.io/phylofactor/)
library(ape)
library(caper)
library(data.table)
library(BiocManager)  ## BiocManager::install(c("Biostrings","ggtree"))
library(phylofactor)  ## devtools::install_github('reptalex/phylofactor')
library(treeio)       ## BiocManager::install("treeio")

## load files
setwd("~/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Fernandez Lab/Projects (Active)/OPV Host Prediction/GitHub/PoxHost/")
data=read.csv('data/cleaned/pox cleaned response and traits.csv')
mtree=readRDS('data/cleaned/mammal phylo trim.rds')

## create pcr+comp var
data$compcr=ifelse(data$pcr==1|data$competence==1,1,0)

## create var labeling observations to keep if genus name in data is also in mtree
data$tree=ifelse(data$treename%in%setdiff(data$treename,mtree$tip.label),'cut','keep')

## trim
bdata=data[which(data$tree=='keep'),]

## subsets bdata based on positions of matches (returned as vector) - redundant?
bdata=bdata[match(mtree$tip.label,bdata$treename),]

## save
bdata$label=bdata$treename
bdata$Species=bdata$treename  #We refer to treenames of genus as "Species" as this is required in later functions

## merge: caper::comparative.data combines phylogenies w/ datasets and ensures consistent structure and ordering
cdata=comparative.data(phy=mtree,data=bdata,names.col=treename,vcv=T,na.omit=F,warn.dropped=T)

## fix
cdata$data$tree=NULL

## proportion detected for each detection method
nrow(data)
round(prop.table(table(data$pcr)),4)*100        #values in each cell are divided by the sum of the 4 cells  
round(prop.table(table(data$competence)),4)*100
round(prop.table(table(data$compcr)),4)*100

## phylogenetic signal in response
## D of 0 = Brownian model, D of 1 = random (no phylogenetic signal)
set.seed(1)
mod1=phylo.d(cdata,binvar=pcr,permut=10000); mod1
set.seed(1)
mod2=phylo.d(cdata,binvar=competence,permut=10000); mod2

## create taxonomy var 
cdata$data$taxonomy=paste(cdata$data$ord,cdata$data$fam,cdata$data$gen,sep='; ')
#cdata$data$taxonomy=paste(cdata$data$fam,cdata$data$gen,sep='; ')

## create data frame of taxonomy
taxonomy=data.frame(cdata$data$taxonomy)
names(taxonomy)="taxonomy"
taxonomy$Species=rownames(cdata$data)
taxonomy=taxonomy[c("Species","taxonomy")]
taxonomy$taxonomy=as.character(taxonomy$taxonomy)

## Holm rejection procedure                      
HolmProcedure <- function(pf,FWER=0.05){         #creates function HolmProcedure with arguments pf and FWER=0.05; FWER=family-wise error rate (alpha level set to 0.05)
  
  ## get split variable
  cs=names(coef(pf$models[[1]]))[-1]             #sets 'cs' as the names of the model coefficients (i.e., variable name) extracted by 'coef' in the 1st list element of 'pf$models' minus the 1st element among those names; double brackets access a list element; coef extracts model coefficients
  split=ifelse(length(cs)>1,cs[3],cs[1])         #returns the 3rd element in 'cs' if the length of the number of elements in 'cs' is greater than 1, or else returns the 1st element
  
  ## obtain p values
  if (pf$models[[1]]$family$family%in%c('gaussian',"Gamma","quasipoisson")){                  #if family$family of the 1st list element of pf$models is in the columns 'gaussian', etc...
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|t|)'])  #then, to each element of pf$models, we apply the summary function with the argument 'fit' and assign the output to 'pvals'; specifically, we use 'summary(fit)' to call the output of 'pf$models', extracting the 'coefficients' section, whereby we index the column named 'Pr(>|t\)' and split the data in that column; see sample output of linear model of R for reference (https://feliperego.github.io/blog/2015/10/23/Interpreting-Model-Output-In-R) 
  } else {
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|z|)'])  #or else, extract p-val based on z statistic
  }
  D <- length(pf$tree$tip.label)                                                              #returns number of elements in pf$tree$tip.label
  
  ## this is the line for Holm's sequentially rejective cutoff                                #HB = Target alpha / (n – rank + 1)
  keepers <- pvals<=(FWER/(2*D-3 - 2*(0:(pf$nfactors-1))))                                    #returns TRUE/FALSE if p-values are <= to 0.05/(n-rank+1)
  
  if (!all(keepers)){                            #if not all pvals were keepers (i.e., all items in keepers were true)...
    nfactors <- min(which(!keepers))-1           #then assign nfactors to the minimum/earliest position of items in keepers that were false, minus 1.
  } else {
    nfactors <- pf$nfactors                      #or else, assign nfactors as the value of pf$nfactors
  }
  return(nfactors)
}

## get species in a clade
cladeget=function(pf,factor){                        #creates function 'cladeget' w/ arguments 'pf' and 'factor'
  spp=pf$tree$tip.label[pf$groups[[factor]][[1]]]    #returns n'th element of the pf$tree$tip.label based on the value of the first component inside the n'th ('factor') component of  'pf$groups'
  return(spp)
}

## summarize pf object                               
pfsum=function(pf){                                  #creates function 'pfsum' w/ argument 'pf'
  
  ## get formula
  chars=as.character(pf$frmla.phylo)[-1]             #returns pf$frmla.phylo minus the first element
  
  ## response                  
  resp=chars[1]                                      #returns 1st element of chars              
  
  ## holm
  hp=HolmProcedure(pf)                               #runs HolmProcedure function using argument pf and assigns to hp
  
  ## save model
  model=chars[2]                                     #returns 2nd element of chars 
  
  ## set key
  setkey(pf$Data,'Species')                          #creates key on the sorted column 'Species' in the datatable pf$Data
  
  ## make data
  dat=data.frame(pf$Data)                            
  
  ## make clade columns in data
  for(i in 1:hp){
    
    dat[,paste0(resp,'_pf',i)]=ifelse(dat$Species%in%cladeget(pf,i),'factor','other')   #paste0 concatenates all elements w/o a separator
    
  }
  
  ## make data frame to store taxa name, response, mean, and other
  results=data.frame(matrix(ncol=6, nrow = hp))                        #creates empty df with 6 columns and hp-number of rows
  colnames(results)=c('factor','taxa','tips','node',"clade",'other')
  
  ## set taxonomy
  taxonomy=dat[c('Species','taxonomy')]                                #creates dataset of Species and taxonomy
  taxonomy$taxonomy=as.character(taxonomy$taxonomy)
  
  ## loop
  for(i in 1:hp){
    
    ## get taxa
    tx=pf.taxa(pf,taxonomy,factor=i)$group1                  #gets taxonomic order
    
    ## get tail
    tx=sapply(strsplit(tx,'; '),function(x) tail(x,1))       #gets taxonomic family as list
    
    ## combine
    tx=paste(tx,collapse=', ')                               #collapses taxonomic family into single string
    
    # save
    results[i,'factor']=i                                    #returns index number in 'factor' column
    results[i,'taxa']=tx                                     #returns string element (tx) in 'taxa' column
    
    ## get node
    tips=cladeget(pf,i)
    node=ggtree::MRCA(pf$tree,tips)                          #MRCA = finds Most Recent Common Ancestor among a vector of tips 
    results[i,'tips']=length(tips)
    results[i,'node']=ifelse(is.null(node) & length(tips)==1,'species',
                             ifelse(is.null(node) & length(tips)!=1,NA,node))
    
    ## get means
    ms=(tapply(dat[,resp],dat[,paste0(resp,'_pf',i)],FUN=mean))   #tapply takes mean of '1 vs. 0' (dat[,resp]) by 'other'/'factor' type (dat[,paste...]
    
    ## add in
    results[i,'clade']=ms['factor']
    results[i,'other']=ms['other']
    
  }
  
  ## return
  return(list(set=dat,results=results))       #returns number of clades with significantly greater propensity of infection adjusting for FWER using Holm rejection procedure
}

## PCR
set.seed(1)
pcr_pf=gpf(Data=cdata$data,tree=cdata$phy,
           frmla.phylo=pcr~phylo,
           family=binomial,algorithm='phylo',nfactors=10,min.group.size=5)

## summarize
HolmProcedure(pcr_pf)
pcr_pf_results=pfsum(pcr_pf)$results        

## competence
set.seed(1)
hc_pf=gpf(Data=cdata$data,tree=cdata$phy,     
          frmla.phylo=competence~phylo,
          family=binomial,algorithm='phylo',nfactors=2,min.group.size=5)

## summarize
HolmProcedure(hc_pf)
hc_pf_results=pfsum(hc_pf)$results       

## save tree
cdata$data$infect=factor(cdata$data$pcr)
cdata$data$comp=factor(cdata$data$competence)
dtree=treeio::full_join(as.treedata(cdata$phy),cdata$data,by="label")

## fix palette
AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
afun=function(x){
  a=AlberColours[1:x]
  return(a)
}

## make low and high
pcols=afun(2)

## set x max
plus=1
pplus=plus+1

## fix taxa
pcr_pf_results$taxa[1]="Rodentia"
hc_pf_results$taxa[1]="italic(Felidae)"

## plot infection w/ ggtree
pcr_gg=ggtree(dtree,size=0.25)+
  geom_tippoint(aes(colour=infect),shape=15)+
  scale_colour_manual(values=c("grey80","black"))+
  guides(colour=F)   

## add clades to plot
for(i in 1:nrow(pcr_pf_results)){
  
  pcr_gg=pcr_gg+
    geom_hilight(node=pcr_pf_results$node[i],
                 alpha=0.25,
                 fill=ifelse(pcr_pf_results$clade>
                               pcr_pf_results$other,pcols[2],pcols[1])[i])+
    geom_cladelabel(node=pcr_pf_results$node[i],
                    label=pcr_pf_results$taxa[i],
                    offset=pplus,
                    hjust=0.75,
                    offset.text=pplus*2,
                    parse=T,
                    angle=90)
}
pcr_gg=pcr_gg


# plot competence
comp_gg=ggtree(dtree,size=0.25)+
  geom_tippoint(aes(colour=comp),shape=15)+
  scale_colour_manual(values=c("grey80","black"))+
  guides(colour=F)

## add clades to plot
for(i in 1:nrow(hc_pf_results)){
  
  comp_gg=comp_gg+
    geom_hilight(node=hc_pf_results$node[i],
                 alpha=0.25,
                 fill=ifelse(hc_pf_results$clade>
                               hc_pf_results$other,pcols[2],pcols[1])[i])+
    geom_cladelabel(node=hc_pf_results$node[i],
                    label=hc_pf_results$taxa[i],
                    offset=pplus,
                    hjust=0.75,
                    offset.text=pplus*2,
                    parse=T,
                    angle=90)
}
comp_gg=comp_gg

## print tree figures for infection and competence
library(ggpubr)
png("figs/Figure 1.png",width=6,height=6,units="in",res=300)
ggarrange(pcr_gg,comp_gg,ncol=2,widths=c(1.2,1),
          labels=c("(a) RT-PCR","(b) virus isolation"),
          label.x=c(-0.1,-0.2),
          font.label=list(face="plain",size=12))
dev.off()

## log1p pubmed cites
cdata$data$logcites=log1p(cdata$data$cites)

## add in pubmed for PCR
set.seed(1)
pcr_pf_pm=gpf(Data=cdata$data,tree=cdata$phy,
                 frmla.phylo=pcr~phylo,
                 weights=cdata$data$logcites,
                 family=binomial,algorithm='phylo',nfactors=10,min.group.size=5)

## summarize
HolmProcedure(pcr_pf_pm)
pcr_pf_pm_results=pfsum(pcr_pf_pm)$results       

## for competence
set.seed(1)
hc_pf_pm=gpf(Data=cdata$data,tree=cdata$phy,
                frmla.phylo=competence~phylo,
                weights=cdata$data$logcites,
                family=binomial,algorithm='phylo',nfactors=10,min.group.size=5)

## summarize
HolmProcedure(hc_pf_pm)
hc_pf_pm_results=pfsum(hc_pf_pm)$results       

## model citations themselves (not log1pm-transformed)
set.seed(1)
pm_pf=gpf(Data=cdata$data,tree=cdata$phy,
             frmla.phylo=cites~phylo,
             family=poisson,algorithm='phylo',nfactors=10,min.group.size=5)
HolmProcedure(pm_pf)
pm_pf_results=pfsum(pm_pf)$results


