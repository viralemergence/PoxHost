
#### HolmProcedure ####

## First, let's try replacing pf argument with pcr_pf
holmTest = HolmProcedure(pcr_pf, FWER = 0.05)               # 0 ???

## Now, let's try line by line
FWER=0.05
cs_test=(names(coef(pcr_pf$models[[1]]))[-1])               # "phyloS"
split_test=ifelse(length(cs_test)>1,cs_test[3],cs_test[1])  # "phyloS"
unique(pcr_pf$models[[1]]$family$family)                    # "binomial"

if (pcr_pf$models[[1]]$family$family%in%c('gaussian',"Gamma","quasipoisson")){                  #if family$family of the 1st list element of pf$models is in the columns 'gaussian', etc...
  pvals_test <- sapply(pcr_pf$models,FUN=function(fit) summary(fit)$coefficients[split_test,'Pr(>|t|)'])  #then, to each element of pf$models, we apply the summary function with the argument 'fit' and assign the output to 'pvals'; specifically, we use 'summary(fit)' to call the output of 'pf$models', extracting the 'coefficients' section, whereby we index the column named 'Pr(>|t\)' and split the data in that column; see sample output of linear model of R for reference (https://feliperego.github.io/blog/2015/10/23/Interpreting-Model-Output-In-R) 
} else {
  pvals_test <- sapply(pcr_pf$models,FUN=function(fit) summary(fit)$coefficients[split_test,'Pr(>|z|)'])  #or else, extract p-val based on z statistic
}

D_test <- length(pcr_pf$tree$tip.label)                     # 87L

keepers_test <- pvals_test<=(FWER/(2*D_test-3 - 2*(0:(pcr_pf$nfactors-1))))    # FALSE FALSE FALSE                           
if (!all(keepers_test)){                                    
  nfactors_test <- min(which(!keepers_test))-1              # 1-1
} else {
  nfactors_test <- pcr_pf$nfactors                          
}
nfactors_test                                               # 0 
# end function


#### cladeget ####
factor==1
spp_test=pcr_pf$tree$tip.label[pcr_pf$groups[[1]][[1]]]   


#### pfsum ####
chars_test=as.character(pcr_pf$frmla.phylo)[-1]             # "pcr" "phylo"
resp_test=chars_test[1]                                     # "pcr"
hp_test=HolmProcedure(pcr_pf)                               # 0
model_test=chars_test[2]                                    # "phylo"
setkey(pcr_pf$Data,'Species')                               
dat_test=data.frame(pcr_pf$Data)                            

# for(i in 1:hp_test){  # returns error b/c hp_test=0; let's assume i=1
    # dat_test[,paste0(resp_test,'_pf',i)]=ifelse(dat_test$Species%in%cladeget(pcr_pf,i),'factor','other')
    # }
    paste0(resp_test,'_pf',1)                                         # "pcr_pf1"
    ifelse(dat_test$Species%in%cladeget(pcr_pf,1),'factor','other')   # lists 'other' or 'factor'
    dat_test[,paste0(resp_test,'_pf',1)]=ifelse(dat_test$Species%in%cladeget(pcr_pf,1),'factor','other')  # creates variable/list in dat_test named "pcr_pf1" set equal to the aforementioned list

results_test=data.frame(matrix(ncol=6,nrow=hp_test))
# returns 0x6 matrix b/c hp_test=0, so let's assume hp_test=1
results_test=data.frame(matrix(ncol=6,nrow=1))
colnames(results_test)=c('factor','taxa','tips','node',"clade",'other')
taxonomy=dat_test[c('Species','taxonomy')]
taxonomy$taxonomy=as.character(taxonomy$taxonomy)


# 'for(i in 1:hp_test)' returns error b/c hp_test=0; let's assume i=1
    tx_test=pf.taxa(pcr_pf,taxonomy,factor=1)$group1                            # gets taxonomic order
    tx_test=sapply(strsplit(tx_test,'; '),function(x) tail(x,1))                # gets taxonomic family as list
    tx_test=paste(tx_test,collapse=', ')                                        # collapses taxonomic family into single string
  
    results_test[1,'factor']=1
    results_test[1,'taxa']=tx_test                                              # fill results table
  
    tips_test=cladeget(pcr_pf,1)
    node_test=ggtree::MRCA(pcr_pf$tree,tips_test)                               # 103L
    results_test[1,'tips']=length(tips_test)                                    # 9
    results_test[1,'node']=ifelse(is.null(node_test) & length(tips_test)==1,'species',
                             ifelse(is.null(node_test) & length(tips_test)!=1,NA,node_test)) 
    
    ms_test=(tapply(dat_test[,resp_test],dat_test[,paste0(resp_test,'_pf',1)],FUN=mean))   # 1 0.628
    
    results_test[1,'clade']=ms_test['factor']                                     # returns 1          
    results_test[1,'other']=ms_test['other']                                      # returns 0.628                  
  
    list(set_test=dat_test,results_test=results_test)
    
   