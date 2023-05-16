# pcrpred_pf_results=pfsum(pcrpred_pf)$results
# comppred_pf_results=pfsum(comppred_pf)$results


#(3) Holm rejection procedure (counteract the problem of multiple comparisons and controls FWER)
HolmProcedure <- function(pf,FWER=0.05){
  
  ## get split variable
  cs=names(coef(pcrpred_pf$models[[1]]))[-1] #returns names of the coefficients in the 1st model minus "Intercept" (the first coefficient name)
  split=ifelse(length(cs)>1,cs[3],cs[1]) #if list has more than 1 name, returns 3rd coefficient name, or else returns first name (e.g., "phy)
  
  ## obtain p values
  if (pcrpred_pf$models[[1]]$family$family%in%c('gaussian',"Gamma","quasipoisson")){  #if c('[etc.]') is in family$family...
    pvals <- sapply(pcrpred_pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|t|)'])  
  } else {
    pvals <- sapply(pcrpred_pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|z|)'])
  }
  D <- length(pcrpred_pf$tree$tip.label) #946L
  
  ## this is the line for Holm's sequentially rejective cutoff
  keepers <- pvals<=(0.05/(2*D-3 - 2*(0:(pcrpred_pf$nfactors-1))))
  
  
  if (!all(keepers)){
    nfactors <- min(which(!keepers))-1
  } else {
    nfactors <- pcrpred_pf$nfactors
  }
  return(nfactors)
}

###nfactors=3

#(4) Summarize pf object 
pfsum=function(pf){
  
  ## get formula
  chars=as.character(pcrpred_pf$frmla.phylo)[-1] #saves all elements except the first
  
  ## response
  resp=chars[1] #saves first item in the array
  
  ## fix
  resp=ifelse(resp=='cbind(pos, neg)','prevalence',resp) #return either 'prevalence' or resp based on conditional
  
  ## holm
  # hp=HolmProcedure(pcrpred_pf) #see line 80
  hp=3
  
  ## save model
  model=chars[2]
  
  ## set key
  setkey(pcrpred_pf$Data,'Species')
  
  ## make data
  dat=data.frame(pcrpred_pf$Data)
  
  ## make clade columns in data
  for(i in 1:hp){
    
    dat[,paste0(resp,'_pf',i)]=ifelse(dat$Species%in%cladeget(pcrpred_pf,i),'factor','other')
    
  }
  
  ## make data frame to store taxa name, response, mean, and other
  results=data.frame(matrix(ncol=6, nrow = hp))
  colnames(results)=c('factor','taxa','tips','node',"clade",'other')
  
  ## set taxonomy
  taxonomy=dat[c('Species','taxonomy')]
  taxonomy$taxonomy=as.character(taxonomy$taxonomy)
  
  ## loop
  for(i in 1:hp){
    
    ## get taxa
    tx=pf.taxa(pcrpred_pf,taxonomy,factor=i)$group1
    
    ## get tail
    tx=sapply(strsplit(tx,'; '),function(x) tail(x,1))
    
    ## combine
    tx=paste(tx,collapse=', ')
    
    # save
    results[i,'factor']=i
    results[i,'taxa']=tx
    
    ## get node
    tips=cladeget(pcrpred_pf,i)
    node=ggtree::MRCA(pcrpred_pf$tree,tips)
    results[i,'tips']=length(tips)
    results[i,'node']=ifelse(is.null(node) & length(tips)==1,'species',
                             ifelse(is.null(node) & length(tips)!=1,NA,node))
    
    ## get means
    ms=(tapply(dat[,resp],dat[,paste0(resp,'_pf',i)],mean))
    
    ## add in
    results[i,'clade']=ms['factor']
    results[i,'other']=ms['other']
    
  }
  
  ## return
  # return(list(set=dat,results=results))
  temp <- list(set=dat,results=results)
  
}



