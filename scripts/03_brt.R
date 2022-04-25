## PoxHost 03: mammal orthopoxvirus brt
## katie.tseng@wsu.edu
## updated 04/20/2022

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
  #packages <- c("gbm","fastDummies","rsample","ROCR","sciplot","ggplot2","pdp",
  #             "PresenceAbsence","tidyr","viridis","caper","phylofactor","ggtree","treeio","caret","InformationValue","mgcv")
  #install.packages(packages)

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
setwd("~/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Fernandez Lab/Projects (Active)/OPV Host Prediction/GitHub/PoxHost")
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

## make binary columns for family
dums=dummy_cols(data["fam"])

## unique
dums=dums[!duplicated(dums$fam),]

## ensure all family vars are factor
for(i in 1:ncol(dums)){
  
  ## column as factor
  dums[,i]=factor(dums[,i])
}

## merge
data=merge(data,dums,by="fam",all.x=T)
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
set$gen=NULL
set$treename=NULL
set$gtip=NULL
set$fam=NULL
set$ord=NULL
set$clade=NULL
set$type=NULL

## covariates
ncol(set)-2

## coverage table s1
ts1=data.frame(apply(set,2,function(x) length(x[!is.na(x)])/nrow(set)))

## get names
ts1$variables=rownames(ts1)
names(ts1)=c("comp","column")
rownames(ts1)=NULL

## trim disease data
ts1=ts1[!ts1$column%in%c("pcr","competence"),]
ts1$feature=ts1$column
ts1$column=NULL
ts1$coverage=ts1$comp
ts1$comp=NULL

## write file
write.csv(ts1, "figs/TableS1.csv")

## check that binary variables are numeric and not factor
str(set)

## hyperparameter tuning ifelse
#hok="ok"
hok="notok"
if(hok!="ok"){
  
  ## hyperparameter grid
  hgrid=expand.grid(n.trees=5000,                              #creates df from all combinations of factor vars (1*3*3*1*10=90 obs & 5 vars)
                    interaction.depth=c(2,3,4),
                    shrinkage=c(0.01,0.001,0.0005),
                    n.minobsinnode=4,
                    seed=seq(1,10,by=1))
  
  ## fix trees
  hgrid$n.trees=ifelse(hgrid$shrinkage<0.001,hgrid$n.trees*3,hgrid$n.trees)
  
  ## trees, depth, shrink, min, prop 
  hgrid$id=with(hgrid,paste(n.trees,interaction.depth,shrinkage,n.minobsinnode))   #creates var 'id' concatenating values from each of the specified columns in hgrid
  
  ## sort by id then seed
  hgrid=hgrid[order(hgrid$id,hgrid$seed),]
  
  ## now add rows
  hgrid$row=1:nrow(hgrid)                                        #adds var 'row' based on row number in hgrid
  
  ## factor id
  hgrid$id2=factor(as.numeric(factor(hgrid$id)))                 #creates 9-level factor var 'id2' 
  
  ## function to assess each hyperpar combination
  hfit=function(row,response){
    
    ## make new data
    ndata=set
    
    ## correct response
    ndata$response=ndata[response][,1]                           #creates var 'response'
    
    ## remove raw
    ndata$pcr=NULL
    ndata$competence=NULL
    
    ## use rsample to split
    set.seed(hgrid$seed[row])                                    #sets seed value of 1-10
    split=initial_split(ndata,prop=0.7,strata="response")        #creates single binary split of data into training set and testing set, where 70% of data is retained for modeling/analysis and resampling is created within the 'response' var
    
    ## test and train
    dataTrain=training(split)
    dataTest=testing(split)
    
    ## yTest and yTrain
    yTrain=dataTrain$response                                    #create array of just response values from training and testing set
    yTest=dataTest$response
    
    ## BRT
    set.seed(1)
    gbmOut=gbm(response ~ . ,data=dataTrain,                     #y~x; gbmOut contains list of 29 elements including train.error and valid.error referenced later in gbm.perf()
               n.trees=hgrid$n.trees[row],                       #total number of trees to fit (number of iterations; default is 100)
               distribution="bernoulli",
               shrinkage=hgrid$shrinkage[row],                   #equiv to learning rate or step-size reduction (smaller learning rate requires more trees, default is 0.1)
               interaction.depth=hgrid$interaction.depth[row],   #max depth of each tree (highest level of variable interactions allowed; default is 1)
               n.minobsinnode=hgrid$n.minobsinnode[row],         #min. number of obs in terminal nodes of trees
               cv.folds=5,class.stratify.cv=TRUE,                #number of cross-val folds to perform; for cv.folds>1, returns estimate of generalization error in 'cv.error'
               bag.fraction=0.5,train.fraction=1,                #fraction of training set obs randomly selected to propose next tree in expansion - this is why we set.seed()
               n.cores=4,                                        #number of CPU cores to use
               verbose=F)
    
    ## performance
    par(mfrow=c(1,1),mar=c(4,4,1,1))                             #sets graphical parameters such that subsequent figure are drawn in a nr-by-nc array by mfrows respectively and gives the number of lines of margin to be specified on the four sides of the plot c(bottom, L, top, R) -> see 'best.iter' plot below 
    best.iter=gbm.perf(gbmOut,method="cv")                       #estimates optimal number of boosting iterations and plots 'training.error' performance measure; cv method extracts this optimal number using cross-validation
    
    ## predict with test data
    preds=predict(gbmOut,dataTest,n.trees=best.iter,type="response")  #number of trees based on the optimal number of boosting iterations as set above (5,352)
    
    ## known
    result=dataTest$response
    
    ## sensitiviy and specificity                                #e.g., test run produced sensitivity of 0 b/c no predictedScores were > 0.5; and specificity of 1 b/c all predictedScores were <0.5
    sen=InformationValue::sensitivity(result,preds)              #calculates sensitivity (# of obs with event AND predicted to have event, divided by # of obs w/ event) for a given logit model where input is the actual binary flag (as numerica vector) for the response variable and the predicted probability scores for each observation; if predicted value is above the threshold (defaults to 0.5), it will be considered an event (1) or else a non-event (0)
    spec=InformationValue::specificity(result,preds)             #calculates specificity (# of obs w/o event AND predicted to not have event, divided by # of obs w/o event)  
    
    ## AUC on train
    auc_train=gbm.roc.area(yTrain,predict(gbmOut,dataTrain,n.trees=best.iter,type="response"))   #compute Information Retrieval measures for pairwise loss for a single group, where input is the observed value and the predicted value
    
    ## AUC on test
    auc_test=gbm.roc.area(yTest,predict(gbmOut,dataTest,n.trees=best.iter,type="response"))
    
    ## print
    print(paste("hpar row ",row," done; test AUC is ",auc_test,sep=""))  #prints "hpar row [x] done; test AUC is []"
    
    ## save outputs
    return(list(best=best.iter,                    #saves optimal number of iterations, AUC on training set, AUC on testing set, specificity, sensitivity, and row number as a list
                trainAUC=auc_train,
                testAUC=auc_test,
                spec=spec,
                sen=sen,
                wrow=row))
  }
  
  ## run the function
  hpars=lapply(1:nrow(hgrid),function(x) hfit(x,response="pcr"))
  
  ## get results
  hresults=data.frame(sapply(hpars,function(x) x$trainAUC),
                      sapply(hpars,function(x) x$testAUC),
                      sapply(hpars,function(x) x$spec),
                      sapply(hpars,function(x) x$sen),
                      sapply(hpars,function(x) x$wrow),
                      sapply(hpars,function(x) x$best))
  names(hresults)=c("trainAUC","testAUC",
                    "spec","sen","row","best")
  
  ## combine and save
  hsearch=merge(hresults,hgrid,by="row")
  
  ## save
  hsearch$type="PCR"
  
  ## rerun for competence
  hpars=lapply(1:nrow(hgrid),function(x) hfit(x,response="competence"))
  
  ## get results
  hresults=data.frame(sapply(hpars,function(x) x$trainAUC),
                      sapply(hpars,function(x) x$testAUC),
                      sapply(hpars,function(x) x$spec),
                      sapply(hpars,function(x) x$sen),
                      sapply(hpars,function(x) x$wrow),
                      sapply(hpars,function(x) x$best))
  names(hresults)=c("trainAUC","testAUC",
                    "spec","sen","row","best")
  
  ## combine and save
  csearch=merge(hresults,hgrid,by="row")
  
  ## assign data type
  csearch$type="competence"
  
  ## combine
  search=rbind.data.frame(csearch,hsearch)
  search$type=factor(search$type,levels=c("PCR","competence"))
  
  ## export
  setwd("~/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Fernandez Lab/Projects (Active)/OPV Host Prediction/GitHub/PoxHost")
  write.csv(search,"figs/par tuning data summary.csv")
  
}else{
  
  ## load
  setwd("~/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Fernandez Lab/Projects (Active)/OPV Host Prediction/GitHub/PoxHost")
  search=read.csv("figs/par tuning data summary.csv")
  
}

## factor parameters
search$shrinkage=factor(search$shrinkage)
lvl=rev(sort(unique(search$shrinkage)))  #sorts unique shrinkage parameters in order of largest to smallest
search$shrinkage=factor(search$shrinkage,levels=lvl); rm(lvl)  #applies as factor

## factor other
search$interaction.depth=factor(search$interaction.depth)

## fix type
search$type=plyr::revalue(search$type,    #replace specified values w/ new values
                          c("PCR"="RT-PCR",
                            "competence"="virus isolation"))

## PCR beta regression for AUC
mod=gam(testAUC~interaction.depth*shrinkage,
        data=search[search$type=="RT-PCR",],method="REML",family=betar)
anova(mod)

## virus isolation
mod=gam(testAUC~interaction.depth*shrinkage,
        data=search[search$type=="virus isolation",],method="REML",family=betar)
anova(mod)

## PCR beta regression for sensitivity
mod=gam(sen~interaction.depth*shrinkage,
        data=search[search$type=="RT-PCR",],method="REML",family=betar)
anova(mod)

## virus isolation
mod=gam(sen~interaction.depth*shrinkage,
        data=search[search$type=="virus isolation",],method="REML",family=betar)
anova(mod)


## PCR beta regression for specificity
mod=gam(spec~interaction.depth*shrinkage,
        data=search[search$type=="RT-PCR",],method="REML",family=betar)
anova(mod)

## virus isolation
mod=gam(spec~interaction.depth*shrinkage,
        data=search[search$type=="virus isolation",],method="REML",family=betar)
anova(mod)

## recast from wide to long
search2=gather(search,measure,value,testAUC:sen)

## revalue and factor (relabel values and change to factor)
search2$measure=plyr::revalue(search2$measure,
                              c("sen"="sensitivity",  
                                "spec"="specificity",
                                "testAUC"="test AUC"))
search2$measure=factor(search2$measure,
                       levels=c("test AUC","sensitivity","specificity"))


## visualize
setwd("~/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Fernandez Lab/Projects (Active)/OPV Host Prediction/GitHub/PoxHost")
png("figs/Figure S2.png",width=5,height=8,units="in",res=600)
set.seed(1)
ggplot(search2,aes(shrinkage,value,
                   colour=interaction.depth,fill=interaction.depth))+
  geom_boxplot(alpha=0.25)+
  geom_point(alpha=0.75,
             position = position_jitterdodge(dodge.width=0.75))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  facet_grid(measure~type,scales="free_y",switch="y")+
  theme(strip.placement="outside",
        strip.background=element_blank())+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        strip.text=element_text(size=12))+
  theme(legend.position="top")+
  scale_color_brewer(palette="Pastel2")+
  scale_fill_brewer(palette="Pastel2")+
  guides(colour=guide_legend(title="interaction depth"),
         fill=guide_legend(title="interaction depth"))+
  labs(y=NULL,
       x="learning rate")+
  scale_y_continuous(n.breaks=4)
dev.off()

## clean
rm(search,search2,hok,mod)

## brt function to use different data partitions
brt_part=function(seed,response){
  
  ## make new data
  ndata=set
  
  ## correct response
  ndata$response=ndata[response][,1]
  
  ## remove raw
  ndata$pcr=NULL
  ndata$competence=NULL
  
  ## fix cites if response
  if(response=="cites"){
    
    ## plus 1 for 0
    ndata$cites=ifelse(ndata$cites==0,1,ndata$cites)
    
  }else{
    
    ndata=ndata
    
  }
  
  ## use rsample to split
  set.seed(seed)
  split=initial_split(ndata,prop=0.7,strata="response")
  
  ## test and train
  dataTrain=training(split)
  dataTest=testing(split)
  
  ## yTest and yTrain
  yTrain=dataTrain$response
  yTest=dataTest$response
  
  ## dist
  dist=ifelse(response=="cites","poisson","bernoulli")
  
  ## n.trees
  nt=ifelse(response=="cites",10000,5000)
  
  ## BRT
  set.seed(1)
  gbmOut=gbm(response ~ . ,data=dataTrain,
             n.trees=nt,
             distribution=dist,
             shrinkage=0.001,
             interaction.depth=3,
             n.minobsinnode=4,
             cv.folds=5,class.stratify.cv=TRUE,
             bag.fraction=0.5,train.fraction=1,
             n.cores=1,
             verbose=F)
  
  ## performance
  par(mfrow=c(1,1),mar=c(4,4,1,1))                         
  best.iter=gbm.perf(gbmOut,method="cv")                   
  
  ## predict with test data
  preds=predict(gbmOut,dataTest,n.trees=best.iter,type="response")
  
  ## known
  result=dataTest$response
  
  ## sensitiviy and specificity
  sen=InformationValue::sensitivity(result,preds)
  spec=InformationValue::specificity(result,preds)
  
  ## AUC on train
  auc_train=gbm.roc.area(yTrain,predict(gbmOut,dataTrain,n.trees=best.iter,type="response"))
  
  ## AUC on test
  auc_test=gbm.roc.area(yTest,predict(gbmOut,dataTest,n.trees=best.iter,type="response"))
  
  ## skip if poisson
  if(response=="cites"){
    
    perf=NA
    
  }else{
    
    ## inner loop if yTest is all 0
    if(var(yTest)==0){
      
      perf=NA
    }else{
      
      ## ROC
      pr=prediction(preds,dataTest$response)                     
      perf=performance(pr,measure="tpr",x.measure="fpr")         #pr=prediction object; measure=performance measure for evaluation; x.measure=second perf measure (2-D)
      perf=data.frame(perf@x.values,perf@y.values)
      names(perf)=c("fpr","tpr")
      
      ## add seed
      perf$seed=seed
      
    }
  }
  
  ## relative importance
  bars=summary(gbmOut,n.trees=best.iter,plotit=F)
  bars$rel.inf=round(bars$rel.inf,2)
  
  ## predict with cites
  preds=predict(gbmOut,data,n.trees=best.iter,type="response")
  pred_data=data[c("gtip",'treename',"fam","ord","pcr","competence")]
  pred_data$pred=preds
  pred_data$type=response
  
  ## predict with mean cites
  pdata=data
  pdata$cites=mean(pdata$cites)
  pred_data$cpred=predict(gbmOut,pdata,n.trees=best.iter,type="response")
  
  ## sort
  pred_data=pred_data[order(pred_data$pred,decreasing=T),]
  
  ## print
  print(paste("BRT ",seed," done; test AUC = ",auc_test,sep=""))
  
  ## save outputs
  return(list(mod=gbmOut,
              best=best.iter,
              trainAUC=auc_train,
              testAUC=auc_test,
              spec=spec,
              sen=sen,
              roc=perf,
              rinf=bars,
              predict=pred_data,
              traindata=dataTrain,
              testdata=dataTest,
              seed=seed))
}

## apply across 10 splits each
smax=101
pcr_brts=lapply(1:smax,function(x) brt_part(seed=x,response="pcr"))
comp_brts=lapply(1:smax,function(x) brt_part(seed=x,response="competence"))

## write to files
setwd("~/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Fernandez Lab/Projects (Active)/OPV Host Prediction/GitHub/PoxHost")
saveRDS(pcr_brts,"data/cleaned/pcr brts.rds")
saveRDS(comp_brts,"data/cleaned/comp brts.rds")

## run wos brts
pm_brts=lapply(1:(smax-1),function(x) brt_part(seed=x,response="cites"))

## write
saveRDS(pm_brts,"data/cleaned/pm brts.rds")