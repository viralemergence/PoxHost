##################################
#Title: Kamiak_BRT_07Sep2022.R
#Description: The following code was copied and pasted from HostPrediction_Code.Rmd
##with the exception of lines 601-603 in which model output is saved into 
##separate files due to file size constraints.
#install.packages(“”, dependencies = TRUE, repos = ‘http://cran.rstudio.com/')
###################################

# # 3. Boosted Regression Trees
# 
# ### Load required packages and set system
# 
# ```{r}

#(1) Libraries for BRT model
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

#(2) Clean environment
rm(list=ls()) 
graphics.off()

# #(3) Set working directory
# setwd("~/Tseng2022")

# ```
# 
# ### Create taxonomic variables as predictors for the model
# 
# ```{r}

#(1) Load data and clean environment
load("Output/HostData_clean.RData")
data <- poxdata
rm(poxdata)

#(2) Classify true negatives
data$type=ifelse(data$pcr==0 & data$competence==0,"true negative","other")

#(3) Which species is competent but no PCR record?
set=data
set$treename[set$pcr==0 & set$competence==1]

#(4) Tabulate PCR/infection and isolation
set$inf=ifelse(set$pcr==0,"PCR negative","PCR positive")
set$iso=ifelse(set$competence==0,"no isolation","isolation")
table(set$inf,set$iso)

#(5) Make binary variables for each taxonomic family; remove any duplicates
dums=dummy_cols(data["fam"])
dums=dums[!duplicated(dums$fam),]

#(6) Ensure all family vars are factor
for(i in 1:ncol(dums)){
  dums[,i]=factor(dums[,i])
}

#(7) Merge family taxa variables with dataset as predictors
data=merge(data,dums,by="fam",all.x=T)

#(8) Drop unnecessary columns and clean environment
data$traitname=NULL
rm(dums,set,ag_columns)

# ```
# 
# ### Assess variation and availability of data
# 
# ```{r}

#(1) Mode function
mode.prop <- function(x) {                
  ux <- unique(x[is.na(x)==FALSE])        # creates array of unique values
  tab <- tabulate(match(na.omit(x), ux))  # creates array of the frequency (number of times) a unique value appears in a column 
  max(tab)/length(x[is.na(x)==FALSE])     # max-frequency / number of elements in each column that are not NA
}

#(2) Assess variation across columns (2 indicates columns)
vars=data.frame(apply(data,2,function(x) mode.prop(x)),
                apply(data,2,function(x) length(unique(x))))    # number of unique elements in each column

#(3) Get names
vars$variables=rownames(vars)
names(vars)=c("var","uniq","column")

# ## round values
# vars$var=round(vars$var,2)

#(4)Label variables "cut" if homogeneous (100%)
vars$keep=ifelse(vars$var<1,"keep","cut")
vars$keep=ifelse(vars$column%in%c('fam','virus','gen','pcr','competence','fam'),'keep',vars$keep) # ensures we keep these columns
vars=vars[order(vars$keep),]

#(5) Trim (creates array of column names to cut and removes from df)
keeps=vars[-which(vars$keep=="cut"),]$column

#(6) Drop if no variation
data=data[keeps]
rm(keeps,vars)

#(7) Assess missing values
mval=data.frame(apply(data,2,function(x) length(x[!is.na(x)])/nrow(data))) # proportion of values that are not NA

#(8) Get names
mval$variables=rownames(mval)
names(mval)=c("comp","column")
# 
# #(9) visualize distribution of NA
# png("Output/Figure S1.png", width=4,height=4,units="in",res=600)
# ggplot(mval[!mval$column%in%c("gen","treename","pcr","competence","tip.label","fam"),],
#        aes(comp))+
#   geom_histogram(bins=50)+
#   geom_vline(xintercept=0.70,linetype=2,size=0.5)+
#   theme_bw()+
#   theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
#   theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
#   theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
#   labs(y="frequency",
#        x="trait coverage across mammal species (genus)")+
#   scale_x_continuous(labels = scales::percent)
# dev.off()

#(10) Label variables "cut" if >30% values are NA
mval$keep=ifelse(mval$comp>=0.70,"keep","cut")
table(mval$keep)
mval=mval[order(mval$keep),]

#(11) Trim (creates array of column names to cut and removes from df)
keeps=mval[-which(mval$keep=="cut"),]$column

#(12) Drop if not well represented
data=data[keeps]
rm(keeps,mval)

#(14) Save list of covariates and their coverage as table S1
set <- subset(data,select=-c(gen,fam,ord,gtip,treename,type,studies,sampled))
ts1=data.frame(apply(set,2,function(x) length(x[!is.na(x)])/nrow(set)))

#(15) Rename and reorder columns
ts1$variables=rownames(ts1)
names(ts1)=c("coverage","feature")
rownames(ts1)=NULL
ts1=ts1[!ts1$feature%in%c("pcr","competence"),]
ts1 <- subset(ts1,select=c(feature,coverage))

#(16) Save Table S1 to results
write.csv(ts1, "Output/TableS1.csv")

#(17) Check that binary variables are numeric and not factor (with the exception of fam_* variables)
str(set)

# ```
# 
# ### Model tuning to asses model performance for each combination of tuning parameters
# 
# ```{r}

#(1) Hyperparameter tuning ifelse
#hok="ok"
hok="notok"
if(hok!="ok"){
  
  ## hyperparameter grid
  hgrid=expand.grid(n.trees=5000,                              #creates df from all combinations of factor vars (1*3*3*1*10=90 obs & 5 vars)
                    interaction.depth=c(2,3,4),
                    shrinkage=c(0.01,0.001,0.0005),
                    n.minobsinnode=4,
                    seed=seq(1,10,by=1))
  
  # fix trees
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
               cv.folds=5,class.stratify.cv=TRUE,                #no. of cross-val folds to perform; for cv.folds>1, returns estimate of generalization error in 'cv.error'
               bag.fraction=0.5,train.fraction=1,                #fraction of training set obs randomly selected to propose next tree in expansion - this is why we set.seed()
               n.cores=5,                                        #no. of CPU cores to use
               verbose=F)
    # par.details=(gbmParallel(num_threads=5)),
    
    ## performance
    par(mfrow=c(1,1),mar=c(4,4,1,1))                             #sets graphical parameters such that subsequent figure are drawn in a nr-by-nc array by mfrows respectively and gives the number of lines of margin to be specified on the four sides of the plot c(bottom, L, top, R) -> see 'best.iter' plot below 
    best.iter=gbm.perf(gbmOut,method="cv")                       #estimates optimal number of boosting iterations and plots 'training.error' performance measure; cv method extracts this optimal number using cross-validation
    
    ## predict with test data
    preds=predict(gbmOut,dataTest,n.trees=best.iter,type="response")  #number of trees based on the optimal number of boosting iterations as set above (5,352)
    
    ## known
    result=dataTest$response
    
    # ##estimate threshold value for classification of predicted probability
    # #library(pROC)
    # analysis <- roc(result,preds)  #roc([actual values],[predicted values])
    # e <- cbind(analysis$thresholds,analysis$sensitivities+analysis$specificities) #pulls each array and binds them into dataframe: 1st column are thresholds, 2nd column are sensitivities + specificities
    # 
    # ##optimum threshold value
    # opt_t <- subset(e,e[,2]==max(e[,2]))[,1] #subsets dataframe and returns the max (sens+spec) value of 2nd column of e 
    # #threshold<-opt_t #set as threshold value
    # #threshold = 0.2
    
    ## sensitivity and specificity                              #e.g., test run produced sensitivity of 0 b/c no predictedScores were > 0.5; and specificity of 1 b/c all predictedScores were <0.5
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
  
  ## run the function for PCR
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
  
  ## rerun the function for competence
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
  write.csv(search,"Output/par_tuning_data_summary.csv")
  
}else{
  
  ## load
  search=read.csv("Output/par_tuning_data_summary.csv")
  
}

# ```
# 
# ### Model tuning results: Figure S2
# ```{r}

#(1) Convert parameters to factor and relabel values
search$shrinkage=factor(search$shrinkage)
lvl=rev(sort(unique(search$shrinkage)))  #sorts unique shrinkage par from large to small
search$shrinkage=factor(search$shrinkage,levels=lvl); rm(lvl)  #applies as factor
search$interaction.depth=factor(search$interaction.depth)
search$type=plyr::revalue(search$type,    #replace specified values w/ new values
                          c("PCR"="RT-PCR",
                            "competence"="virus isolation"))

#(2) PCR beta regression for AUC
mod=gam(testAUC~interaction.depth*shrinkage,   #gen additive models (gam) w/ integrated smoothness estimation
        data=search[search$type=="RT-PCR",],method="REML",family=betar)
anova(mod)

#(3) Competence beta regression for AUC
mod=gam(testAUC~interaction.depth*shrinkage,
        data=search[search$type=="virus isolation",],method="REML",family=betar)
anova(mod)

#(4) PCR beta regression for sensitivity
mod=gam(sen~interaction.depth*shrinkage,
        data=search[search$type=="RT-PCR",],method="REML",family=betar)
anova(mod)

#(5) Competence beta regression for sensitivity
mod=gam(sen~interaction.depth*shrinkage,
        data=search[search$type=="virus isolation",],method="REML",family=betar)
anova(mod)


#(6) PCR beta regression for specificity
mod=gam(spec~interaction.depth*shrinkage,
        data=search[search$type=="RT-PCR",],method="REML",family=betar)
anova(mod)

#(7) Competence beta regression for specificity
mod=gam(spec~interaction.depth*shrinkage,
        data=search[search$type=="virus isolation",],method="REML",family=betar)
anova(mod)

#(8) Recast from wide to long
search2=gather(search,measure,value,testAUC:sen)

#(9) Relabel values and convert to factor 
search2$measure=plyr::revalue(search2$measure,
                              c("sen"="sensitivity",  
                                "spec"="specificity",
                                "testAUC"="test AUC"))
search2$measure=factor(search2$measure,
                       levels=c("test AUC","sensitivity","specificity"))

#(10) Visualize - Figure S2
png("Output/FigureS2.png",width=5,height=8,units="in",res=600)
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

#(11) To determine optimal parameters for model training, subset tuning results by number of trees
search_nt5000 <- search[search$n.trees==5000,]
search_nt15000 <- search[search$n.trees==15000,]
search_nt5000_sh0.01 <- search_nt5000[search_nt5000$shrinkage==0.010,]  #subset models with shrinakge==0.010

#(12) Plot best.iter by type (pcr/competence) to see max number of trees to include
search_nt5000 %>%
  ggplot( aes(x=best, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) 
search_nt15000 %>%
  ggplot( aes(x=best, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) 
search_nt5000_sh0.01 %>%
  ggplot( aes(x=best, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) 

#(13) Clean
rm(search,search2,hok,mod,search_nt5000,search_nt15000,search_nt5000_sh0.01)

# ```
# 
# ### BRT function for applying across multiple data partitions
# 
# ```{r}

#(1) BRT function to use different data partitions
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
  nt=ifelse(response=="cites",10000,
            ifelse(response=="pcr",4500,5000)) #see plots of best.iter 
  
  ## BRT
  set.seed(1)
  gbmOut=gbm(response ~ . ,data=dataTrain,
             n.trees=nt,
             distribution=dist,
             shrinkage=0.01, #see plots of best.iter 
             interaction.depth=3,
             n.minobsinnode=4,
             cv.folds=5,class.stratify.cv=TRUE,
             bag.fraction=0.5,train.fraction=1,
             n.cores=5,
             verbose=F)
  # par.details=(gbmParallel(num_threads=5)),
  
  ## performance
  par(mfrow=c(1,1),mar=c(4,4,1,1))                         
  best.iter=gbm.perf(gbmOut,method="cv")  #estimates optimal number of boosting iterations for a gbm object                 
  
  ## predict with test data
  preds=predict(gbmOut,dataTest,n.trees=best.iter,type="response")
  
  ## known
  result=dataTest$response
  
  ## sensitivity and specificity
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
  #pred_data=data[c("gtip",'treename',"fam","ord","pcr","competence")]
  pred_data=data[c("virus","gtip",'treename',"fam","ord","pcr","competence")]
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

# ```
# 
# ### Apply BRT function across 100 partitions to generate ensemble
# 
# ```{r}

#(1) Apply across 100 splits each
# smax=101
smax=100
pcr_brts=lapply(1:smax,function(x) brt_part(seed=x,response="pcr"))
comp_brts=lapply(1:smax,function(x) brt_part(seed=x,response="competence"))

#(2) Run wos brts
pm_brts=lapply(1:(smax-1),function(x) brt_part(seed=x,response="cites"))

#(3) Save results to wd
# save(pcr_brts,comp_brts,pm_brts,file="Output/HostData_results.RData")
save(pcr_brts,file="pcr_brts.RData")
save(comp_brts,file="comp_brts.RData")
save(pm_brts,file="pm_brts.RData")

# ```
