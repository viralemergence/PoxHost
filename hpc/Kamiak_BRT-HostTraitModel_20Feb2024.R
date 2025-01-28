
# Libraries for BRT model
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
# library(phylofactor)
library(ggtree)
library(treeio)
library(caret) 
library(InformationValue)
library(mgcv) #for beta regression on performance metrics

# Clean environment
rm(list=ls()) 
graphics.off()

# # Set working directory
# setwd("~/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Fernandez Lab/Projects (Active)/OPV Host Prediction/GitHub/PoxHost/Tseng2022/Host Trait Model")


# Load data and clean environment
load("hosttrait_cleandata.RData")
data <- poxdata
rm(poxdata)

# Classify true negatives
data$type=ifelse(data$pcr==0 & data$competence==0,"true negative","other")

# Which species is competent but not PCR positive?
set=data
set$treename[set$pcr==0 & set$competence==1]

# Tabulate PCR/infection and isolation
set$inf=ifelse(set$pcr==0,"PCR negative","PCR positive")
set$iso=ifelse(set$competence==0,"no isolation","isolation")
table(set$inf,set$iso)

# Make binary variables for each taxonomic family; remove any duplicates
dums=dummy_cols(data["fam"])
dums=dums[!duplicated(dums$fam),]

# Ensure all family vars are factor
for(i in 1:ncol(dums)){
  dums[,i]=factor(dums[,i])
}

# Merge family taxa variables with dataset as predictors
data=merge(data,dums,by="fam",all.x=T)

# Drop unnecessary columns and clean environment
data$traitname=NULL
rm(dums,set,i)


# Create mode function where for each variable, we extract the frequency of the most frequently occurring value for that variable and divide it by the number of non-NA elements in that variable
mode.prop <- function(x) {                
  ux <- unique(x[is.na(x)==FALSE])        # creates array of unique values
  tab <- tabulate(match(na.omit(x), ux))  # creates array of the frequency (number of times) a unique value appears in a column 
  max(tab)/length(x[is.na(x)==FALSE])     # max-frequency / number of elements in each column that are not NA
}

# Assess variation across columns (2 indicates columns)
vars=data.frame(apply(data,2,function(x) mode.prop(x)),
                apply(data,2,function(x) length(unique(x))))    # number of unique elements in each column

# Get names
vars$variables=rownames(vars)
names(vars)=c("var","uniq","column")

# # Round values
# vars$var=round(vars$var,2)

# Label variables "cut" if homogeneous (100%)
vars$keep=ifelse(vars$var<1,"keep","cut")
vars$keep=ifelse(vars$column%in%c('fam','virus','gen','pcr','competence','fam'),'keep',vars$keep) # ensures we keep these columns
vars=vars[order(vars$keep),]

# Trim (creates array of column names to cut and removes from df)
keeps=vars[-which(vars$keep=="cut"),]$column

# Drop if no variation
data=data[keeps]
rm(keeps,vars)

# Assess missing values
mval=data.frame(apply(data,2,function(x) length(x[!is.na(x)])/nrow(data))) # proportion of values that are not NA

# Get names
mval$variables=rownames(mval)
names(mval)=c("comp","column")

# Exclude observations of non-predictor variables from mval dataframe
mval_hist <- mval[-(1:9),]
nrow(mval_hist)

# Plot frequency distribution of coverage among 185 traits
png("Kamiak_Results_HostTraitModel/histogram_trait_coverage.png", width=4,height=4,units="in",res=600)
ggplot(mval_hist[!mval_hist$column%in%c("gen","treename","pcr","competence","tip.label","fam"),],
       aes(comp))+
  geom_histogram(bins=50)+
  geom_vline(xintercept=0.60,linetype=2,linewidth=0.5)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  labs(y="Frequency",
       x="Trait coverage across mammal genera")+
  scale_x_continuous(labels = scales::percent)
dev.off()

# Label variables "cut" if >40% values are NA
# mval$keep=ifelse(mval$comp>=0.70,"keep","cut")
mval$keep=ifelse(mval$comp>=0.60,"keep","cut")
table(mval$keep)
mval=mval[order(mval$keep),]

# Trim (creates array of column names to cut and removes from df)
keeps=mval[-which(mval$keep=="cut"),]$column

# Drop if not well represented
data=data[keeps]

# Subset data to include only predictor and response variables
set <- subset(data,select=-c(gen,fam,ord,gtip,treename,type,studies,sampled))
## this dataframe will also be used in the next code block

# Create dataframe of trait coverage
table1=data.frame(apply(set,2,function(x) length(x[!is.na(x)])/nrow(set)))
## apply(X, MARGIN,...): "2" indicates applying the function over columns (as opposed to "1")

# Generate 'variable' column from rownames
table1$variables=rownames(table1)

# Rename columns and drop rownames
names(table1)=c("Coverage","Feature")
rownames(table1)=NULL

# Drop predictor variables from list
table1=table1[!table1$Feature%in%c("pcr","competence"),]
table1 <- subset(table1,select=c(Feature,Coverage))

# Save table of trait coverage
write.csv(table1, "Kamiak_Results_HostTraitModel/table_trait_coverage.csv")

# Check that binary variables are numeric and not factor (with the exception of fam_* variables)
str(set)

# Clean environment
rm(keeps,mval,mval_hist,table1)


# Hyperparameter tuning ifelse
#hok="ok"
hok="notok"
if(hok!="ok"){
  
  ## hyperparameter grid
  hgrid=expand.grid(n.trees=5000,   #creates df from all combinations of factor vars (1*3*3*1*10=90 obs & 5 vars)
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
  write.csv(search,"Kamiak_Results_HostTraitModel/par_tuning_data_summary.csv")
  
}else{
  
  ## load
  search=read.csv("Kamiak_Results_HostTraitModel/par_tuning_data_summary.csv")
  
}


# Load model tuning results
search=read.csv("Kamiak_Results_HostTraitModel/par_tuning_data_summary.csv")

# Convert parameters to factor and relabel values
search$shrinkage=factor(search$shrinkage)  #search=read.csv("Kamiak_Results_HostTraitModel/par_tuning_data_summary.csv")
lvl=rev(sort(unique(search$shrinkage)))    #sorts unique shrinkage par from large to small
search$shrinkage=factor(search$shrinkage,levels=lvl); rm(lvl)  #applies as factor
search$interaction.depth=factor(search$interaction.depth)
search$type=plyr::revalue(search$type,     #replace specified values w/ new values
                          c("PCR"="RT-PCR",
                            "competence"="virus isolation"))

# Fit beta regression model to the test AUCs of our PCR-based BRTs; ANOVA
mod=gam(testAUC~interaction.depth*shrinkage,
        data=search[search$type=="RT-PCR",],method="REML",family=betar)
anova(mod)

# Fit beta regression model to test AUCs of our competence-based BRTs; ANOVA
mod=gam(testAUC~interaction.depth*shrinkage,
        data=search[search$type=="virus isolation",],method="REML",family=betar)
anova(mod)

# Fit beta regression model to sensitivities of our PCR-based BRTs; ANOVA
mod=gam(sen~interaction.depth*shrinkage,
        data=search[search$type=="RT-PCR",],method="REML",family=betar)
anova(mod)

# Fit beta regression model to sensitivities of our competence-based BRTs; ANOVA
mod=gam(sen~interaction.depth*shrinkage,
        data=search[search$type=="virus isolation",],method="REML",family=betar)
anova(mod)

# Fit beta regression model to specificities of our PCR-based BRTs; ANOVA
mod=gam(spec~interaction.depth*shrinkage,
        data=search[search$type=="RT-PCR",],method="REML",family=betar)
anova(mod)

# Fit beta regression model to specificities of our competence-based BRTs; ANOVA
mod=gam(spec~interaction.depth*shrinkage,
        data=search[search$type=="virus isolation",],method="REML",family=betar)
anova(mod)

# To plot model tuning performance, transform dataframe from wide to long
search_long=gather(search,measure,value,testAUC:sen)

# Relabel values and convert to factor 
search_long$measure=plyr::revalue(search_long$measure,
                                  c("sen"="sensitivity",  
                                    "spec"="specificity",
                                    "testAUC"="test AUC"))
search_long$measure=factor(search_long$measure,
                           levels=c("test AUC","sensitivity","specificity"))

# Boxplot performance of model tuning w/ various parameter combinations
png("Kamiak_Results_HostTraitModel/boxplot_brt_tuning.png",width=5,height=8,units="in",res=600)
set.seed(1)
ggplot(search_long,aes(shrinkage,value,
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
  scale_color_viridis(discrete=TRUE,option="D")+
  scale_fill_viridis(discrete=TRUE,option="D")+
  guides(colour=guide_legend(title="interaction depth"),
         fill=guide_legend(title="interaction depth"))+
  labs(y=NULL,
       x="learning rate")+
  scale_y_continuous(n.breaks=4)
dev.off()

# To determine optimal parameters for model training, subset tuning results by number of trees
search_nt5000 <- search[search$n.trees==5000,]
search_nt15000 <- search[search$n.trees==15000,]
search_nt5000_sh0.01 <- search_nt5000[search_nt5000$shrinkage==0.010,]  #subset models with shrinkage==0.010
search_nt5000_sh0.001 <- search_nt5000[search_nt5000$shrinkage==0.001,]  #subset models with shrinkage==0.010

# Plot best.iter by evidence type (pcr vs. competence) to see max number of trees to include
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
search_nt5000_sh0.001 %>%
  ggplot( aes(x=best, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) 

# Clean
rm(search,search_long,hok,mod,search_nt5000,search_nt15000,search_nt5000_sh0.01, search_nt5000_sh0.001)


# BRT function to use different data partitions
brt_part=function(seed,response){
  
  ## Make new dataset
  ndata=set
  
  ## Correct response variable
  ndata$response=ndata[response][,1]
  
  ## Remove raw response variables
  ndata$pcr=NULL
  ndata$competence=NULL
  
  ## For BRT where cites is the response variable...
  if(response=="cites"){
    
    ## We add 1 to cites if cites equals 0
    ndata$cites=ifelse(ndata$cites==0,1,ndata$cites)
    
  }else{
    
    ndata=ndata
    
  }
  
  ## Use rsample package to split data
  set.seed(seed)
  split=initial_split(ndata,prop=0.7,strata="response")
  
  ## Create test and train datasets
  dataTrain=training(split)
  dataTest=testing(split)
  
  ## Get response variable for test dataset and train dataset
  yTrain=dataTrain$response
  yTest=dataTest$response
  
  ## Save distribution
  dist=ifelse(response=="cites","poisson","bernoulli")
  
  ## Save number of trees based on previous plots of optimal iterations
  nt=ifelse(response=="cites",10000,
            ifelse(response=="pcr",4500,5000)) #see plots of best.iter 
  
  ## Run BRTs using gbm package
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
  
  ## Get optimal number of iterations using gbm.perf & set plotting parameters for performance chart generated by gbm.perf
  par(mfrow=c(1,1),mar=c(4,4,1,1))                         
  best.iter=gbm.perf(gbmOut,method="cv")  #estimates optimal number of boosting iterations for a gbm object     
  
  ## Predict with test data, applying the optimal number of iterations as n.trees
  preds=predict(gbmOut,dataTest,n.trees=best.iter,type="response")
  
  ## Save known associations
  result=dataTest$response
  
  ## Get sensitivity and specificity
  sen=InformationValue::sensitivity(result,preds)
  spec=InformationValue::specificity(result,preds)
  
  ## Get AUC from model training
  auc_train=gbm.roc.area(yTrain,predict(gbmOut,dataTrain,n.trees=best.iter,type="response"))
  
  ## Get AUC from model test
  auc_test=gbm.roc.area(yTest,predict(gbmOut,dataTest,n.trees=best.iter,type="response"))
  
  ## Skip if poisson
  if(response=="cites"){
    
    perf=NA
    
  }else{
    
    ## Inner loop if yTest is all 0
    if(var(yTest)==0){
      
      perf=NA
    }else{
      
      ## To construct an ROC curve, create a prediction object using the predicted probabilities and the true class labels (known responses)
      pr=prediction(preds,dataTest$response)    
      
      ## Next, calculate the desired performance measures specified by the 'measures' argument
      perf=performance(pr,measure="tpr",x.measure="fpr")         #pr=prediction object; measure=performance measure for evaluation; x.measure=second perf measure (2-D)
      
      ## Create a dataframe of those performance values
      perf=data.frame(perf@x.values,perf@y.values)
      
      ## Rename columns
      names(perf)=c("fpr","tpr")
      
      ## Add seed
      perf$seed=seed
      
    }
  }
  
  ## Get relative importance
  bars=summary(gbmOut,n.trees=best.iter,plotit=F)
  bars$rel.inf=round(bars$rel.inf,2)
  
  ## Predict with cites
  preds=predict(gbmOut,data,n.trees=best.iter,type="response")
  pred_data=data[c("gtip",'treename',"fam","ord","pcr","competence")]
  pred_data$pred=preds
  pred_data$type=response
  
  ## Predict with mean cites
  pdata=data
  pdata$cites=mean(pdata$cites)
  pred_data$cpred=predict(gbmOut,pdata,n.trees=best.iter,type="response")
  
  ## Sort by decreasing predicted probability
  pred_data=pred_data[order(pred_data$pred,decreasing=T),]
  
  ## Print
  print(paste("BRT ",seed," done; test AUC = ",auc_test,sep=""))
  
  ## Save outputs
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


# Apply across 100 splits each
smax=100
pcr_brts=lapply(1:(smax),function(x) brt_part(seed=x,response="pcr"))
comp_brts=lapply(1:(smax),function(x) brt_part(seed=x,response="competence"))

# Run wos brts
pm_brts=lapply(1:(smax),function(x) brt_part(seed=x,response="cites"))

# Save results to output
saveRDS(pcr_brts, "Kamiak_Results_HostTraitModel/pcr_brts.rds")
saveRDS(comp_brts, "Kamiak_Results_HostTraitModel/comp_brts.rds")
saveRDS(pm_brts, "Kamiak_Results_HostTraitModel/pm_brts.rds")
