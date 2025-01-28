# ##################################
# Title: "Kamiak_BRT_30Jan2024.R"
# Description: The following code was copied and pasted from LinkPredictionModel_Code.Rmd.
# install.packages(“”, dependencies = TRUE, repos = ‘http://cran.rstudio.com/')
# ##################################
# 
# 
# 3. BRT Model
# ============
# This chapter builds boosted regression tree models to predict host-virus links. The following section is coded to run on your local computer and will likely take ~72+ hours (with parallel processing on 5 cores), outputting an RData file size of 2GB. For tips on running the code on an available HPC node, see 'Tseng2022/HPC Example' on the PoxHost GitHub repository: https://github.com/viralemergence/PoxHost/tree/0a1effef83dbd5f6f3d88c6d0c15c563eb499452/Tseng2022/HPC%20Example_01Jun2023.
# 
# ### *Load required packages and set system*

# ```{r brt_load}

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
###to install ggtree, need to first install BiocManager:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("ggtree")
#library(phylofactor)
library(ggtree)
library(treeio)
library(caret) 
library(InformationValue)
library(mgcv) #for beta regression on performance metrics
# library(cowplot) #for combining plots of best.iter

# Clean environment
rm(list=ls()) 
graphics.off()

# # Set working directory
# setwd("~/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Fernandez Lab/Projects (Active)/OPV Host Prediction/GitHub/PoxHost/Tseng2022/Link Prediction Model")

# ```
# 
# ### *Reclassify (drop) feline poxvirus as cowpox virus ###
# 
# ```{r brt_}

# Load data and clean environment
load("linkpred_cleandata.RData")
data <- linkdata
rm(linkdata)

# Identify known interactions between feline poxvirus and mammal genera
data[data$virus=="feline poxvirus ita2_bc" & data$link==1, c('source','genus')]

# Identify whether known interaction between cowpox virus and Homo exists
data[data$virus=="cowpox virus" & data$link==1 & data$genus=="Homo", c('source','genus')]

# Drop all feline poxvirus pairings since cowpoxvirus-homo link is already in our dataset
nrow(data[data$virus=="feline poxvirus ita2_bc",])
data <- data[data$virus!="feline poxvirus ita2_bc",]

# ```
# 
# ### *Create variables of taxonomic family as predictors for the model*
# 
# ```{r brt_taxo}

# Calculate number of pseudoabsences
length(which(data$link==0))

# Ensure all accessory gene variables are numeric
PC_columns <- colnames(data[which(grepl("PC",names(data)))])
data[,c(PC_columns)] <- lapply(data[c(PC_columns)],as.numeric)
str(data)

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
rm(dums, i, PC_columns)

# ```
# 
# ### *Assess variation and availability of data*
# We explore the variation and availability of our data and drop variables if more than 40% of observations are missing/NA (vs. 30% for Host Trait Model; this is because of low coverage, ~60%, among PC/accessory gene variables). We also remove variables that are not needed for the BRT analysis and save the dataset as *data_LinkBRT.RData*.
# Figure: *histogram_trait_coverage.png*
#   Table: *table_trait_coverage.csv*
#   
# ```{r brt_var}

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
vars$keep=ifelse(vars$column%in%c('fam','source','sequence','link','virus','genus','ord','gtip','treename','cites','ed_equal'),'keep',vars$keep) # ensures we keep these columns
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

# Plot frequency distribution of coverage among traits
png("Kamiak_Results_Model3/histogram_trait_coverage.png", width=4,height=4,units="in",res=600)
ggplot(mval[!mval$column%in%c("gen","treename","pcr","competence","tip.label","fam"),],
       aes(comp))+
  geom_histogram(bins=50)+
  geom_vline(xintercept=0.70,linetype=2,size=0.5)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  labs(y="frequency",
       x="Trait coverage across mammal genera")+
  scale_x_continuous(labels = scales::percent)
dev.off()

# Label variables "cut" if >40% values are NA
mval$keep=ifelse(mval$comp>=0.60,"keep","cut")
table(mval$keep)
mval=mval[order(mval$keep),]

# Trim (creates array of column names to cut and removes from df)
keeps=mval[-which(mval$keep=="cut"),]$column

# Drop if not well represented
data=data[keeps]
rm(keeps,mval)

# Subset data to include only covariates
set <- subset(data,select=-c(fam, source, sequence, virus, genus, ord, gtip, treename))

#  Get trait coverage
trait_coverage=data.frame(apply(set,2,function(x) length(x[!is.na(x)])/nrow(set)))

# Rename and reorder columns
trait_coverage$variables=rownames(trait_coverage)
names(trait_coverage)=c("coverage","feature")
rownames(trait_coverage)=NULL
trait_coverage=trait_coverage[!trait_coverage$feature%in%c("pcr","competence"),]
trait_coverage <- subset(trait_coverage,select=c(feature,coverage))

# Save table
write.csv(trait_coverage, "Kamiak_Results_Model3/table_trait_coverage.csv")

# Check that binary variables are numeric and not factor (except for fam vars)
str(set)

# Remove vars not needed for BRT analysis
data$source=NULL
data$sequence=NULL
# save(data, file='Output/data_LinkBRT.RData')

# Clean environment
rm(keeps, mval, trait_coverage)
# 
# ```
# 
# ### *Alternative datasets Models 2 (excluding vaccinia virus sequences) and Model 3 (host traits only)*
# For Model 2, we prepare an alternative dataset that excludes vaccinia virus sequences (host associations with vaccinia virus) from BRT analysis. For Model 3, we prepare an alternative dataset that excludes viral traits from BRT analysis and is trained only on host traits (includes vaccinia virus).
# 
# ```{r brt_alt}
# 
# # Model 2: Subset data to exclude host-vaccinia virus links
# # data <- data[data$virus!="vaccinia virus",]
# # set <- subset(data,select=-c(fam, virus, genus, ord, gtip, treename))
# 
# Model 3: Subset data to exclude viral traits
# We prepare an alternative dataset that excludes viral traits from BRT analysis so that model is trained only on host traits.
#
# # Subset data to exclude viral traits
data <- subset(data,select=-c(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))
set <- subset(data,select=-c(fam, virus, genus, ord, gtip, treename))
# 
# ```
# 
# ### *Tuning to assess model performance for each combination of tuning parameters*
# We create a hyperparameter 'grid' that represents different combinations of parameter values to which we tune the model. We then save a table of the resulting model performance measures for 'Model 1', where we predict on all possible host-virus combinations. If model tuning on the alternative dataset (see above), {r brt_alt}, we recommend you modify the results table name to *par...Model2.csv*.
# Table: *par_tuning_data_summary_Model1.csv*
#   
#   ```{r brt_tuning}

# Hyperparameter tuning ifelse
#hok="ok"
hok="notok"
if(hok!="ok"){
  
  ## hyperparameter grid
  hgrid=expand.grid(n.trees=5000,                              #creates df from all combinations of factor vars (1*3*3*1*10=90 obs & 5 vars)
                    interaction.depth=c(2,3,4),
                    shrinkage=c(0.01,0.001,0.0005),
                    n.minobsinnode=4,
                    seed=seq(1,10,by=1))
  # hgrid=expand.grid(n.trees=500,                              #creates df from all combinations of factor vars (1*3*3*1*10=90 obs & 5 vars)
  #                   interaction.depth=c(2,3,4),
  #                   shrinkage=c(0.1,0.01,0.005),
  #                   n.minobsinnode=4,
  #                   seed=seq(1,10,by=1))
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
    # ndata$pcr=NULL
    # ndata$competence=NULL
    ndata$link=NULL
    
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
  
  ## run the function for link
  hpars=lapply(1:nrow(hgrid),function(x) hfit(x,response="link"))
  
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
  
  # ## save
  # hsearch$type="PCR"
  
  # ## rerun the function for competence
  # hpars=lapply(1:nrow(hgrid),function(x) hfit(x,response="competence"))
  # 
  # ## get results
  # hresults=data.frame(sapply(hpars,function(x) x$trainAUC),
  #                     sapply(hpars,function(x) x$testAUC),
  #                     sapply(hpars,function(x) x$spec),
  #                     sapply(hpars,function(x) x$sen),
  #                     sapply(hpars,function(x) x$wrow),
  #                     sapply(hpars,function(x) x$best))
  # names(hresults)=c("trainAUC","testAUC",
  #                   "spec","sen","row","best")
  # 
  # ## combine and save
  # csearch=merge(hresults,hgrid,by="row")
  # 
  # ## assign data type
  # csearch$type="competence"
  # 
  # ## combine
  # search=rbind.data.frame(csearch,hsearch)
  # search$type=factor(search$type,levels=c("PCR","competence"))
  
  search=hsearch
  
  ## export
  # write.csv(search,"Output/par_tuning_data_summary_Model1.csv")
  # write.csv(search,"Output/par_tuning_data_summary_Model2.csv")
  write.csv(search,"Kamiak_Results_Model3/par_tuning_data_summary_Model3.csv")
  
}else{
  
  ## 
  # search=read.csv("Output/par_tuning_data_summary_Model1.csv")
  # search=read.csv("Output/par_tuning_data_summary_Model2.csv")
  search=read.csv("Kamiak_Results_Model3/par_tuning_data_summary_Model3.csv")

}

# ```
# 
# ### *Assess model tuning results*
# We fit a beta regression model using mgcv::gam to explore the main interaction effects of tuning parameters, interaction depth and shrinkage, on performance metrics (AUC, sensitivity, and specificity). This model is appropriate for analyzing continuous response variables bounded b/w 0-1 (i.e,, beta distribution). Using ANOVA, we determine whether the coefficients b/w the two interactions depths, the two shrinkage rates, and the four possible interactions b/w them are significantly different. We then plot the performance of BRTs based on various parameter combinations. Finally, we assess the distribution of best.iter (optimal # of iterations) to determine the max # of trees to use for model training.
# Figure: *boxplot_brt_tuning.png*
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
# ```{r brt_tuning_results}

# search=read.csv("Output/par_tuning_data_summary_Model1.csv")
# search=read.csv("Output/par_tuning_data_summary_Model2.csv")
search=read.csv("Kamiak_Results_Model3/par_tuning_data_summary_Model3.csv")


# Convert parameters to factor and relabel values
search$shrinkage=factor(search$shrinkage)
lvl=rev(sort(unique(search$shrinkage)))  #sorts unique shrinkage par from large to small
search$shrinkage=factor(search$shrinkage,levels=lvl); rm(lvl)  #applies as factor
search$interaction.depth=factor(search$interaction.depth)

# Fit beta regression model to the test AUCs; ANOVA
mod=gam(testAUC~interaction.depth*shrinkage,  
        data=search,method="REML",family=betar)
anova(mod)

# Fit beta regression model to the sensitivities; ANOVA
mod=gam(sen~interaction.depth*shrinkage,
        data=search,method="REML",family=betar)
anova(mod)

# Fit beta regression model to the specificities; ANOVA
mod=gam(spec~interaction.depth*shrinkage,
        data=search,method="REML",family=betar)
anova(mod)

# To plot model tuning performance, transform dataframe from wide to long
search2=gather(search,measure,value,testAUC:sen)

# Relabel values and convert to factor 
search2$measure=plyr::revalue(search2$measure,
                              c("sen"="sensitivity",  
                                "spec"="specificity",
                                "testAUC"="test AUC"))
search2$measure=factor(search2$measure,
                       levels=c("test AUC","sensitivity","specificity"))

# Boxplot performance fo model tuning w/ various parameter combinations
png("Kamiak_Results_Model3/boxplot_brt_tuning.png",width=5,height=8,units="in",res=600)
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
  facet_grid(search2$measure,scales="free_y",switch="y")+
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

# To determine optimal parameters for model training, subset tuning results by number of trees
search_nt5000 <- search[search$n.trees==5000,]
search_nt15000 <- search[search$n.trees==15000,]
search_nt5000_sh0.01 <- search_nt5000[search_nt5000$shrinkage==0.010,]  #subset models with shrinkage==0.010

# Plot best.iter to see max number of trees to include
p1 <- search_nt5000 %>%
  # ggplot( aes(x=best, fill=type)) +
  ggplot( aes(x=best)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  ggtitle("(A) ntrees=5000 and shrinkage=0.001")
p2 <- search_nt15000 %>%
  # ggplot( aes(x=best, fill=type)) +
  ggplot( aes(x=best)) +  
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080"))  +
  ggtitle("(B) ntrees=15000 and shrinkage=0.001")
p3 <- search_nt5000_sh0.01 %>%
  # ggplot( aes(x=best, fill=type)) +
  ggplot( aes(x=best)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  ggtitle("(C) ntrees=5000 and shrinkage=0.01")

# Combine plots and save to output (using cowplot) - commenting out for HPC due to issues installing cowplot
# png("Output/histogram_brt_tuning_bestiter.png",width=12,height=4,units="in",res=600)
# plot_grid(p1,p2,p3) +
#   ggtitle("Histogram of best iteration")
# dev.off()

# Clean
rm(search,search2,hok,mod,search_nt5000,search_nt15000,search_nt5000_sh0.01,p1,p2,p3)

# ```
# 
# ### *BRT function for applying across multiple data partitions*
# We create our BRT function for model training, where we apply the function across multiple data partitions.
# 
# ```{r brt_partition}

# BRT function to use different data partitions
brt_part=function(seed,response){
  
  ## Make new dataset
  ndata=set
  
  ## Correct response variable
  ndata$response=ndata[response][,1]
  
  ## Remove raw response variable
  ndata$link=NULL
  
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
            ifelse(response=="link",4500,5000)) #see plots of best.iter 
  
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
  
  ## Get optimal number of iterations using gbm.perf & set plotting parameters for permance chart generated by gbm.perf
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
  #pred_data=data[c("gtip",'treename',"fam","ord","pcr","competence")]
  pred_data=data[c("virus","gtip",'treename',"fam","ord","link")]
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

# ```
# 
# ### *Apply BRT function across 100 partitions to generate ensemble*
# 
# ```{r brt_ensemble}

# Apply across 100 splits each
# smax=101
smax=100
brts=lapply(1:smax,function(x) brt_part(seed=x,response="link"))

# Run wos brts
pm_brts=lapply(1:(smax),function(x) brt_part(seed=x,response="cites"))

# Save results to working directory
# save(brts,pm_brts,file="Output/brts_Model1.RData")
# save(brts,pm_brts,file="Output/brts_Model2.RData")
save(brts,pm_brts,file="Kamiak_Results_Model3/brts_Model3.RData")

# ```



