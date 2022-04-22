hgrid=expand.grid(n.trees=5000,                   
                  interaction.depth=c(2,3,4),
                  shrinkage=c(0.01,0.001,0.0005),
                  n.minobsinnode=4,
                  seed=seq(1,10,by=1))

## fix trees
hgrid$n.trees=ifelse(hgrid$shrinkage<0.001,hgrid$n.trees*3,hgrid$n.trees) 

## trees, depth, shrink, min, prop
hgrid$id=with(hgrid,paste(n.trees,interaction.depth,shrinkage,n.minobsinnode)) 

## sort by id then seed
hgrid=hgrid[order(hgrid$id,hgrid$seed),]

## now add rows
hgrid$row=1:nrow(hgrid)  #adds var 'row' based on row number in hgrid

## factor id
hgrid$id2=factor(as.numeric(factor(hgrid$id))) 

## function to assess each hyperpar combination
hfit=function(row,response){
  
  ## make new data
  ndata=set
  
  ## correct response
  ndata$response=ndata["pcr"][,1]  #creates var 'response'
  
  ## remove raw
  ndata$pcr=NULL
  ndata$competence=NULL
  
  ## use rsample to split
  set.seed(hgrid$seed[1])          #sets seed value of 1-10
  split=initial_split(ndata,prop=0.7,strata="response") #creates single binary split of data into training set and testing set, where 70% of data is retained for modeling/analysis and resampling is created within the 'response' var
  
  ## test and train
  dataTrain=training(split)
  dataTest=testing(split)
  
  ## yTest and yTrain
  yTrain=dataTrain$response  #create array of just response values from training and testing set
  yTest=dataTest$response
  
  ## BRT
  set.seed(1)
  gbmOut=gbm(response ~ . ,data=dataTrain,                     #y~x; dataframe containing vars in model
             n.trees=hgrid$n.trees[1],                       #total number of trees to fit (number of iterations; default is 100)
             distribution="bernoulli",
             shrinkage=hgrid$shrinkage[1],                   #equiv to learning rate or step-size reduction (smaller learning rate requires more trees, default is 0.1)
             interaction.depth=hgrid$interaction.depth[1],   #max depth of each tree (highest level of variable interactions allowed; default is 1)
             n.minobsinnode=hgrid$n.minobsinnode[1],         #min. number of obs in terminal nodes of trees
             cv.folds=5,class.stratify.cv=TRUE,                #number of cross-val folds to perform; for cv.folds>1, returns estimate of generalization error in 'cv.error'
             bag.fraction=0.5,train.fraction=1,                #fraction of training set obs randomly selected to propose next tree in expansion - this is why we set.seed()
             n.cores=1,                                        #number of CPU cores to use
             verbose=F)
  
  ## performance
  par(mfrow=c(1,1),mar=c(4,4,1,1))
  best.iter=gbm.perf(gbmOut,method="cv")                       #estimates optimal number of boosting iterations and plots 'training.error' performance measure; cv method extracts this optimal number using cross-validation
  
  ## predict with test data
  preds=predict(gbmOut,dataTest,n.trees=best.iter,type="response") #number of trees based on the optimal number of boosting iterations as set above
  
  ## known
  result=dataTest$response
  
  ## sensitiviy and specificity
  sen=InformationValue::sensitivity(result,preds, threshold=0.1)
  spec=InformationValue::specificity(result,preds, threshold=0.1)
  
  ## AUC on train
  auc_train=gbm.roc.area(yTrain,predict(gbmOut,dataTrain,n.trees=best.iter,type="response"))
  
  ## AUC on test
  auc_test=gbm.roc.area(yTest,predict(gbmOut,dataTest,n.trees=best.iter,type="response"))
  
  ## print
  print(paste("hpar row ",1," done; test AUC is ",auc_test,sep=""))
  
  ## save outputs
  #return(list(best=best.iter,
   return_test <- list(best=best.iter,            #saves optimal number of iterations, AUC on training set, AUC on testing set, specificity, sensitivity, and row number as a list
              trainAUC=auc_train,         
              testAUC=auc_test,           
              spec=spec,                  
              sen=sen,                    
              wrow=row)                   
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