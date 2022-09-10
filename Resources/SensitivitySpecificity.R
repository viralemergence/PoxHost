
# BRT Figures

### Load required packages and set system

```{r}

#(1) Libraries for BRT figures
library(ROCR)
library(tidyr)
library(ggplot2)
library(sciplot)
library(fastDummies)
library(caper)
library(ape)
library(phylofactor)
library(treeio)
library(ggtree)
library(plotrix)
library(rstatix)
library(ggrepel)
library(ggpubr)
library(plyr)

#(2) Clean environment
rm(list=ls()) 
graphics.off()

#(3) Set working directory
setwd("/Users/katietseng/Downloads/Results/Kamiak")

```

### Evaluate performance measures

```{r}

#### If needed, increase vector memory in R environment and reboot R before proceeding (https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached)

#(1) Load BRT results for PCR data
load("pcr_brts.RData")
#load("comp_brts.RData")

#(2) Explore PCR prediction results
pcr_brts[[1]]$predict
pcr_brts[[1]]$predict$pred

#These roc performance measures were calculated for each iteration using ROCR::prediction/performance functions
head(pcr_brts[[5]]$roc$fpr)
pcr_brts[[1]]$roc$tpr

#(3) Create list of ROC fpr and tpr values from PCR-model results
roc_out <- list()
roc_out$fpr <- list()
roc_out$tpr <- list()
for (i in 1:length(pcr_brts)) {
  roc_out$fpr[[i]] <- pcr_brts[[i]]$roc$fpr
  roc_out$tpr[[i]] <- pcr_brts[[i]]$roc$tpr
}
plot(roc_out$fpr, roc_out$tpr)
plot(do.call(cbind, roc_out), type = "l", 
     colorize=TRUE,lwd=2,main='ROC curves')
roc_fpr <- reshape2::melt(roc_out$fpr, value.name="fpr")
roc_tpr <- reshape2::melt(roc_out$tpr, value.name="tpr")
roc <- merge(roc_fpr,roc_tpr,by="L1")
roc <- cbind(roc_fpr$fpr, roc_tpr)
colnames(roc) <- c('fpr','tpr','seed')
ggplot(roc) + geom_line(aes(fpr, tpr, group = seed), 
                        alpha = 0.3)
plot(roc$fpr, roc$tpr, col=factor(roc$seed), type="l",colorize=TRUE,lwd=2,main='ROC curves')

#(3) Create list of PCR predictions and true values/labels from PCR-model results
pcr_out <- list()
pcr_out$predictions <- list()
pcr_out$labels <- list()
for (i in 1:length(pcr_brts)) {  
  pcr_out$predictions[[i]] <- pcr_brts[[i]]$predict$pred
  pcr_out$labels[[i]] <- pcr_brts[[i]]$predict$pcr
} 

#(4) Transform predictions & labels into standardized format using ROCR::prediction
pred = prediction(pcr_out$predictions,pcr_out$labels)
class(pred)
slotNames(pred)

#(5) ROC curves (x-axis: fpr, y-axis: tpr) using ROC::performance 
perf = performance(pred,'tpr','fpr')
slotNames(perf)
plot(perf,colorize=TRUE,lwd=2,main='ROC curves')
plot(perf,lwd=2,avg="vertical",spread.estimate="boxplot",main='ROC curves', add=TRUE)
# shows average vertical curve, boxplots and spread

#(6) Find optimal threshold (cutoff) for calculating sensitivity (tpr) and specificity (1-fpr)
cutoffs <- data.frame(cut=perf@alpha.values[[1]], fpr=perf@x.values[[1]], tpr=perf@y.values[[1]])
opt_threshold <- cutoffs[which.max((1-cutoffs$fpr) + cutoffs$tpr), "cut"]

#(7) Precision/recall curve (x-axis: recall, y-axis: precision)
perf <- performance(pred, "prec", "rec")
plot(perf)

#(8) Sensitivity/specificity curve (x-axis: specificity, y-axis: sensitivity)
perf <- performance(pred, "sens", "spec")
plot(perf)

#(9) AUC
perf <- performance(pred, "auc")
```

# ## get net AUC
# mean(c(sapply(pcr_brts,function(x) x$testAUC),sapply(comp_brts,function(x) x$testAUC)))
# se(c(sapply(pcr_brts,function(x) x$testAUC),sapply(comp_brts,function(x) x$testAUC)))
# 
# ## get net sensitivity
# mean(c(sapply(pcr_brts,function(x) x$sen),sapply(comp_brts,function(x) x$sen)))
# se(c(sapply(pcr_brts,function(x) x$sen),sapply(comp_brts,function(x) x$sen)))

## independent auc
mean(sapply(pcr_brts,function(x) x$testAUC))
se(sapply(pcr_brts,function(x) x$testAUC))
# mean(sapply(comp_brts,function(x) x$testAUC))
# se(sapply(comp_brts,function(x) x$testAUC))