setwd("~/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Fernandez Lab/Projects (Active)/OPV Host Prediction/GitHub/PoxHost")
data=read.csv('data/cleaned/pox cleaned response and traits.csv')
data=subset(data,select=-c(X.1))

data=subset(data,select=c(pcr,competence,studied,gen,fam,ord))

studied_data=data[data$studied==1,]
unique(studied_data$ord)
  
studied_data %>% count(gen, sort = TRUE)
studied_data %>% count(fam, sort = TRUE)

pos_data=data[data$studied==1&(data$pcr==1|data$competence==1),]
temp %>% count(gen, sort = TRUE)
unique(studied_data$fam[!studied_data$gen %in% temp$gen])


---
#gdata = ran up to line 244 of "01_cleaning.R" 
  gdata %>% count(fam,ord,sort=TRUE)
  gdata %>% count(ord, sort=TRUE)

sampled=gdata[gdata$studied==1,]
pos=gdata[gdata$studied==1&(data$pcr==1|data$competence==1),]

sampled %>% count(gen,fam,sort=TRUE)
sampled %>% count(fam,ord,sort=TRUE)
gen_sum_sampled <- subset(sampled_df, select=c(gen,pcr,competence))
fam_sum_sampled <- aggregate(cbind(pcr,competence) ~ fam, data=sampled_df, FUN = sum)

pos %>% count(gen,fam,sort=TRUE)
pos %>% count(fam,ord,sort=TRUE)
gen_sum_pos <- subset(pos_df, select=c(gen,pcr,competence))
fam_sum_pos <- aggregate(cbind(pcr,competence) ~ fam, data=pos_df, FUN = sum)




