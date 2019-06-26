###  Perform anaova analysis and exclude everything above and equal p.values =>0.9

# Side step 

# We have four Groups: 

# 1: control-high score

# 2: case-high score

# 3: case-low score 

# 4: control-low score

#  and I will  join the 1-2 and 3-4 groups together 

# In thid first step I am interested control versus cases

rm(list = ls())
output.path <- "~/Documents/PhD_2016-2019/projects//Baltica_Folder/Baltica/Baltica_project/results/"

# load the data
load("~/Documents/PhD_2016-2019/projects/Baltica_Folder/Baltica/Baltica_project/results/data/BALTICA_rawdata_18052019.rda")
load("~/Documents/PhD_2016-2019/projects/Baltica_Folder/Baltica/Baltica_project/results/data/CLEAN_DRIFT_allmodes.Rdata")
load("~/Documents/PhD_2016-2019/projects/Baltica_Folder/Baltica/Baltica_project/results/data/covariates_data.rda")
load("~/Documents/PhD_2016-2019/projects/Baltica_Folder/Baltica/Baltica_project/results/data/full_compound_data.Rdata")
load("~/Documents/PhD_2016-2019/projects/Baltica_Folder/Baltica/Baltica_project/results/model/MetaboliteRanks_MUVRmodel16042019.rda")

# apply glog 
library(FitAR)
log_baltica <-glog(baltica, a = 1, InverseQ = FALSE)



# 
# # apply normalization
library(KODAMA)

# Normalization 
# I will perform PQN normalization in order to do that I will need a reference sample and for that 
# I will use a QC sample the log.average.QC which the mean of all QCs in the Baltica data
# Extract only qcs 
only.qcs <- clean.data.all.modes.with.qcs[which(clean.data.all.modes.with.qcs$GROUP=='QC'), ]
# calculate the average of the QCs in the data
log.only.qcs <- glog(only.qcs[,-c(1,2)], a = 1, InverseQ = FALSE)


log.average.QC <- as.data.frame(colMeans(log.only.qcs))



data_norm_temp <- KODAMA::normalization(log_baltica, method = "pqn",ref=log.average.QC)
data_norm <- data_norm_temp$newXtrain



new_data <- cbind.data.frame(GROUP,data_norm)


# combine all the controls together using the  subset function 
Western <-  subset(new_data, GROUP=="1" | GROUP=="2" )


# combine all the cases together using the  subset function 

Healthy <-  subset(new_data, GROUP=="3" | GROUP=="4" )



# Create new labels 0 for controls and 1 for cases

Healthy_label <-t(rep(0,dim(Healthy)[1]))
Western_label <- t(rep(1,dim(Western)[1]))
New_Group <-c(Healthy_label,Western_label)

# combine Case and control data with the new labels

new_data2 <- rbind.data.frame(Healthy,Western)
new_data2 <-cbind.data.frame(New_Group,new_data2)

# extract the data by removing the first two columns which contain the GROUP info and the RUN_ORDER




new.covariates.data <- cbind.data.frame(GROUP,covariates.data)

# combine all the controls together using the  subset function 
West <-  subset(new.covariates.data, GROUP=="1" | GROUP=="2" )


# combine all the cases together using the  subset function 

Health <-  subset(new.covariates.data, GROUP=="3" | GROUP=="4" )


final_covariates <- rbind.data.frame(Health,West)

# Western Versus Healthy 

# 1st model without adjustements

 # Control Versus Healthy 
  
  # 1st model without adjustements
  
RF_metabolites <- subset(new_data2[,-c(1,2)], select=as.character(colnames(new_data2[,-c(1,2)])))
responses.factor_w <-factor(new_data2$New_Group)
# p.glm <-NULL
# odds.CI <- NULL
# 
# for(i in 1:length(RF_metabolites)){
#   
#   fit <- glm(responses.factor_w~RF_metabolites[,i],family=binomial(link='logit'))
#   p.glm[i] <- summary(fit)$coefficients[2, 4]
#   tmp <-data.frame(Metabolite=colnames(RF_metabolites)[i],exp(cbind(OR = coef(fit), confint(fit))))
#   odds.CI <- rbind(odds.CI, tmp)
#   
# }
# 
# 
# adj.pvalues_healthy_vs_western <- p.adjust(p.glm, method = "fdr", n = length(p.glm))
# Healthy_Western_model_all_pvalues <- cbind.data.frame(colnames(new_data2[,-c(1,2)]),p.glm,adj.pvalues_healthy_vs_western)
# write.csv(Healthy_Western_model_all_pvalues, file = paste0(output.path, "model/","Baltica_logisticregression_pvalues_healthyvswestern_LOG_NORM_17062019.csv"), row.names = FALSE, na = "")
# write.csv(odds.CI, file = paste0(output.path, "model/","Baltica_logisticregression_OR_healthyvswestern_LOG_NORM_17062019.csv"), row.names = FALSE, na = "")
# 
# 
# summary(fit) # display results
# confint(fit) # 95% CI for the coefficients
# exp(coef(fit)) # exponentiated coefficients
# exp(confint(fit)) # 95% CI for exponentiated coefficients
# predict(fit, type="response") # predicted values
# residuals(fit, type="deviance") # residuals
# 
# # 2st adjustment logistic regression  model for age, waist to hip ratio, physical activity, alcohol consumption, total cholesterol
# 
# 
# p.glm2 <-NULL
# odds.CI2 <- NULL
# 
# for(i in 1:length(RF_metabolites)){
#   
#   fit <- glm(responses.factor_w~RF_metabolites[,i]+final_covariates$`Age years`+final_covariates$`S-Tot-chol mmol/l`+final_covariates$`Alcohol /week g (JK)`+final_covariates$WaistHipRatio+final_covariates$`Ener. exp. CLTPA (LIMETC2/365xp) kcal/d)`,family=binomial(link='logit'))
#   p.glm2[i] <- summary(fit)$coefficients[2, 4]
#   tmp <-data.frame(Metabolite=colnames(RF_metabolites)[i],exp(cbind(OR = coef(fit), confint(fit))))
#   odds.CI2 <- rbind(odds.CI2, tmp)
# }
# 
# adj.pvalues_healthy_vs_western2 <- p.adjust(p.glm2, method = "fdr", n = length(p.glm2))
# Healthy_Western_model_all_pvalues2 <- cbind.data.frame(rownames(allMetabolites_rank),p.glm2,adj.pvalues_healthy_vs_western2)
# write.csv(Healthy_Western_model_all_pvalues2, file = paste0(output.path, "model/","Baltica_logisticregression_pvalues_healthyvswestern_model2_16042019.csv"), row.names = FALSE, na = "")
# write.csv(odds.CI2, file = paste0(output.path, "model/","Baltica_logisticregression_OR_healthyvswestern_model2_16042019.csv"), row.names = FALSE, na = "")





# Western Versus Healthy linear regression 

# 1st model without adjustements

# Control Versus Healthy 

# 1st model without adjustements

RF_metabolites <- subset(new_data2, select=as.character(colnames(new_data2[,-c(1,2)])))
responses.factor_w <-factor(new_data2$New_Group)
p.lm <-NULL
lm.CI <- NULL
estimate1 <-NULL
for(i in 1:length(RF_metabolites)){
  
  fit <- lm(RF_metabolites[,i]~responses.factor_w)
  p.lm[i] <- anova(fit)$'Pr(>F)'[1]
  #lm.CI[i] <-confint(fit, 'responses.factor_w1',level=0.95)
  estimate1[i] <- summary(fit)$coefficients[2,1]
  tmp <-data.frame(Metabolite=colnames(RF_metabolites)[i], p.lm[i] ,estimate1[i],confint(fit, 'responses.factor_w1',level=0.95))
  lm.CI <- rbind(lm.CI, tmp)

}


adj.pvalues_healthy_vs_western <- p.adjust(p.lm, method = "fdr", n = length(p.lm))
Healthy_Western_model_all_pvalues <- cbind.data.frame(lm.CI,adj.pvalues_healthy_vs_western)
write.csv(Healthy_Western_model_all_pvalues, file = paste0(output.path, "model/","Baltica_linearregression_pvaluesCI_healthyvswestern_26062019.csv"), row.names = FALSE, na = "")

# 2st adjustment logistic regression  model for age, waist to hip ratio, physical activity, alcohol consumption, total cholesterol


p.lm2 <-NULL
lm.CI2 <- NULL
estimate2 <-NULL
for(i in 1:length(RF_metabolites)){
  
  fit <- lm(RF_metabolites[,i]~responses.factor_w+final_covariates$`Age years`+final_covariates$`S-Tot-chol mmol/l`+final_covariates$`Alcohol /week g (JK)`+final_covariates$WaistHipRatio+final_covariates$`Ener. exp. CLTPA (LIMETC2/365xp) kcal/d)`)
  p.lm2[i] <- summary(fit)$coefficients[2, 4]
  estimate2[i] <- summary(fit)$coefficients[2,1]
  tmp <-data.frame(Metabolite=colnames(RF_metabolites)[i],estimate2[i],p.lm2[i],confint(fit,'responses.factor_w1',level=0.95))
  lm.CI2 <- rbind(lm.CI2, tmp)
}

adj.pvalues_healthy_vs_western2 <- p.adjust(lm.CI2$p.lm2, method = "fdr", n = length(p.lm2))
Healthy_Western_model_all_pvalues2 <- cbind.data.frame(lm.CI2,adj.pvalues_healthy_vs_western2)
write.csv(Healthy_Western_model_all_pvalues2, file = paste0(output.path, "model/","Baltica_linearregression_pvaluesCI_healthyvswestern_26062019_model2.csv"), row.names = FALSE, na = "")



