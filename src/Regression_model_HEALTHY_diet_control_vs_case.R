
#### SECOND STEP OF THE ANALYSIS CONTROL VERSUS HEALTHY #### 
### FILTER THE BAD SIGNALS 

###  Perform anaova analysis and exclude everything above and equal p.values =>0.9

# Side step 

# We have four Groups: 

# 1: control-high score

# 2: case-high score

# 3: case-low score 

# 4: control-low score

#  and I will  join the 1-4 and 2-3 groups together 

# In thid first step I am interested control versus cases



rm(list = ls())
output.path <- "~/Documents/PhD_2016-2019/projects//Baltica_Folder/Baltica/Baltica_project/results/"

# load the data
load("~/Documents/PhD_2016-2019/projects/Baltica_Folder/Baltica/Baltica_project/results/data/BALTICA_rawdata_18052019.rda")
load("~/Documents/PhD_2016-2019/projects/Baltica_Folder/Baltica/Baltica_project/results/data/CLEAN_DRIFT_allmodes.Rdata")
load("~/Documents/PhD_2016-2019/projects/Baltica_Folder/Baltica/Baltica_project/results/data/covariates_data.rda")
load("~/Documents/PhD_2016-2019/projects/Baltica_Folder/Baltica/Baltica_project/results/data/full_compound_data.Rdata")
load("~/Documents/PhD_2016-2019/projects/Baltica_Folder/Baltica/Baltica_project/results/model/MetaboliteRanks_MUVRmodel16042019.rda")

# combine the pareto scale data with the GROUP vector 

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

#  the controls with high  together using the  subset function 
Controls_H <-  subset(new_data, GROUP=="4")


# combine all the cases together using the  subset function 

Cases_H <-  subset(new_data, GROUP=="3" )



# Create new labels 0 for controls and 1 for cases
dim(Controls_H)
dim(Cases_H)
Control_label <-t(rep(0,dim(Controls_H)[1]))
Case_label <- t(rep(1,dim(Cases_H)[1]))
New_Group <-c(Control_label,Case_label)

# combine Case and control data with the new labels

new_data2 <- rbind.data.frame(Controls_H,Cases_H)
new_data2 <-cbind.data.frame(New_Group,new_data2)

#################  REGRESSION MODELS ###############

GROUP <- clean.data.all.modes.no.qcs$GROUP

new.covariates.data <- cbind.data.frame(GROUP,covariates.data)

# combine all the controls together using the  subset function 
Controls <-  subset(new.covariates.data, GROUP=="4" )


# combine all the cases together using the  subset function 

Cases <-  subset(new.covariates.data, GROUP=="3"  )


final_covariates <- rbind.data.frame(Controls,Cases)
#save(round_allMetabolites_rank,allMetabolites_rank,listMetabolites,responses.factor,final_covariates,file = paste0(output.path, "model/",  "Baltica_Metabolites_control_case21032019.rda"))

# Control Versus Healthy 

# 1st model without adjustements

#save(responses,responses.factor,allMetabolites_rank,listMetabolites,file = paste0(output.path,"data/","list_metabolites_control_vs_case_MUVRmodel.rda"))

RF_metabolites <- subset(new_data2, select=as.character(rownames(allMetabolites_rank)))
responses.factor_c <-factor(new_data2$New_Group)


p.glm <-NULL
odds.CI <-NULL
estimate1 <-NULL
CI1 <- NULL

for(i in 1:length(RF_metabolites)){
  
  fit <- glm(responses.factor_c~RF_metabolites[,i],family=binomial(link='logit'))
  p.glm[i] <- summary(fit)$coefficients[2, 4]
  estimate1[i] <- summary(fit)$coefficients[2,1]
  CI1 <- confint(fit)[2,]
  tmp <-data.frame(Metabolite=colnames(RF_metabolites)[i],exp(cbind(OR = estimate1[i], CI=t(CI1))))
  odds.CI <- rbind(odds.CI, tmp)
  
  
}


adj.pvalues_Control_vs_case <- p.adjust(p.glm, method = "fdr", n = length(p.glm))
Control_Case_model_all_pvalues <- cbind.data.frame(rownames(allMetabolites_rank),p.glm,adj.pvalues_Control_vs_case)
all_results <-c(Control_Case_model_all_pvalues,odds.CI)

write.csv(all_results, file = paste0(output.path, "model/","ALL_Results_Baltica_logisticregression_pvalues_HEALTHYDIET_controlvscase_26062019.csv"), row.names = FALSE, na = "")

#write.csv(Cont
#write.csv(Control_Case_model_all_pvalues, file = paste0(output.path, "model/","Baltica_logisticregression_pvalues_HEALTHYDIET_controlvscase_16042019.csv"), row.names = FALSE, na = "")
#write.csv(odds.CI, file = paste0(output.path, "model/","Baltica_logisticregression_OR_HEALTHYDIET_controlvscase_16042019.csv"), row.names = FALSE, na = "")


summary(fit) # display results
confint(fit) # 95% CI for the coefficients
exp(coef(fit)) # exponentiated coefficients
exp(confint(fit)) # 95% CI for exponentiated coefficients
predict(fit, type="response") # predicted values
residuals(fit, type="deviance") # residuals

