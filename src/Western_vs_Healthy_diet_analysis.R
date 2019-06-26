
rm(list = ls())

# list of libraries
library(doParallel)     # Parallel processing
library(dplyr)
library(doMC)
library(MUVR)           # Multivariate modelling
library(mdatools)
require(corrplot)
require(car) 


# Path settings
output.path <- "~/Documents/PhD_2016-2019/projects//Baltica_Folder/Baltica/Baltica_project/results/"

# load data
load("~/Documents/material_PhD/projects/baltica_preliminary/Baltica/Baltica_project/results/data/CLEAN_DRIFT_allmodes.Rdata")

###### FIRST PART OF THE PIPELINE : THE PREPROCESSING STEP

#### 1. LOG TRANSFORMATION
#### 2. Normalization
#### 3. PARETO SCALED

# extract the data by removing the first two columns which contain the GROUP info and the RUN_ORDER

baltica <- clean.data.all.modes.no.qcs[,-c(1,2)]
# put in a seperate vectror the group info 
GROUP <- clean.data.all.modes.no.qcs$GROUP


#
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


## NB quicky check the data by doing PCA

## 
# 
# m2 = mdatools::pca(data_norm, 7, scale = T, info = "PCA modelfor quality control")
# m2 = selectCompNum(m2, 5)
# print(m2)
# 
# # PCA plots 
# par(mfrow = c(1, 1))
# plotScores(m2, c(1, 2), cgroup = GROUP, show.labels = F)
# summary(m)




#### SECOND STEP OF THE ANALYSIS CONTROL VERSUS HEALTHY #### 
### FILTER THE BAD SIGNALS 

###  Perform anaova analysis and exclude everything above and equal p.values =>0.9

# Side step 

# We have four Groups: 

# 1: control-high score

# 2: case-high score

# 3: case-low score 

# 4: control-low score

#  and I will  join the 1-2 and 3-4 groups together 

# In thid first step I am interested control versus cases

# combine the pareto scale data with the GROUP vector 

new_data <- cbind.data.frame(GROUP,data_norm)

# combine all the controls together using the  subset function 
Western <-  subset(new_data, GROUP=="1" | GROUP=="2" )


# combine all the cases together using the  subset function 

Healthy <-  subset(new_data, GROUP=="3" | GROUP=="4" )



# Create new labels 0 for controls and 1 for cases

Healthy_label <-t(rep(0,182))
Western_label <- t(rep(1,182))
New_Group <-c(Healthy_label,Western_label)

# combine Case and control data with the new labels

new_data2 <- rbind.data.frame(Healthy,Western)
new_data2 <-cbind.data.frame(New_Group,new_data2)

# mean center
mu <-apply(new_data2[,-c(1,2)],2,mean)

# pareto scaled the data 
library(MetabolAnalyze)
pareto_scale <- scaling(new_data2[,-c(1,2)], type = "pareto")

# Perform anova in the new_data2 as a prefiltering procedure

# Here I will use nonparametric anova because the metabolites are abit skewed.

# Kruskal–Wallis Test

# http://rcompanion.org/handbook/F_08.html
# https://rcompanion.org/rcompanion/d_06.html
pvalues_BALTICA <- NULL
# remove the columns of GROUP and RUN_ORDER

tmp.data <- pareto_scale

# some statistics

# multiple testing with wilcoxon test and ttest 
for(i in 1:dim(tmp.data)[2]) {
  result_ANOVA <- kruskal.test(tmp.data[,i]~ New_Group)
  pvalues_BALTICA[i] <-result_ANOVA$p.value
  
}

# FDR correction for control versus helathy from wilcoxon test 

adj.pvalues_BALTICA <- p.adjust(pvalues_BALTICA, method = "fdr", n = length(pvalues_BALTICA))

responses <- colnames(tmp.data)
View(responses)
# combine p values and and adjusted pvalues into one dataframe
BALTICA_all_pvalues <- cbind.data.frame(responses,pvalues_BALTICA,adj.pvalues_BALTICA)

# select the ones that are below the significant level 0.05
index <- which(BALTICA_all_pvalues$adj.pvalues_BALTICA <= 0.9)
smallpvalues <-BALTICA_all_pvalues[index,]

smallpvalues_WC <- smallpvalues
write.csv(smallpvalues_WC, file = paste0(output.path, "data/","Baltica_one_way_anova_prefiltering_pvalues_healthy_vesrus_western.csv"), row.names = FALSE, na = "")
write.csv(BALTICA_all_pvalues, file = paste0(output.path, "data/","Baltica_one_way_anova_pvalues_healthy_vesrus_western.csv"), row.names = FALSE, na = "")

save(smallpvalues_WC,file = paste0(output.path, "data/","BALTICA_smallpvalues_westen_vesrus_healthy.rda"))

# data with the most significant metabolites from pvalues
data.prefiltering <- subset(tmp.data, select=as.character(smallpvalues_WC$responses))
datat <- t.data.frame(data.prefiltering)
write.csv(datat, file = paste0(output.path, "data/","Dataprefiltering_WC.csv"), row.names = TRUE, na = "")



# a good site for multiple comparisons 
# http://www.nonlinear.com/support/progenesis/comet/faq/v2.0/pq-values.aspx


# m2 = mdatools::pca(data.prefiltering, 7, scale = T, info = "PCA modelfor quality control")
# m2 = selectCompNum(m2, 5)
# print(m2)
# 
# # PCA plots 
# par(mfrow = c(1, 1))
# plotScores(m2, c(1, 2), cgroup = New_Group, show.labels = F)
# summary(m)

#### THIRD STEP OF THE ANALYSIS #### 

# Perform random forest Classification trees

responses.factor  <- factor(New_Group)
responses <-colnames(data.prefiltering)
save(data.prefiltering,responses.factor,responses,file = paste0(output.path, "data/","BALTICA_MUVRdata_22032019_healthy_versus_western.rda"))

#  ANALYSIS
# # Set method parameters
nCore=detectCores()-1   # Number of processor threads to use
nRep=2*nCore   # Number of repetitions per actual model and permutations
nOuter=6       # Number of validation segments
varRatio=0.6  # Proportion of variables to keep per iteration during variable selection
method='RF'    # Core modelling technique
model=1        # 1 for min, 2 for mid and 3 for max
nPerm=25
permFit=numeric(nPerm)
#

# 
# Set up parallel processing using doParallel 
cl=makeCluster(nCore)
registerDoParallel(cl)
# Perform modelling

classModel = MUVR(X=data.prefiltering, Y=responses.factor, nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method,scale = F)
# # Stop parallel processing
stopCluster(cl)
save (classModel,paste0(output.path, "model/", "BALTICA_MUVR_model_14032019_healthy_vs_western.rda"))
# Compute permuted models and extract fitness metrics; Approx 12 mins
YPerm <- sample(responses.factor)



for (p in 1:nPerm) {
  cat('\nPermutation',p,'of',nPerm)
  YPerm=sample(responses.factor)
  perm=MUVR(X=data.prefiltering,Y=YPerm,nRep=nRep,nOuter=nOuter,varRatio=varRatio,method=method,scale = F)
  permFit[p]=perm$miss[model]
}
stopCluster(cl)

save(perm,file = paste0(output.path, "model/","full.permutation_healthy_vs_western.rda"))




# Examine model performance and output
classModel$miss                   # Number of misclassifications for min, mid and max models

classModel$nVar                   # Number of variables for min, mid and max models
#
cbind(responses.factor, classModel$yClass)    # Actual class side-by-side with min, mid and max predictions



confusionMatrix=function(MVObj,model='mid'){
  if(!any(class(MVObj)=='Classification')) stop ('The MUVR object needs to be from a classification analysis')
  nMod=ifelse(model=='min',1,ifelse(model=='mid',2,3))
  actual=MVObj$inData$Y
  predicted=MVObj$yClass[,nMod]
  table(actual=actual,predicted=predicted)
}
confusionMatrix(classModel,model='max')

# …  …    …   …   …   …  Plot MODEL
plotVAL(classModel) 
plotMV(classModel, model='max')         # Look at the model of choice: min, mid or max
plotStability(classModel, model='max') 
plotVIP(classModel, model='max')  

listMetabolites <-getVIP(classModel, model='max')
write.csv(listMetabolites,  file = paste0(output.path, "model/",  "Baltica_VIP_maxmodel_healthy_western.csv"), row.names = FALSE, na = "")
allMetabolites_rank <- classModel$VIP
round_allMetabolites_rank <- round(allMetabolites_rank,digits = 2)
write.csv(round_allMetabolites_rank, file =paste0(output.path, "model/",  "Baltica_allMetabolites_rank_healthy_vs_western_21032019.csv"), row.names = TRUE, na = "")



#################  REGRESSION MODELS ###############

rm(list = ls())
#load data 
load("~/Documents/material_PhD/projects/baltica_preliminary/Baltica/Baltica_project/results/data/covariates_data.rda")
load("~/Documents/material_PhD/projects/baltica_preliminary/Baltica/Baltica_project/results/data/CLEAN_DRIFT_allmodes.Rdata")
load("~/Documents/material_PhD/projects/baltica_preliminary/Baltica/Baltica_project/results/data/BALTICA_MUVRdata_22032019_healthy_versus_western.rda")

load("~/Documents/material_PhD/projects/baltica_preliminary/Baltica/Baltica_project/results/model/BALTICA_MUVR_model_22032019_healthy_vs_western.rda")

save(round_allMetabolites_rank,responses,responses.factor,allMetabolites_rank,listMetabolites,file = paste0(output.path,"data/","list_metabolites_healthy_vs_western_MUVRmodel.rda"))

# extract the data by removing the first two columns which contain the GROUP info and the RUN_ORDER


# put in a seperate vectror the group info 
GROUP <- clean.data.all.modes.no.qcs$GROUP


new.covariates.data <- cbind.data.frame(GROUP,covariates.data)

# combine all the controls together using the  subset function 
Controls <-  subset(new.covariates.data, GROUP=="1" | GROUP=="2" )


# combine all the cases together using the  subset function 

Cases <-  subset(new.covariates.data, GROUP=="3" | GROUP=="4" )


final_covariates <- rbind.data.frame(Controls,Cases)

# Control Versus Healthy 

# 1st model without adjustements

RF_metabolites <- subset(data.prefiltering, select=as.character(listMetabolites$name))

p.glm <-NULL
for(i in 1:length(RF_metabolites)){
  
  fit <- glm(responses.factor~RF_metabolites[,i],family=binomial())
  p.glm[i] <- summary(fit)$coefficients[2, 4]
}


adj.pvalues_healthy_vs_western <- p.adjust(p.glm, method = "fdr", n = length(p.glm))
Healthy_Western_model_all_pvalues <- cbind.data.frame(listMetabolites$name,p.glm,adj.pvalues_healthy_vs_western)
write.csv(Healthy_Western_model_all_pvalues, file = paste0(output.path, "model/","Baltica_logisticregression_pvalues_healthyvswestern.csv"), row.names = FALSE, na = "")


summary(fit) # display results
confint(fit) # 95% CI for the coefficients
exp(coef(fit)) # exponentiated coefficients
exp(confint(fit)) # 95% CI for exponentiated coefficients
predict(fit, type="response") # predicted values
residuals(fit, type="deviance") # residuals

# 2st adjustment logistic regression  model for age, waist to hip ratio, physical activity, alcohol consumption, total cholesterol


RF_metabolites <- data.prefiltering[,listMetabolites$name]
p.glm2 <-NULL
for(i in 1:length(RF_metabolites)){
  
  fit <- glm(responses.factor~RF_metabolites[,i]+final_covariates$`Age years`+final_covariates$`S-Tot-chol mmol/l`+final_covariates$`Alcohol /week g (JK)`+final_covariates$WaistHipRatio+final_covariates$`Ener. exp. CLTPA (LIMETC2/365xp) kcal/d)`,family=binomial())
  p.glm2[i] <- summary(fit)$coefficients[2, 4]
}


adj.pvalues_healthy_vs_western2 <- p.adjust(p.glm2, method = "fdr", n = length(p.glm2))
Healthy_Western_model_all_pvalues2 <- cbind.data.frame(listMetabolites$name,p.glm2,adj.pvalues_healthy_vs_western2)
write.csv(Healthy_Western_model_all_pvalues, file = paste0(output.path, "model/","Baltica_logisticregression_pvalues_healthyvswestern_model2.csv"), row.names = FALSE, na = "")

