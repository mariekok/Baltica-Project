
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
load("~/Documents/PhD_2016-2019/projects//Baltica_Folder/Baltica/Baltica_project/results/data/CLEAN_DRIFT_allmodes.Rdata")

###### FIRST PART OF THE PIPELINE : THE PREPROCESSING STEP

#### 1. LOG TRANSFORMATION
#### 2. Normalization
#### 3. PARETO SCALED

# extract the data by removing the first two columns which contain the GROUP info and the RUN_ORDER

baltica <- clean.data.all.modes.no.qcs[,-c(1,2)]
# put in a seperate vectror the group info 
GROUP <- clean.data.all.modes.no.qcs$GROUP
responses.factor  <- factor(GROUP)

responses <-colnames(baltica)

save(responses.factor,baltica,GROUP,responses,file = paste0(output.path, "data/","BALTICA_rawdata_18052019.rda"))


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

# mean center
mu <-apply(data_norm,2,mean)

# pareto scaled the data 
library(MetabolAnalyze)
pareto_scale <- scaling(data_norm, type = "pareto")

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




#### SECOND STEP OF THE ANALYSIS RF #### 


# Side step 

# We have four Groups: 

# 1: control-high score

# 2: case-high score

# 3: case-low score 

# 4: control-low score

#  and I will  join the 1-2 and 3-4 groups together 




# Perform random forest Classification trees

responses.factor  <- factor(GROUP)
responses <-colnames(data_norm)
save(responses.factor,pareto_scale,data_norm,GROUP,baltica,responses,file = paste0(output.path, "data/","BALTICA_MUVRdata_16042019.rda"))
# 
# #  ANALYSIS
# # # Set method parameters
# nCore=detectCores()-1   # Number of processor threads to use
# nRep=2*nCore   # Number of repetitions per actual model and permutations
# nOuter=6       # Number of validation segments
# varRatio=0.6  # Proportion of variables to keep per iteration during variable selection
# method='RF'    # Core modelling technique
# model=1        # 1 for min, 2 for mid and 3 for max
# nPerm=25
# permFit=numeric(nPerm)
# #
# 
# # 
# # Set up parallel processing using doParallel 
# cl=makeCluster(nCore)
# registerDoParallel(cl)
# # Perform modelling
# 
# classModel = MUVR(X=data.prefiltering, Y=responses.factor, nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method,scale = F)
# # # Stop parallel processing
# stopCluster(cl)
# save (classModel,paste0(output.path, "model/", "BALTICA_MUVR_model_14032019_healthy_vs_western.rda"))
# # Compute permuted models and extract fitness metrics; Approx 12 mins
# YPerm <- sample(responses.factor)
# 
# 
# 
# for (p in 1:nPerm) {
#   cat('\nPermutation',p,'of',nPerm)
#   YPerm=sample(responses.factor)
#   perm=MUVR(X=data.prefiltering,Y=YPerm,nRep=nRep,nOuter=nOuter,varRatio=varRatio,method=method,scale = F)
#   permFit[p]=perm$miss[model]
# }
# stopCluster(cl)
# 
# save(perm,file = paste0(output.path, "model/","full.permutation_healthy_vs_western.rda"))

load("~/Documents/PhD_2016-2019/projects/Baltica_Folder/Baltica/Baltica_project/results/model/BALTICA_16042019_paretoscaled.rda")


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
write.csv(listMetabolites,  file = paste0(output.path, "model/",  "Baltica_VIP_maxmodel_4groups.csv"), row.names = FALSE, na = "")
allMetabolites_rank <- classModel$VIP
round_allMetabolites_rank <- round(allMetabolites_rank,digits = 4)

write.csv(round_allMetabolites_rank, file =paste0(output.path, "model/",  "Baltica_allMetabolites_rank_4groups.csv"), row.names = TRUE, na = "")
save(allMetabolites_rank,listMetabolites,round_allMetabolites_rank,file = paste0(output.path, "model/","MetaboliteRanks_MUVRmodel16042019.rda"))
