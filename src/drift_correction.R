
rm(list=ls())

# list of libraries
library(doMC)
library(dplyr)
library(doParallel)     # Parallel processing
library(MUVR)
library(dplyr)
library(doMC)
library(MUVR)           # Multivariate modelling
library(mdatools)
require(corrplot)
require(car) 
require(glmnet)
library(readr)
library(tidyr)

# Path settings
input.path <- "~/Documents/PhD_2016-2019/projects/Baltica_Folder/Baltica/data/"
output.path <- "~/Documents/PhD_2016-2019/projects//Baltica_Folder/Baltica/Baltica_project/results/"

# Load external functions
# contains the vizualization functions
source("~/Documents/PhD_2016-2019/projects/visualizations/src/functions.R")
source("~/Documents/PhD_2016-2019/projects/visualizations/src/visualization_functions.R")
# contains the drift correction and quality check functions
source("~/Documents/PhD_2016-2019/projects/MetaboQC-master/src/new_functions.R")


# load data
load("~/Documents/PhD_2016-2019/projects/Baltica_Folder/Baltica/Baltica_project/results/data/full_compound_data.Rdata")
# How many traits will be analyzed, NULL for all
n.responses <- NULL

# Should everything be run from the scratch
initialize <- TRUE
# Log settings
logging <- TRUE
log.file <- paste(output.path, "analysis_log.txt", sep="")
init.log(log.file=log.file, logging=logging)


# hilic negative  # hilic negative  # hilic negative  # hilic negative  # hilic negative  # hilic negative  # hilic negative
# the metabolites in one mode
original.responses <- select.responses(data=combine.hilic.neg, pattern = "^COMPOUND_", n=n.responses)
responses <- original.responses
# start quality check
qm <- compute.quality.measures(data = combine.hilic.neg, responses = responses, group.col = "GROUP")
qm$Stage <- "Original"
qm.full <-data.frame()
qm.full <- rbind(qm.full,qm)
viz.output.path <- paste(output.path, "figures/", "fig_hilicneg", "_", sep="")


# drift correction 
drift.correction.results <- correct.drift.spline(data = combine.hilic.neg, responses = responses, group.col = "GROUP", order.col = "RUN_ORDER", smooth.par.low = 0.5,
                                                  plotting = FALSE, plot.path = paste0(viz.output.path, "drift_correction.pdf"), log.file = log.file, logging = logging)
# extract drift corrected data
drift.corrected.data.hilic.neg <- drift.correction.results$drift.corrected.data
operations <- drift.correction.results$operations

qm <- compute.quality.measures(data = drift.corrected.data.hilic.neg, responses, group.col = "GROUP")
qm$Stage <- "Drift_corrected"
qm.full <- rbind(qm.full,qm)
#
# # seperate good responses  from low quality signals 
response.stage <- data.frame(response = responses, stage = "kept", stringsAsFactors = FALSE)
low.quality.responses <- operations[!operations$kept,"response"]
responses <- setdiff(responses, low.quality.responses)

cleaned.data.hilic.neg <- drift.corrected.data.hilic.neg[,responses]
cleaned.data.hilic.neg <- cbind.data.frame(combine.hilic.neg[,1:2],cleaned.data.hilic.neg)

qm <- compute.quality.measures(data = cleaned.data.hilic.neg, responses, group.col = "GROUP")
qm$Stage <- "Cleaned"
qm.full <- rbind(qm.full,qm)
#save.quality.measure.summary.plot(qm.full, paste0(output.path, "figures/Quality_summary_", "full_compound_data_hilic_neg", ".pdf"))

no.qc.data.hilic.neg <- droplevels(cleaned.data.hilic.neg[(cleaned.data.hilic.neg$GROUP != "QC"),])

#########################################################################################################################################################

# hilic positive # hilic positive # hilic positive  # hilic positive # hilic positive # hilic positive # hilic positive # hilic positive # hilic positive
# the metabolites in one mode
original.responses <- select.responses(data=combine.hilic.pos, pattern = "^COMPOUND_", n=n.responses)
responses <- original.responses
# start quality check
qm <- compute.quality.measures(data = combine.hilic.pos, responses = responses, group.col = "GROUP")
qm$Stage <- "Original"
qm.full <-data.frame()
qm.full <- rbind(qm.full,qm)
viz.output.path <- paste(output.path, "figures/", "fig_hilicpos", "_", sep="")


# drift correction 
drift.correction.results <- correct.drift.spline(data = combine.hilic.pos, responses = responses, group.col = "GROUP", order.col = "RUN_ORDER", smooth.par.low = 0.5,
                                                  plotting = FALSE, plot.path = paste0(viz.output.path, "drift_correction.pdf"), log.file = log.file, logging = logging)



# extract drift corrected data
drift.corrected.data.hilic.pos <- drift.correction.results$drift.corrected.data
operations <- drift.correction.results$operations
#compute quality measures for dfit corrected
qm <- compute.quality.measures(data = drift.corrected.data.hilic.pos, responses, group.col = "GROUP")
qm$Stage <- "Drift_corrected"
qm.full <- rbind(qm.full,qm)

# keep the good signals
response.stage <- data.frame(response = responses, stage = "kept", stringsAsFactors = FALSE)
low.quality.responses <- operations[!operations$kept,"response"]
responses <- setdiff(responses, low.quality.responses)

# data cleaned from low quality signals
cleaned.data.hilic.pos <- drift.corrected.data.hilic.pos[,responses]
cleaned.data.hilic.pos <- cbind.data.frame(combine.hilic.pos[,1:2],cleaned.data.hilic.pos)
qm <- compute.quality.measures(data = cleaned.data.hilic.pos, responses, group.col = "GROUP")
qm$Stage <- "Cleaned"
qm.full <- rbind(qm.full,qm)

# data cleaned and no qcs
no.qc.data.hilic.pos <- droplevels(cleaned.data.hilic.pos[(cleaned.data.hilic.pos$GROUP != "QC"),])


#########################################################################################################################################################

# rp positive # rp positive # rp positive  # rp positive # rp positive # rp positive # rp positive # rp positive # rp positive
# the metabolites in one mode
original.responses <- select.responses(data=combine.rp.pos, pattern = "^COMPOUND_", n=n.responses)
responses <- original.responses
# start quality check
qm <- compute.quality.measures(data = combine.rp.pos, responses = responses, group.col = "GROUP")
qm$Stage <- "Original"
qm.full <-data.frame()
qm.full <- rbind(qm.full,qm)
viz.output.path <- paste(output.path, "figures/", "fig_rppos", "_", sep="")


# drift correction 
drift.correction.results <- correct.drift.spline(data = combine.rp.pos, responses = responses, group.col = "GROUP", order.col = "RUN_ORDER", smooth.par.low = 0.5,
                                                  plotting = FALSE, plot.path = paste0(viz.output.path, "drift_correction.pdf"), log.file = log.file, logging = logging)



# extract drift corrected data
drift.corrected.data.rp.pos <- drift.correction.results$drift.corrected.data
operations <- drift.correction.results$operations
# write drift correct data in a csv file

qm <- compute.quality.measures(data = drift.corrected.data.rp.pos, responses, group.col = "GROUP")
qm$Stage <- "Drift_corrected"
qm.full <- rbind(qm.full,qm)

# keep the good signals
response.stage <- data.frame(response = responses, stage = "kept", stringsAsFactors = FALSE)
low.quality.responses <- operations[!operations$kept,"response"]
responses <- setdiff(responses, low.quality.responses)
# data cleaned from low quality signals
cleaned.data.rp.pos <- drift.corrected.data.rp.pos[,responses]
cleaned.data.rp.pos <- cbind.data.frame(combine.rp.pos[,1:2],cleaned.data.rp.pos)
qm <- compute.quality.measures(data = cleaned.data.rp.pos, responses, group.col = "GROUP")
qm$Stage <- "Cleaned"
qm.full <- rbind(qm.full,qm)
#save.quality.measure.summary.plot(qm.full, paste0(output.path, "figures/Quality_summary_", "full_compound_data_CSF_rp_pos", ".pdf"))

# data cleaned and no qcs
no.qc.data.rp.pos <- droplevels(cleaned.data.rp.pos[(cleaned.data.rp.pos$GROUP != "QC"),])
#########################################################################################################################################################

# rp negative # rp negative # rp negative  # rp negative # rp negative # rp negative # rp negative # rp negative # rp negative
# the metabolites in one mode
original.responses <- select.responses(data=combine.rp.neg, pattern = "^COMPOUND_", n=n.responses)
responses <- original.responses
# start quality check
qm <- compute.quality.measures(data = combine.rp.neg, responses = responses, group.col = "GROUP")
qm$Stage <- "Original"
qm.full <-data.frame()
qm.full <- rbind(qm.full,qm)
viz.output.path <- paste(output.path, "figures/", "fig_rpneg", "_", sep="")


# drift correction 
drift.correction.results <- correct.drift.spline(data = combine.rp.neg, responses = responses, group.col = "GROUP", order.col = "RUN_ORDER", smooth.par.low = 0.5,
                                                  plotting = FALSE, plot.path = paste0(viz.output.path, "drift_correction.pdf"), log.file = log.file, logging = logging)
# extract drift corrected data
drift.corrected.data.rp.neg <- drift.correction.results$drift.corrected.data
operations <- drift.correction.results$operations

# quality control  in drift corrected
qm <- compute.quality.measures(data = drift.corrected.data.rp.neg, responses, group.col = "GROUP")
qm$Stage <- "Drift_corrected"
qm.full <- rbind(qm.full,qm)

# keep the good signals
response.stage <- data.frame(response = responses, stage = "kept", stringsAsFactors = FALSE)
low.quality.responses <- operations[!operations$kept,"response"]
responses <- setdiff(responses, low.quality.responses)

# data cleaned from low quality signals
cleaned.data.rp.neg <- drift.corrected.data.rp.neg[,responses]
cleaned.data.rp.neg <- cbind.data.frame(combine.rp.neg[,1:2],cleaned.data.rp.neg)
qm <- compute.quality.measures(data = cleaned.data.rp.neg, responses, group.col = "GROUP")
qm$Stage <- "Cleaned"
qm.full <- rbind(qm.full,qm)
#save.quality.measure.summary.plot(qm.full, paste0(output.path, "figures/Quality_summary_", "full_compound_data_CSF_rp_neg", ".pdf"))

# data cleaned and no qcs
no.qc.data.rp.neg <- droplevels(cleaned.data.rp.neg[(cleaned.data.rp.neg$GROUP != "QC"),])

# combine all modes together 
clean.data.all.modes.no.qcs <- cbind.data.frame(no.qc.data.hilic.neg,no.qc.data.hilic.pos[,-c(1,2)],no.qc.data.rp.neg[,-c(1,2)],no.qc.data.rp.pos[,-c(1,2)])
# save data
clean.data.all.modes.with.qcs <- cbind.data.frame(cleaned.data.hilic.neg,cleaned.data.hilic.pos[,-c(1,2)],cleaned.data.rp.neg[,-c(1,2)],cleaned.data.rp.pos[,-c(1,2)])



save(clean.data.all.modes.with.qcs,clean.data.all.modes.no.qcs,no.qc.data.hilic.neg,no.qc.data.hilic.pos,no.qc.data.rp.neg,no.qc.data.rp.pos,group.info,file = paste0(output.path, "data/","CLEAN_DRIFT_allmodes.Rdata"))
# 
