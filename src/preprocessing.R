
rm(list=ls())

library(doMC)
library(dplyr)

# Path settings
# Path settings
input.path <- "~/Documents/PhD_2016-2019/projects/Baltica_Folder/Baltica/Baltica_project/data/"
output.path <- "~/Documents/PhD_2016-2019/projects/Baltica_Folder/Baltica/Baltica_project/results/"

# Load external functions
source("~/Documents/PhD_2016-2019/projects/visualizations/src/functions.R")
source("~/Documents/PhD_2016-2019/projects/visualizations/src/visualization_functions.R")
# How many traits will be analyzed, NULL for all
n.responses <- NULL

# Should everything be run from the scratch
initialize <- TRUE
# Log settings
logging <- TRUE
log.file <- paste(output.path, "analysis_log.txt", sep="")
init.log(log.file=log.file, logging=logging)

# Parallel processing settings
#registerDoMC(cores=30)



# load the data 
library(readr)
baltica <- read_csv(paste0(input.path,"20181121_BALTICA_newalignment_raw_filtered_forstat.csv"))
group.info <- read_delim("Baltica/Baltica_project/data/20181120_BALTICA_groupinfo_raw_forstat (1).csv", 
                        ";", escape_double = FALSE, trim_ws = TRUE)
GROUP <- group.info$Group
RUN_ORDER <- group.info$DataFileNo


# divide the data into 4 modes and per different part of the brain 

## Hilic Negative  ## Hilic Negative  ## Hilic Negative  

hilic.neg_t<- filter(baltica, Column == "HILIC" & Mode == "neg" )
# extract compound names 
compound_names_tmp <-hilic.neg_t[,1]$Compound
# add the prefix COMPOUND infront of every compound name 
compound_names <- sprintf('COMPOUND_%s',  compound_names_tmp)
# extract compounds and transpose the matrix
hilic.neg_tmp <- hilic.neg_t[,10:407]
hilic.neg <- t.data.frame(hilic.neg_tmp)
hilic.neg <-data.frame(hilic.neg[2:398,])
# name the columns of the matrix with the compound names
names(hilic.neg) <- compound_names
# combine GROUP and RUN_ORDER with data
combine.hilic.neg <- cbind.data.frame(GROUP,RUN_ORDER,hilic.neg)


### Hilic Positive ### Hilic Positive ### Hilic Positive
hilic.pos_t <- filter(baltica, Column == "HILIC" & Mode == "pos" )

# extract compound names 
compound_names_tmp <-hilic.pos_t[,1]$Compound
# add the prefix COMPOUND infront of every compound name 
compound_names <- sprintf('COMPOUND_%s',  compound_names_tmp)
# extract compounds and transpose the matrix
hilic.pos_tmp <- hilic.pos_t[,10:407]
hilic.pos <- t.data.frame(hilic.pos_tmp)
hilic.pos <-data.frame(hilic.pos[2:398,])
# name the columns of the matrix with the compound names
names(hilic.pos) <- compound_names
# combine GROUP and RUN_ORDER with data
combine.hilic.pos <- cbind.data.frame(GROUP,RUN_ORDER,hilic.pos)

### RP pos ### RP pos  ### RP pos  
rp.pos_t<- filter(baltica, Column == "RP" & Mode == "pos" )
# extract compound names 
compound_names_tmp <-rp.pos_t[,1]$Compound
# add the prefix COMPOUND infront of every compound name 
compound_names <- sprintf('COMPOUND_%s',  compound_names_tmp)
# extract compounds and transpose the matrix
rp.pos_tmp <- rp.pos_t[,10:407]
rp.pos <- t.data.frame(rp.pos_tmp)
rp.pos <-data.frame(rp.pos[2:398,])
# name the columns of the matrix with the compound names
names(rp.pos) <- compound_names
# combine GROUP and RUN_ORDER with data
combine.rp.pos <- cbind.data.frame(GROUP,RUN_ORDER,rp.pos)

## RP Negative  ## RP Negative  ## RP Negative
rp.neg_t <- filter(baltica, Column == "RP" & Mode == "neg")
# extract compound names 
compound_names_tmp <-rp.neg_t[,1]$Compound
# add the prefix COMPOUND infront of every compound name 
compound_names <- sprintf('COMPOUND_%s',  compound_names_tmp)
# extract compounds and transpose the matrix
rp.neg_tmp <- rp.neg_t[,10:407]
rp.neg <- t.data.frame(rp.neg_tmp)
rp.neg <-data.frame(rp.neg[2:398,])
# name the columns of the matrix with the compound names
names(rp.neg) <- compound_names
# combine GROUP and RUN_ORDER with data
combine.rp.neg <- cbind.data.frame(GROUP,RUN_ORDER,rp.neg)


# save data
save(combine.hilic.neg,combine.hilic.pos,combine.rp.neg,combine.rp.pos,group.info,file = paste0(output.path, "data/","full_compound_data.Rdata"))
