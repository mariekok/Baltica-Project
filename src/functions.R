
#' Title baltica_prefiltering
#'
#' @param data 
#' @param group.info 
#'
#' @return
#' @export
#'
#' @examples
baltica_prefiltering <- function(data ,GROUP){
  
  require(MUVR)
  # Perform anova in the X as a prefiltering procedure
  
  # Here I will use nonparametric anova because the metabolites are a bit skewed.
  
  # Kruskal–Wallis Test
  
  # http://rcompanion.org/handbook/F_08.html
  # https://rcompanion.org/rcompanion/d_06.html
  pvalues_BALTICA <- NULL
  # remove the columns of GROUP and RUN_ORDER
  
  tmp.data <- data
  
  # some statistics
  
  # multiple testing with wilcoxon test and ttest 
  for(i in 1:dim(tmp.data)[2]) {
    result_ANOVA <- kruskal.test(tmp.data[,i]~ GROUP)
    pvalues_BALTICA[i] <-result_ANOVA$p.value
    
  }
  
  # FDR correction for control versus helathy from wilcoxon test 
  
  adj.pvalues_BALTICA <- p.adjust(pvalues_BALTICA, method = "fdr", n = length(pvalues_BALTICA))
  
  responses <- colnames(tmp.data)
  #View(responses)
  # combine p values and and adjusted pvalues into one dataframe
  BALTICA_all_pvalues <- cbind.data.frame(responses,pvalues_BALTICA,adj.pvalues_BALTICA)

  # select the ones that are below the significant level of 0.9 in the FDR (this number is a bit arbitary)
  index <- which(BALTICA_all_pvalues$adj.pvalues_BALTICA <= 0.9)
  smallpvalues <-BALTICA_all_pvalues[index,]

  # data with the most significant metabolites from pvalues
  data.prefiltering <- subset(tmp.data, select = as.character(smallpvalues$responses))
  
  data.prefiltering
  
}


#' Title baltica_glog
#' Apply glog tranformation in the data
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
baltica_glog <-function(data){ 
  
  #
  # apply glog 
  require(FitAR)
  log_baltica <-glog(data, a = 1, InverseQ = FALSE)
  
  log_baltica
  
}





#' Title baltica_pqnNormalization
#' Apply PQN normalization in the data based on the QCs
#' @param data data
#' @param QC   only QC samples 

#'
#' @return
#' @export
#'
#' @examples
baltica_pqnNormalization <- function(data , QC ){
  # # apply normalization
  require(KODAMA)

  # Normalization 
  # I will perform PQN normalization in order to do that I will need a reference sample and for that 
  # I will use a QC sample the log.average.QC which the mean of all QCs in the Baltica data
  # Extract only qcs 
  # X should be only QC samples 
  # calculate the average of the QCs in the data

  average.QC <- as.data.frame(colMeans(QC))
  
  
  
  data_norm_temp <- KODAMA::normalization(data, method = "pqn",ref= average.QC)
  data_norm <- data_norm_temp$newXtrain
  
  data_norm
}



#' Title baltica_preprocessing
#' Mean cenetring and pareto scaling the data
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
baltica_preprocessing <- function(data){
  

  # mean center
  mu <-apply(data,2,mean)
  
  # pareto scaled the data 
  require(MetabolAnalyze)
  pareto_scale <- scaling(data, type = "pareto")
  
  return(list(meancenter = mu,paretoscaled = pareto_scale))
}


