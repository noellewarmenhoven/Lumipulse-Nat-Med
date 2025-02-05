#Bootstrapping thresholds for plasma biomarkers

library(boot)
library(cutpointr) 
library(pROC)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(readxl)

##Read datafile ----
lumi <- read_xlsx("BFMC_data.xlsx")
lumi$csf_status <- as.factor(as.numeric(lumi$csf_status))


##Thresholds for Lumi plasma p-tau217----

f_thresholds <- function(data, indices){
  d <- data[indices, ]
  roc_curve <- pROC::roc(csf_status~pl_ptau217, data=d, quiet=T)
  coordinates <- coords(roc_curve, x="all", transpose=F)
  
  #Specificity at 90%
  find_threshold <- coordinates[coordinates$specificity >= 0.895 & coordinates$specificity <= 0.904, ]
  threshold_spec90 <- if (nrow(find_threshold) > 0) {
    find_threshold[which.max(find_threshold$sensitivity), 1]
  } else {
    NA
  }
  
  # #Specificity at 95%
  find_threshold_spec95 <- coordinates[coordinates$specificity >= 0.946 & coordinates$specificity <= 0.954, ]
  threshold_spec95 <- if (nrow(find_threshold_spec95) > 0) {
    find_threshold_spec95[which.max(find_threshold_spec95$sensitivity), 1]
  } else {
    NA
  }
  
  #Sensitivity at 95
  find_threshold_sens95 <- coordinates[coordinates$sensitivity >= 0.94555 & coordinates$sensitivity <= 0.954, ]
  threshold_sens95 <- if (nrow(find_threshold_sens95) > 0) {
    find_threshold_sens95[which.max(find_threshold_sens95$specificity), 1]
  } else {
    NA
  }
  #Youden
  threshold_youden1 <- coords(roc_curve, "best", best.method = "youden", ret = "threshold")
  threshold_youden <- threshold_youden1$threshold
  return(c(threshold_spec90, threshold_spec95, threshold_sens95, threshold_youden))
}

set.seed(123)
boot_res <- boot(data = lumi, statistic = f_thresholds, R = 1000)

###mean of bootstraps and 95%CIs and median ----

#90% Specificity
lumi_90spec <- mean(boot_res$t[,1], na.rm=T)
# lumi_ci_spec90 <- boot.ci(boot_res, type = "perc",seed=12345,index = 1)
# lumi_90spec_ci_low <- round(lumi_ci_spec90$percent[4],3)
# lumi_90spec_ci_upper <-round(lumi_ci_spec90$percent[5],3)
# lumi_ci_spec90_formatted <- paste0(round(lumi_90spec,3), " (", lumi_90spec_ci_low, "-",lumi_90spec_ci_upper, ")")

#95% Specificity
lumi_95spec <- mean(boot_res$t[,2], na.rm=T)
# lumi_ci_spec95 <- boot.ci(boot_res, type = "perc",seed=12345,index = 2)
# lumi_95spec_ci_low <- round(lumi_ci_spec95$percent[4],3)
# lumi_95spec_ci_upper <-round(lumi_ci_spec95$percent[5],3)
# lumi_ci_spec95_formatted <- paste0(round(lumi_95spec,3), " (", lumi_95spec_ci_low, "-",lumi_95spec_ci_upper, ")")

#95% Sensitivity
lumi_95sens <-mean(boot_res$t[,3], na.rm=T)
# lumi_ci_sens95 <- boot.ci(boot_res, type = "perc",seed=12345,index = 3)
# lumi_95sens_ci_low <- round(lumi_ci_sens95$percent[4],3)
# lumi_95sens_ci_upper <-round(lumi_ci_sens95$percent[5],3)
# lumi_ci_sens95_formatted <- paste0(round(lumi_95sens,3), " (", lumi_95sens_ci_low, "-",lumi_95sens_ci_upper, ")")

#Youden
lumi_youden <-mean(boot_res$t[,4], na.rm=T)
# lumi_ci_youden <- boot.ci(boot_res, type = "perc",seed=12345,index = 4)
# lumi_youden_ci_low <- round(lumi_ci_youden$percent[4],3)
# lumi_youden_ci_upper <-round(lumi_ci_youden$percent[5],3)
# lumi_ci_youden_formatted <- paste0(round(lumi_youden,3), " (", lumi_youden_ci_low, "-",lumi_youden_ci_upper, ")")


##Thresholds for C2N plasma p-tau217----

library(pROC)
f_thresholds_c2n <- function(data, indices){
  d <- data[indices, ]
  roc_curve <- pROC::roc(csf_status~p_tau217_v2, data=d, quiet=TRUE)
  coordinates <- coords(roc_curve, x="all", transpose=F)
  
  #Specificity at 90%
  find_threshold <- coordinates[coordinates$specificity >= 0.895 & coordinates$specificity <= 0.904, ]
  threshold_spec90 <- if (nrow(find_threshold) > 0) {
    find_threshold[which.max(find_threshold$sensitivity), 1]
  } else {
    NA
  }
  
  # #Specificity at 95%
  find_threshold_spec95 <- coordinates[coordinates$specificity >= 0.946 & coordinates$specificity <= 0.954, ]
  threshold_spec95 <- if (nrow(find_threshold_spec95) > 0) {
    find_threshold_spec95[which.max(find_threshold_spec95$sensitivity), 1]
  } else {
    NA
  }
  
  #Sensitivity at 95
  find_threshold_sens95 <- coordinates[coordinates$sensitivity >= 0.94555 & coordinates$sensitivity <= 0.954, ]
  threshold_sens95 <- if (nrow(find_threshold_sens95) > 0) {
    find_threshold_sens95[which.max(find_threshold_sens95$specificity), 1]
  } else {
    NA
  }

  #Youden
  threshold_youden1 <- coords(roc_curve, "best", best.method = "youden", ret = "threshold")
  threshold_youden <- threshold_youden1$threshold
  return(c(threshold_spec90, threshold_spec95, threshold_sens95, threshold_youden))
}

set.seed(123)
boot_res_c2n <- boot(data = lumi, statistic = f_thresholds_c2n, R = 2000)

###mean of bootstraps and 95%CIs and median ----

#90% Specificity
c2n_90spec <-mean(boot_res_c2n$t[,1], na.rm=T)
# c2n_ci_spec90 <- boot.ci(boot_res_c2n, type = "perc",seed=12345,index = 1)
# c2n_90spec_low <- round(c2n_ci_spec90$percent[4],3)
# c2n_90spec_high <-round(c2n_ci_spec90$percent[5],3)
# c2n_90spec_formatted <- paste0(round(c2n_90spec,3), " (", c2n_90spec_low, "-",c2n_90spec_high, ")")

#95% Specificity
c2n_95spec <-mean(boot_res_c2n$t[,2], na.rm=T)
# c2n_ci_spec95 <- boot.ci(boot_res_c2n, type = "perc",seed=12345,index = 2)
# c2n_95spec_low <- round(c2n_ci_spec95$percent[4],3)
# c2n_95spec_high <-round(c2n_ci_spec95$percent[5],3)
# c2n_95spec_formatted <- paste0(round(c2n_95spec,3), " (", c2n_95spec_low, "-",c2n_95spec_high, ")")

#95% Sensitivity
c2n_95sens <-mean(boot_res_c2n$t[,3], na.rm=T)
# c2n_ci_sens95 <- boot.ci(boot_res_c2n, type = "perc",seed=12345,index = 3)
# c2n_95sens_low <- round(c2n_ci_sens95$percent[4],3)
# c2n_95sens_high <-round(c2n_ci_sens95$percent[5],3)
# c2n_95sens_formatted <- paste0(round(c2n_95sens,3), " (", c2n_95sens_low, "-",c2n_95sens_high, ")")

#Youden
c2n_youden <-mean(boot_res_c2n$t[,4], na.rm=T)
# c2n_ci_youden <- boot.ci(boot_res_c2n, type = "perc",seed=12345,index = 4)
# c2n_youden_low <- round(c2n_ci_youden$percent[4],3)
# c2n_youden_high <-round(c2n_ci_youden$percent[5],3)
# c2n_youden_formatted <- paste0(round(c2n_youden,3), " (", c2n_youden_low, "-",c2n_youden_high, ")")


##Thresholds for C2N plasma %p-tau217----
lumi  <- lumi %>% drop_na(p_tau217_ratio_v2, csf_status)
f_thresholds_c2nratio <- function(data, indices){
  d <- data[indices, ]
  roc_curve <- pROC::roc(csf_status~p_tau217_ratio_v2, data=d, quiet=T)
  coordinates <- coords(roc_curve, x="all", transpose=F)

    #Specificity at 90%
  find_threshold <- coordinates[coordinates$specificity >= 0.895 & coordinates$specificity <= 0.904, ]
  threshold_spec90 <- if (nrow(find_threshold) > 0) {
    find_threshold[which.max(find_threshold$sensitivity), 1]
  } else {
    NA
  }
  
   #Specificity at 95%
  find_threshold_spec95 <- coordinates[coordinates$specificity >= 0.946 & coordinates$specificity <= 0.954, ]
  threshold_spec95 <- if (nrow(find_threshold_spec95) > 0) {
    find_threshold_spec95[which.max(find_threshold_spec95$sensitivity), 1]
  } else {
    NA
  }
  
  #Sensitivity at 95
  find_threshold_sens95 <- coordinates[coordinates$sensitivity >= 0.94555 & coordinates$sensitivity <= 0.954, ]
  threshold_sens95 <- if (nrow(find_threshold_sens95) > 0) {
    find_threshold_sens95[which.max(find_threshold_sens95$specificity), 1]
  } else {
    NA
  }
  #Youden
  threshold_youden1 <- coords(roc_curve, "best", best.method = "youden", ret = "threshold")
  threshold_youden <- threshold_youden1$threshold
  return(c(threshold_spec90, threshold_spec95, threshold_sens95, threshold_youden))
}

set.seed(123)
boot_res_c2nratio <- boot(data = lumi, statistic = f_thresholds_c2nratio, R = 2000)

###mean of bootstraps and 95%CIs and median ----

#90% Specificity
c2nratio_90spec <- mean(boot_res_c2nratio$t[,1], na.rm=T)
# c2nratio_ci_spec90 <- boot.ci(boot_res_c2nratio, type = "perc",seed=12345,index = 1)
# c2n_ratio_spec90_low <- round(c2nratio_ci_spec90$percent[4],3)
# c2n_ratio_spec90_high <-round(c2nratio_ci_spec90$percent[5],3)
# c2n_ratio_spec90_formatted <- paste0(round(c2nratio_90spec,3), " (", c2n_ratio_spec90_low, "-",c2n_ratio_spec90_high, ")")

#95% Specificity
c2nratio_95spec <-mean(boot_res_c2nratio$t[,2], na.rm=T)
# c2nratio_ci_spec95 <- boot.ci(boot_res_c2nratio, type = "perc",seed=12345,index = 2)
# c2n_ratio_spec95_low <- round(c2nratio_ci_spec95$percent[4],3)
# c2n_ratio_spec95_high <-round(c2nratio_ci_spec95$percent[5],3)
# c2n_ratio_spec95_formatted <- paste0(round(c2nratio_95spec,3), " (", c2n_ratio_spec95_low, "-",c2n_ratio_spec95_high, ")")

#95% Sensitivity
c2nratio_95sen <- mean(boot_res_c2nratio$t[,3], na.rm=T)
# c2nratio_ci_sens95 <- boot.ci(boot_res_c2nratio, type = "perc",seed=12345,index = 3)
# c2n_ratio_sens95_low <- round(c2nratio_ci_sens95$percent[4],3)
# c2n_ratio_sens95_high <-round(c2nratio_ci_sens95$percent[5],3)
# c2n_ratio_sens95_formatted <- paste0(round(c2nratio_95sen,3), " (", c2n_ratio_sens95_low, "-",c2n_ratio_sens95_high, ")")

#Youden
c2nratio_youden <- mean(boot_res_c2nratio$t[,4], na.rm=T)
# c2nratio_ci_youden <- boot.ci(boot_res_c2nratio, type = "perc",seed=12345,index = 4)
# c2n_ratio_youden_low <- round(c2nratio_ci_youden$percent[4],3)
# c2n_ratio_youden_high <-round(c2nratio_ci_youden$percent[5],3)
# c2n_ratio_youden_formatted <- paste0(round(c2nratio_youden,3), " (", c2n_ratio_youden_low, "-",c2n_ratio_youden_high, ")")


#Thresholds for p-tau217ab42 Lumipulse----
#lumi <- read_xlsx("BFMC_data.xlsx") 

lumi$ptau217_ab42 <- (lumi$pl_ptau217)/(lumi$ab42)
lumi <- lumi %>% drop_na(ptau217_ab42)

f_thresholds <- function(data, indices){
  d <- data[indices, ]
  roc_curve <- pROC::roc(csf_status~ptau217_ab42, data=d, quiet=T)
  coordinates <- coords(roc_curve, x="all", transpose=F)
  
  #Specificity at 90%
  find_threshold <- coordinates[coordinates$specificity >= 0.895 & coordinates$specificity <= 0.904, ]
  threshold_spec90 <- if (nrow(find_threshold) > 0) {
    find_threshold[which.max(find_threshold$sensitivity), 1]
  } else {
    NA
  }
  
  # #Specificity at 95%
  find_threshold_spec95 <- coordinates[coordinates$specificity >= 0.946 & coordinates$specificity <= 0.954, ]
  threshold_spec95 <- if (nrow(find_threshold_spec95) > 0) {
    find_threshold_spec95[which.max(find_threshold_spec95$sensitivity), 1]
  } else {
    NA
  }
  
  #Sensitivity at 95
  find_threshold_sens95 <- coordinates[coordinates$sensitivity >= 0.94555 & coordinates$sensitivity <= 0.954, ]
  threshold_sens95 <- if (nrow(find_threshold_sens95) > 0) {
    find_threshold_sens95[which.max(find_threshold_sens95$specificity), 1]
  } else {
    NA
  }
  #Youden
  threshold_youden1 <- coords(roc_curve, "best", best.method = "youden", ret = "threshold")
  threshold_youden <- threshold_youden1$threshold
  return(c(threshold_spec90, threshold_spec95, threshold_sens95, threshold_youden))
}

set.seed(123)
boot_res <- boot(data = lumi, statistic = f_thresholds, R = 2000)

###mean of bootstraps and 95%CIs and median ----

#90% Specificity
lumi_90spec <- mean(boot_res$t[,1], na.rm=T)
# lumi_ci_spec90 <- boot.ci(boot_res, type = "perc",seed=12345,index = 1)
# lumi_90spec_ci_low <- round(lumi_ci_spec90$percent[4],3)
# lumi_90spec_ci_upper <-round(lumi_ci_spec90$percent[5],3)
# lumi_ci_spec90_formatted <- paste0(round(lumi_90spec,3), " (", lumi_90spec_ci_low, "-",lumi_90spec_ci_upper, ")")

#95% Specificity
lumi_95spec <- mean(boot_res$t[,2], na.rm=T)
# lumi_ci_spec95 <- boot.ci(boot_res, type = "perc",seed=12345,index = 2)
# lumi_95spec_ci_low <- round(lumi_ci_spec95$percent[4],3)
# lumi_95spec_ci_upper <-round(lumi_ci_spec95$percent[5],3)
# lumi_ci_spec95_formatted <- paste0(round(lumi_95spec,3), " (", lumi_95spec_ci_low, "-",lumi_95spec_ci_upper, ")")

#95% Sensitivity
lumi_95sens <-mean(boot_res$t[,3], na.rm=T)
# lumi_ci_sens95 <- boot.ci(boot_res, type = "perc",seed=12345,index = 3)
# lumi_95sens_ci_low <- round(lumi_ci_sens95$percent[4],3)
# lumi_95sens_ci_upper <-round(lumi_ci_sens95$percent[5],3)
# lumi_ci_sens95_formatted <- paste0(round(lumi_95sens,3), " (", lumi_95sens_ci_low, "-",lumi_95sens_ci_upper, ")")

