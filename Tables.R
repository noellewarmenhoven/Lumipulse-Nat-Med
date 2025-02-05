#This file contains the codes used to create Table 1 & Supplementary tables

library(dplyr)
library(pROC)
library(boot)
library(cutpointr)
library(readxl)
library(table1)
library(gtsummary)
library(tidyverse)

#Table 1----
#list  datasets

BFMC <- read_xlsx("BFMC_data.xlsx")

Goth <- read_xlsx("Gothenburg_data.xlsx")

BBRC <- read_xlsx("BBRC_data.xlsx")

Bres <- read_xlsx("Brescia_data.xlsx")

BFPC <- read_xlsx("BFPC_data.xlsx")

lumi <- rbind(BFMC, Goth, BBRC, Bres, BFPC)

lumi$ptau217_spec90_lumi <- ifelse(lumi$pl_ptau217 > 0.27, 1, 0) 

foldchange <- median(lumi$pl_ptau217[lumi$ptau217_spec90_lumi == 1]) / median(lumi$pl_ptau217[lumi$ptau217_spec90_lumi == 0])

#alter variable types 
lumi$mmse_score <- as.numeric(as.character(lumi$mmse_score))
lumi$sex <- as.factor(as.numeric(lumi$sex))
lumi$cognitive_status <- as.factor(as.numeric(lumi$cognitive_status))
lumi$apoe <- as.factor(as.numeric(lumi$apoe))
lumi$CKD <- as.factor(as.numeric(lumi$CKD))
lumi$diabetes <- as.factor(as.numeric(lumi$diabetes))
lumi$csf_status <- as.factor(as.numeric(lumi$csf_status))
lumi$ab42_40 <- as.numeric(as.character(lumi$ab42_40))
lumi$ptau181 <- as.numeric(as.character(lumi$ptau181))
#add cohort variable for stratification
lumi$Cohort <- factor(lumi$Cohort, levels=c("Malmo", "Gothenburg", "Barcelona", "Brescia", "Sweden"))

#Create table 1
library(table1)
table1(~ sex + age + education + mmse_score + cognitive_status + apoe + CKD + diabetes +  csf_status + ab42_40 + ptau181 + csf_ab42_ptau181 + pl_ptau217 + p_tau217_v2 + p_tau217_ratio_v2 | Cohort, data=lumi,
       overall=c(left="Total"))


#eTable 2----
#list  datasets

BFMC <- read_xlsx("BFMC_data_LP.xlsx")

BBRC <- read_xlsx("BBRC_data_LP.xlsx")

Bres <- read_xlsx("Brescia_data_LP.xlsx")

Goth <- read_xlsx("Gothenburg_data_LP.xlsx")

BFPC <- read_xlsx("BFPC_data_LP.xlsx")

#create pooled group
pooled <- rbind(BFMC, BBRC, Bres, Goth)
pooled$Cohort2 <- "Pooled"
BFPC$Cohort2 <- "BFPC"

#bind data together
lumi <- rbind(pooled, BFPC)

#create p-tau217/AB42 ratio
lumi$ptau217_ab42 <- (lumi$pl_ptau217)/(lumi$ab42)
lumi$ptau217_ab42[is.infinite(lumi$ptau217_ab42)] <- NA 
lumi <- lumi %>% drop_na(ptau217_ab42)

#change variable types
lumi$ab42 <- as.numeric(lumi$ab42)
lumi$mmse_score <- as.numeric(as.character(lumi$mmse_score))
lumi$sex <- as.factor(as.numeric(lumi$sex))
lumi$cognitive_status <- as.factor(as.numeric(lumi$cognitive_status))
lumi$apoe <- as.factor(as.numeric(lumi$apoe))
lumi$CKD <- as.factor(as.numeric(lumi$CKD))
lumi$diabetes <- as.factor(as.numeric(lumi$diabetes))
lumi$csf_status <- as.factor(as.numeric(lumi$csf_status))
lumi$ab42_40 <- as.numeric(as.character(lumi$ab42_40))
lumi$ptau181 <- as.numeric(as.character(lumi$ptau181))
lumi$Cohort2 <- factor(lumi$Cohort2, levels=c("Pooled","BFPC"))
lumi$education <- as.numeric(lumi$education)
  
#Create table
library(table1)
table1(~ sex + age + education + mmse_score + cognitive_status + apoe + CKD + diabetes +  csf_status + ab42_40 + ptau181 + csf_ab42_ptau181 + pl_ptau217 + ptau217_ab42 + Cohort | Cohort2, data=lumi,
       overall=c(left="Total"))


#eTable 3----
#list datasets

BFMC <- read_xlsx("BFMC_data.xlsx")
BFPC <- read_xlsx("BFPC_data.xlsx")
Got <- read_xlsx("Gothenburg_data_C2N.xlsx")
Bres <- read_xlsx("Brescia_data_C2N.xlsx")

#create pooled group
pooled <- rbind(BFMC, Got, Bres)
pooled$Cohort2 <- "Pooled"
BFPC$Cohort2 <- "BFPC"

#bind data together
lumi <- rbind(pooled, BFPC)

lumi <- lumi %>% drop_na(p_tau217_ratio_v2)

#change variable types
lumi$mmse_score <- as.numeric(as.character(lumi$mmse_score))
lumi$sex <- as.factor(as.numeric(lumi$sex))
lumi$cognitive_status <- as.factor(as.numeric(lumi$cognitive_status))
lumi$apoe <- as.factor(as.numeric(lumi$apoe))
lumi$CKD <- as.factor(as.numeric(lumi$CKD))
lumi$diabetes <- as.factor(as.numeric(lumi$diabetes))
lumi$csf_status <- as.factor(as.numeric(lumi$csf_status))
lumi$ab42_40 <- as.numeric(as.character(lumi$ab42_40))
lumi$ptau181 <- as.numeric(as.character(lumi$ptau181))
lumi$Cohort2 <- factor(lumi$Cohort2, levels=c("Pooled","BFPC"))

#create table
table1(~ sex + age + education + mmse_score + cognitive_status + apoe + CKD + diabetes +  csf_status + ab42_40 + ptau181 + csf_ab42_ptau181 + pl_ptau217 + p_tau217_v2 + p_tau217_ratio_v2 | Cohort2, data=lumi,
       overall=c(left="Total"))
