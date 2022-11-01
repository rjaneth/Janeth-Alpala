rm(list = ls(all = TRUE))
graphics.off()
shell("cls")

getwd()



library(dplyr) 
library(pROC) 
library(MASS)
library(broom) 
library(AUC)
library(pROC)
library(tidyverse)
dplyr::roc()