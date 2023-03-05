
# std libs
library(dplyr)
library(mice)

base::source("scripts/data_preprocessing.R")


# use mice for imputation -------------------------------------------------

sort(apply(dat.observe,2,pMiss), decreasing=TRUE )



