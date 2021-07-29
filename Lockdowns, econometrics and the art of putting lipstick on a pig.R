library(lfe)
library(stargazer)
library(knitr)
library(plm)
library(latex2exp)
library(ggplot2)
library(ggthemes)
library(estimatr)
library(gridExtra)
library(grid)
library(mvtnorm)
library(kableExtra)
library(data.table)
library(reshape2)
library(AER)
library(hdm)
library(randomForest)
library(glmnet)
library(magrittr)
library(readr)
library(stringr)
library(latexpdf)
library(cowplot)
library(rlist)

# allow parallel processing
#future::plan(future::multisession)

# this variable is needed by some of the files that are run below, so make sure that you have set
# your working directory to the root of the directory
rootdir <- paste(getwd(), "Original repository with modifications", sep = "/")

# this variable is used to choose which version of the dataset to use (alternative datasets with up to date
# epidemic data and/or data on policies were created on 2021/06/24)
dataset_version <- "Original dataset"
# dataset_version <- "Up to date epidemic data but data on policies as in paper"
# dataset_version <- "Up to date epidemic data and data on policies as in latest CUSP database"

# set to TRUE in order to use data on cases and deaths from JHU instead of the NYT as in the paper
JHU <- FALSE

# set to FALSE if you don't want to smooth the policy dummy variables
policy_smoothing <- TRUE

# those variables are needed for some of the code in regprep.R that is needed for the double machine learning
# sensitivity analysis, which I don't use but that I keep because I don't want to make unnecessary changes to
# the original files to limit the risk that I might inadvertently break something
L.c <- 14
L.d <- 21

source(paste(rootdir,"cases_and_policies/R/utils.R",sep="/"))
source(paste(rootdir,"cases_and_policies/R/dataprep - updated.R", sep="/"), local = TRUE)
source(paste(rootdir,"cases_and_policies/rmd/varlabels.R", sep="/"))
source(paste(rootdir,"cases_and_policies/rmd/regprep - updated.R", sep="/"), local = TRUE)
# the only changes I made in this file are to avoid running unnecessary regressions and speed things up
source(paste(rootdir,"cases_and_policies/rmd/generatetables - updated.R", sep="/"), local = TRUE)
source(paste(rootdir,"cases_and_policies/rmd/bootstrap_felm.R", sep="/"), local = TRUE)
source(paste(rootdir,"cases_and_policies/rmd/smoothtests.R", sep="/"), local = TRUE)

dataset_version <- ifelse(!policy_smoothing, paste0(dataset_version, " - No policy smoothing"), dataset_version)
dataset_version <- ifelse(JHU, paste0(dataset_version, " - JHU"), dataset_version)

# this doesn't work the first time on my computer, which is probably due to a
# bug with furrr (see comments in file), but it does after I try again
source("Reanalysis of the results/Reanalysis of the results.R", local = TRUE)

# I must run "Model on simulated data.R" after I have run all the regressions I want on the actual data
# because the code for the simulation loads dplyr and there is a weird conflict with lfe that improperly
# changes the results of the regressions if dplyr has been loaded before I run the regressions
source("Model on simulated data/Model on simulated data.R", local = TRUE)

source("Comparison of the timing of policies in different versions of the dataset/Comparison of the timing of policies in different versions of the dataset.R", local = TRUE)
