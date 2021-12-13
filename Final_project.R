##### R script for the data analysis and models for final project #####

#Install libraries
if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(ggplot2)) install.packages("ggplot2", repos = "http://cran.us.r-project.org")
if(!require(ggpubr)) install.packages("ggpubr", repos = "http://cran.us.r-project.org")
#if(!require(groom)) install.packages("ggpubr", repos = "http://cran.us.r-project.org")

#load libraries 
library(tidyverse)
library(caret)
library(ggplot2)

#load data 
raw_data <- read.csv(file = "raw_data.csv", check.names = FALSE)

