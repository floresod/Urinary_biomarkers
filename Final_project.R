##### R script for the data analysis and models for final project #####

#### Libraries ####
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

#adjut R options 
options(digits = 4, scipen = 999)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#### data ####
#load data 
raw_data <- read.csv(file = "raw_data.csv", check.names = FALSE)
#correct col1: sample_id
colnames(raw_data)[1] <- "sample_id"

#adjust levels of diagnosis levels 
raw_data$diagnosis <- factor(raw_data$diagnosis) %>% 
  fct_recode("healthy" = "1", "benign" = "2", "cancer" = "3" )

#tidy format - put biomarkers in a column
tidy_data <- raw_data %>%
  pivot_longer(cols = c(creatinine, LYVE1, REG1B, TFF1, REG1A, plasma_CA19_9), 
               names_to = "biomarker", values_to = "levels" )
#sort biomarkers 
tidy_data$biomarker <- factor(tidy_data$biomarker, 
                              levels = c("creatinine", "LYVE1", "REG1B", "TFF1", "plasma_CA19_9", "REG1A"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#### data exploration ####
#data features 
summary(raw_data)

# biomarkers in each diagnosis group 
tidy_data %>%
  ggplot(aes(x = biomarker, y = levels, color = diagnosis )) + 
  geom_boxplot() + 
  scale_y_continuous(trans = "log10") + 
  theme_bw()



