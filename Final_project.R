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
               names_to = "biomarker", values_to = "levels" ) %>% 
  filter(!is.na(levels))
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
  facet_grid(~sex) + 
  theme_bw()
#levels of some biomarkers seem to be different in groups with cancer

#distribution of diagnosis by sex 
#by sex
table(raw_data$sex, raw_data$diagnosis) %>% prop.table(1)
#by diagnosis 
table(raw_data$sex, raw_data$diagnosis) %>% prop.table(2)

# it seems that male are more prone to get cancer 

## distribution of diagnosis by age 
tidy_data %>% group_by(age, diagnosis, biomarker) %>%
  summarise(avg_lvs = mean(levels)) %>% 
  ggplot(aes(x = as.numeric(age), y = avg_lvs, color = diagnosis)) + 
  geom_line() +
  labs(title = "Figure Biomarkers concentration in diagnosis group vs age", 
       x = "age", y = "concentration") +
  scale_x_continuous(breaks = seq(25,90,5)) + 
  scale_y_continuous(trans = "log10") +
  facet_wrap(~biomarker, ncol = 1) +
  theme_bw()
#some biomarkers seem to be more impotant over time 

## proportion of diagnosis by cohort 
table(raw_data$sample_origin, raw_data$diagnosis)


