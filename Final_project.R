##### R script for the data analysis and models for final project #####

#### Libraries ####
#Install libraries
if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(ggplot2)) install.packages("ggplot2", repos = "http://cran.us.r-project.org")
if(!require(scales)) install.packages("scales", repos = "http://cran.us.r-project.org")
if(!require(ggpubr)) install.packages("ggpubr", repos = "http://cran.us.r-project.org")

#load libraries 
library(tidyverse)
library(caret)
library(ggplot2)
library(scales)
library(ggpubr)

#adjust R options 
options(digits = 4, scipen = 999)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#### data ####

## data was downloaded from https://www.kaggle.com/johnjdavisiv/urinary-biomarkers-for-pancreatic-cancer ## 
# file is saved as 'raw_data.csv' 

#load data 
raw_data <- read.csv(file = "data/Debernardi et al 2020 data.csv", check.names = FALSE)

#correct col1 name: sample_id
colnames(raw_data)[1] <- "sample_id"

## defined level of different factor ## 
# patient_cohort 
raw_data$patient_cohort <- factor(raw_data$patient_cohort)
# sample origin 
raw_data$sample_origin <-  factor(raw_data$sample_origin)
# sex 
raw_data$sex <- factor(raw_data$sex)
# diagnosis
raw_data$diagnosis <- factor(raw_data$diagnosis) %>% 
  fct_recode("healthy" = "1", "benign" = "2", "cancer" = "3" )
# stage 
raw_data$stage <- factor(raw_data$stage)
# benign diagnosis 
raw_data$benign_sample_diagnosis <- factor(raw_data$benign_sample_diagnosis)

## ## ## ## ## 

#save as .rda file for PDF report 
save(raw_data, file = "rda/raw_data.rda")

## load documentation file with column's details ##
column_details <- read.csv("data/Debernardi et al 2020 documentation.csv")[,-2] #remove 2nd column with not necessary information
#change columns names 
names(column_details) <- c("Column name", "Details")
#save for PDF report 
save(column_details, file = "rda/columns_details.rda")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#### data exploration ####

### Features of the raw data ### 
summary(raw_data)

# proportion of healthy and cancer 
table(raw_data$diagnosis) %>% prop.table()

### Shape data for more analysis  ###
#tidy format - put biomarkers in a column
tidy_data <- raw_data %>% dplyr::select(-stage, -patient_cohort, -sample_origin, -stage, -benign_sample_diagnosis) %>% 
  pivot_longer(cols = c(creatinine, LYVE1, REG1B, TFF1, REG1A, plasma_CA19_9), 
               names_to = "biomarker", values_to = "levels" ) %>% 
  filter(!is.na(levels))

#sort biomarkers 
tidy_data$biomarker <- factor(tidy_data$biomarker, 
                              levels = c("creatinine", "LYVE1", "REG1B", "TFF1", "plasma_CA19_9", "REG1A"))
## ## ## 

### distribution of diagnosis by sex ### 
#by sex
table(raw_data$sex, raw_data$diagnosis) %>% prop.table(1)
#by diagnosis 
table(raw_data$sex, raw_data$diagnosis) %>% prop.table(2)

# statistical analysis #
sex_diag <- table(raw_data$sex, raw_data$diagnosis)
chisq.test(sex_diag)
# it seems that male are more prone to get cancer 

### distribution of diagnosis by age ### 
my_colors <- c("#56B4E9", "#0072B2", "#D55E00")

#Figure 1 Diagnosis by age 
fig.1 <- table(raw_data$age, raw_data$diagnosis) %>% data.frame() %>%
  mutate(age = as.character(Var1) %>% as.numeric()) %>%
  ggplot(aes(x = age, y = Freq, fill = Var2)) +  
  geom_bar(stat = "identity", position=position_dodge()) +
  labs(title = "Figure 1. Diagnosis by age", 
       x = "age (years)", 
       y = "number of cases") +
  scale_y_continuous(breaks = seq(0,12,2)) +
  scale_x_continuous(breaks = seq(25, 90, 5)) + 
  scale_fill_manual("diagnosis", values = my_colors ) + 
  theme_bw() + 
  theme(legend.position = "bottom")
#print
fig.1

#save for PDF report 
save(fig.1, file = "figs/fig.1")

# proportion of cancer cases > 50 years 
raw_data %>% filter(diagnosis == "cancer") %>% summarise(prop.50 = mean(age >= 50))


### Distribution of biomarkers' concentrations by diagnosis group ###
# figure 3 biomarkers by diagnosis nad sex 
fig.2 <- tidy_data %>%
  mutate(sex = ifelse(sex == "F", "Female", "Male")) %>% 
  ggplot(aes(x = biomarker, y = levels, color = diagnosis )) + 
  geom_boxplot() +
  labs(title = "Figure 2. Biomarkers levels by diagnosis group and sex", 
       x = "biomarker", 
       y = "concentration") + 
  scale_y_continuous(trans = "log10", breaks = c(0, 0.01, 10, 10000) ) + 
  scale_x_discrete(guide = guide_axis(n.dodge= 2)) +
  scale_color_manual(values = my_colors) + 
  facet_grid(~sex) + 
  theme_bw() + 
  theme(legend.position = "bottom")

#print fig 2
fig.2

#save for PDF report 
save(fig.2, file = "figs/fig.2")


## distribution of diagnosis by age ##
fig.3 <- tidy_data %>% 
  group_by(age, diagnosis, biomarker) %>%
  summarise(avg_lvs = mean(levels)) %>% 
  ggplot(aes(x = as.numeric(age), y = avg_lvs, color = diagnosis)) + 
  geom_line() +
  labs(title = "Figure 3. Average biomarker concentration by age and diagnosis group", 
       x = "age", y = "concentration") +
  scale_x_continuous(breaks = seq(25,90,5)) + 
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = my_colors) +
  facet_wrap(~biomarker, ncol = 1) +
  theme_bw()

#print figure 3 
fig.3
#save fro PDF report 
save(fig.3, file = "figs/fig.3")


###### Multidimensional scaling: PCoA to evaluate similarities #####
#install package 
if(!require(vegan)) install.packages("vegan", repos = "http://cran.us.r-project.org")
library(vegan)


#Using protein panel (biomarkers) only 
#create matrix: biomarkers only 
matrix_pcoa <- raw_data %>% 
  dplyr::select(creatinine, LYVE1, REG1B, TFF1, REG1A, plasma_CA19_9) %>%
  `rownames<-`(raw_data$sample_id) %>%
  as.matrix()

#calculate distance 
dist.res <-  vegdist(matrix_pcoa, method = "bray", na.rm = TRUE)
mds.res = cmdscale(dist.res, eig = TRUE, x.ret = TRUE)
mds_varperc = round(mds.res$eig/sum(mds.res$eig)*100, 1)

##### Convert data into dataframe for ggplot #####
#dataframe with sample's information
sample_info <- raw_data %>% 
  dplyr::select(sample_id, sex, diagnosis) %>% 
  mutate(y = ifelse(diagnosis == "cancer", "cancer", "benign"))

#table with PCoA and sample's information
mds_data = data.frame(sample_id = row.names(mds.res$points), 
                             X=mds.res$points[,1], 
                             Y=mds.res$points[,2]) %>% 
  left_join(sample_info) # include sample's details 

## plot PCoA biomarkers only## 
fig.4a <- mds_data %>% 
  ggplot(aes(x = X, y = Y, color = y, shape = sex)) + 
  geom_point(alpha = 0.8) +
  labs(title = "Biomarkers only", 
       x = paste("PC1 - ", mds_varperc[1], "%", sep = ""), 
       y = paste("PC2 -", mds_varperc[2], "%", sep = "")) + 
  scale_color_manual("diagnosis", values = c("#56B4E9","#D55E00")) + 
  guides(shape = guide_legend(override.aes = list(1)), 
         color = guide_legend(override.aes = list(1))) +
  theme_bw() +
  theme(legend.title = element_text(size = 9), 
        legend.text = element_text(size = )) 
# print
fig.4a
#save for PDF report 
save(fig.4a, file = "figs/fig.4a")

#Using protein panel + age
#create matrix: biomarkers only 
matrix_pcoa <- raw_data %>% 
  dplyr::select(age, creatinine, LYVE1, REG1B, TFF1, REG1A, plasma_CA19_9) %>%
  `rownames<-`(raw_data$sample_id) %>%
  as.matrix()

#calculate distance 
dist.res <-  vegdist(matrix_pcoa, method = "bray", na.rm = TRUE)
mds.res = cmdscale(dist.res, eig = TRUE, x.ret = TRUE)
mds_varperc = round(mds.res$eig/sum(mds.res$eig)*100, 1)

##### Convert data into dataframe for ggplot #####
mds_data = data.frame(sample_id = row.names(mds.res$points), 
                      X=mds.res$points[,1], 
                      Y=mds.res$points[,2]) %>% 
  left_join(sample_info) # include sample's details 

## plot PCoA ## 
fig.4b <- mds_data %>% 
  ggplot(aes(x = X, y = Y, color = y, shape = sex)) + 
  geom_point(alpha = 0.8) +
  labs(title = "Biomarkers + age", 
       x = paste("PC1 - ", mds_varperc[1], "%", sep = ""), 
       y = paste("PC2 -", mds_varperc[2], "%", sep = "")) + 
  scale_color_manual("diagnosis", values = c("#56B4E9", "#D55E00")) +
  guides(shape = guide_legend(override.aes = list(1)), 
         color = guide_legend(override.aes = list(1))) +
  theme_bw() +
  theme(legend.title = element_text(size = 9), 
        legend.text = element_text(size = ))
  
#print
fig.4b
#save for PDF report 
save(fig.4b, file = "figs/fig.4b")

#print fig.4 
ggarrange(fig.4a, fig.4b, 
          labels = c("a", "b"), 
          ncol = 1, nrow = 2) %>%
  annotate_figure("Figure 4. PCoA using Bray-Curtis distance")

### delete objects to free memory ### 
rm(fig.1, fig.2, fig.3, fig.4a, fig.4b, column_details, matrix_hm, matrix_pcoa, mds_data, mds.res)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### MODELS ####

## data for models ##
#create a data frame with center data and with two outcomes only: benign or cancer 
data <-  tidy_data %>% 
  mutate(y = ifelse(diagnosis == "cancer", "cancer", "benign") %>% factor(levels = c("cancer", "benign"))) %>% 
  group_by(biomarker) %>%
  mutate(norm_lvs = log(levels),
         norm_lvs = ifelse(norm_lvs == "-Inf", -15, norm_lvs), #set -inf values to -15 (below detection levels)
         avg = mean(norm_lvs, na.rm = TRUE), #average of each biomarker to center data
         norm_lvs = norm_lvs-avg, #center observations by subtracting group mean (diagnosis and biomarker)
         age = as.numeric(age)) %>%  
  dplyr:: select(-levels, -avg) %>% ungroup() %>%
  pivot_wider(names_from = biomarker, values_from = norm_lvs)

## extra modification ##
#make sex factor 
data$sex <-factor(data$sex)

### create train and test set ###
set.seed(14, sample.kind = "Rounding")
ind_test <- createDataPartition(data$y, times = 1, p = 0.30, list = FALSE) #20% of the data for testing 
test_set <- data[ind_test,]
train_set <- data[-ind_test,]


#### Models NOT using REG1A and plasma_A19_9 ####

#set x's and y in train set : remove unnecessary columns 
train_x <- train_set %>%dplyr::select (-REG1A, -plasma_CA19_9, -sample_id, -y, -sex, -diagnosis)
train_y <- train_set$y

#### GLM model: was used by authors #### 

#cross validation 
train_control <- trainControl(method = "cv", number = 5)
set.seed(112, sample.kind = "Rounding")
glm_model <- train(x = train_x, y = train_y, # in the original paper sex was not used as predictor 
                   method = "glm", 
                   trControl = train_control, 
                   family = binomial())

#glm prediction
pred_glm <- predict(glm_model, test_set)

#glm accuracy 
acc_glm <- mean(pred_glm == test_set$y)
acc_glm

#other performance parameters 
cm_glm <- confusionMatrix(pred_glm, reference = test_set$y)
cm_glm


#### ### ### ### #### 

#### Random forest ####
set.seed(107, sample.kind = "Rounding")
rf_model <- train(x = train_x, y = train_y, 
                  method = "rf", 
                  tuneGrid = data.frame(mtry = seq(1,20,1)), 
                  importance = TRUE)

#tuning 
ggplot(rf_model) +
  scale_x_continuous(breaks = seq(0,30,2)) + 
  theme_bw()
#best value 
rf_model$bestTune

#predictions 
pred_rf <- predict(rf_model, test_set)
#accuracy 
acc_rf <- mean(pred_rf == test_set$y)
acc_rf

#other performance parameters 
cm_rf <- confusionMatrix(pred_rf, reference = test_set$y)
cm_rf

# predictors importance 
rf_model$finalModel$importance %>% data.frame() %>%
  arrange(desc(MeanDecreaseAccuracy))
#this indicates sex is a bad predictor, so we will remove it 


#### ### ### ###

#### Knn model #### 
train_control <- trainControl(method = "repeatedcv", repeats = 5)
set.seed(7, sample.kind = "Rounding")
knn_model <- train(x = train_x, y = train_y,#remove sex
                   method = "knn", 
                   tuneGrid = data.frame(k = seq(3,50,1)), 
                   trControl = train_control)

#tuning 
ggplot(knn_model, highlight = TRUE) + 
  scale_x_continuous(breaks = seq(0, 50, 5 )) + 
  theme_bw()

knn_model$bestTune

#knn predictions
pred_knn <- predict(knn_model, test_set)
#accuracy 
acc_knn <- mean(pred_knn == test_set$y)
acc_knn
# other performance parameters 
cm_knn <- confusionMatrix(pred_knn, reference = test_set$y)
cm_knn


#### ### ### ### 

### LDA model ####
if(!require(MASS)) install.packages("MASS", repos = "http://cran.us.r-project.org")
library(MASS)
set.seed(145, sample.kind = "Rounding")
lda_model <- train(x = train_x, y = train_y, 
                   method = "lda", 
                   trControl = trainControl(method = "cv"))

#predictions
pred_lda <- predict(lda_model, test_set)
# accuracy 
acc_lda <- mean(pred_lda == test_set$y)
acc_lda

#other performance parameters 
cm_lda <- confusionMatrix(pred_lda, reference = test_set$y)
cm_lda


#### ### ### ### 

#### QDA model #### 
set.seed(150, sample.kind = "Rounding")
qda_model <- train(x = train_x, y = train_y, 
                   method = "qda",
                   trControl = trainControl(method = "cv")) #remove sex

#predictions
pred_qda <- predict(qda_model, test_set)

#accuracy 
acc_qda <- mean(pred_qda == test_set$y)
acc_qda

#other perforamnce parameters 
cm_qda <- confusionMatrix(pred_qda, reference = test_set$y)
cm_qda


#### ### ### ### 

#### Bayesian Generalized Linear Model  #### 
if(!require(arm)) install.packages("arm", repos = "http://cran.us.r-project.org")
library(arm)

set.seed(277, sample.kind = "Rounding")
bayesglm_model <- train(x = train_x, y = train_y, 
                        method = "bayesglm", 
                        trControl = trainControl(method = "repeatedcv"))

#predictions 
pred_bglm <- predict(bayesglm_model, test_set)
#accuracy 
acc_bglm <- mean(pred_bglm == test_set$y)
acc_bglm

#other performance parameters 
cm_bglm <- confusionMatrix(pred_bglm, reference = test_set$y)
cm_bglm


#### ### ### ### 
#### Support Vector Machines with Linear Kernel Model #### 
if(!require(kernlab)) install.packages("kernlab", repos = "http://cran.us.r-project.org")
library(kernlab)
temp_train <- train_x %>% 
  mutate(y = train_y)

set.seed(290, sample.kind = "Rounding")
svm_model <- train(y~., data =temp_train, 
                   method = "svmLinear", 
                   tuneGrid = data.frame(C = seq(1,10,1)))
#tuning 
ggplot(svm_model) +
  theme_bw()
#best 
svm_model$bestTune

#predictions 
temp_test <- test_set %>% dplyr:: select(age, creatinine, LYVE1, REG1B, TFF1)
pred_svm <- predict(svm_model, test_set)

#accuracy
acc_svm <- mean(pred_svm == test_set$y)
acc_svm
#other performance parameters 
cm_svm <- confusionMatrix(pred_svm, reference = test_set$y)
cm_svm


#### ### ### ### #### 

#### ENSEMBLE #### 
pred_ens <- data.frame(pred_glm, pred_bglm, pred_lda, pred_qda, pred_rf, pred_knn, pred_svm) %>% 
  mutate(votes = rowMeans(.=="cancer"), 
         pred_ens = ifelse(votes >= 0.5, "cancer", "benign") %>% 
           factor(levels = c("cancer", "benign"))) %>%
  dplyr::select(-votes)

#accuracy 
acc_ens <- mean(pred_ens$pred_ens == test_set$y)
acc_ens

#other performanc parameters 
cm_ens <- confusionMatrix(pred_ens$pred_ens, reference = test_set$y)
cm_ens

#### results ####
#model's names 
model_names = c("GLM", "Bayesian GLM", "LDA", "QDA", "Random Forest",
                "KNN", "SVM Linear Kernel", "Ensemble:Consensus")
#final results ##
model_perf <- data.frame(Model = character(), 
                         Accuracy = numeric(), 
                         Sensitivity = numeric(), 
                         Specificity = numeric())

#fill table 
for(i in 1:length(model_names)) {
  model_perf[i,] <- data.frame(Model = model_names[i], 
                               Accuracy = confusionMatrix(pred_ens[,i], reference = test_set$y)$overall[1] %>% as.numeric(),
                               Sensitivity = confusionMatrix(pred_ens[,i], reference = test_set$y)$byClass[1] %>% as.numeric(), 
                               Specificity = confusionMatrix(pred_ens[,i], reference = test_set$y)$byClass[2] %>% as.numeric())
}

#Arrange results 
model_perf <- model_perf %>%
  arrange(desc(Accuracy))
#print
model_perf

#save table for PDF report 
save(model_perf, file = "rda/model_perf.rda")


#### ### ### ### #### ### ### ### ### #### ### ### ### #### ### ### ### ###

#### ### ### ### #### ### ### ### ### #### ### ### ### #### ### ### ### ###

#### MODEL USING ALL biomarkers  ########

#filter out entries with NA values 
data2 <- data %>% 
  filter(!is.na(REG1A) & !is.na(plasma_CA19_9))

### split data and create test and train data ### 
set.seed(533, sample.kind = "Rounding")
ind_test <- createDataPartition(data2$y, times = 1, p = 0.30, list = FALSE)
test_set <- data2[ind_test,] %>% dplyr::select(-sex, -sample_id)
train_x <- data2[-ind_test,] %>% dplyr::select(-sex, -sample_id, -y, -diagnosis)
train_y <- data2[-ind_test,]$y

#### ### ### ### #### 

#### GLM MODEL #### 
set.seed(377, sample.kind = "Rounding")
model_glm <- train(x = train_x, y =train_y,
                   method = "glm", 
                   trControl = trainControl(method = "cv", number = 5), 
                   family = binomial()) 
#predictions
pred_glm <- predict(model_glm, test_set)

#performance 
cm_glm <- confusionMatrix(pred_glm, reference = test_set$y)
cm_glm

#### ### ### ### ####

#### Bayesian GLM #### 
set.seed(289, sample.kind = "Rounding")
model_bglm <- train(x = train_x, y = train_y,
                    method = "bayesglm")
#predictions 
pred_bglm <- predict(model_bglm, test_set)
#performance 
cm_bglm <- confusionMatrix(pred_bglm, reference = test_set$y)
cm_bglm

#### Random Forest #### 
set.seed(398, sample.kind = "Rounding")
model_rf <- train(x = train_x, y = train_y, 
                  method = "rf", 
                  tuneGrid = data.frame(mtry = seq(1,20,1)), 
                  importance = TRUE)
#tuning 
ggplot(model_rf) + 
  theme_bw()

#predictions 
pred_rf <- predict(model_rf, test_set)
#performance 
cm_rf <- confusionMatrix(pred_rf, reference = test_set$y)
cm_rf

#### ### ### ### ####

#### LDA model #### 
set.seed(414, sample.kind = "Rounding")
model_lda <- train(x = train_x, y = train_y, 
                   method = "lda")
#prediction 
pred_lda <- predict(model_lda, test_set)

#performance 
cm_lda <- confusionMatrix(pred_lda, reference = test_set$y)
cm_lda

#### ### ### ### ####

#### QDA MODEL #### 
set.seed(425, sample.kind = "Rounding")
model_qda <- train(x =train_x, y = train_y, 
                   method = "qda")
#predictions 
pred_qda <- predict(model_qda, test_set)

#performance 
cm_qda <- confusionMatrix(pred_qda, test_set$y)
cm_qda

#### ### ### ### ####

#### KNN model ####
set.seed(436, sample.kind = "Rounding")
model_knn <- train(x = train_x, y =train_y, 
                   method = "knn", 
                   tuneGrid = data.frame(k = seq(1,50, 1)))
#tuning 
ggplot(model_knn, highlight = TRUE) + 
  theme_bw()
#predictions
pred_knn <- predict(model_knn, test_set)
#performance 
cm_knn <- confusionMatrix(pred_knn, test_set$y)
cm_knn

#### ### ### ### ####

#### SVM linear kernel model ####
temp_train <- train_x %>% 
  mutate(y = train_y)

set.seed(450, sample.kind = "Rounding")
model_svm <- train(y~., data =temp_train, 
                   method = "svmLinear", 
                   tuneGrid = data.frame(C = seq(1,10,1)))
#tuning 
ggplot(model_svm, highlight = TRUE)
#prediction
pred_svm <- predict(model_svm, test_set)
#performance
cm_svm <- confusionMatrix(pred_svm, test_set$y)
cm_svm


#### ### ### ### ####

#### ENSEMBLE ### 
#predictions
pred_ens_i <- data.frame(pred_glm, pred_bglm, pred_lda, pred_qda, pred_rf, pred_knn, pred_svm) %>%
  mutate(votes = rowMeans(.=="cancer"), 
         pred_ens = ifelse(votes >= 0.5, "cancer", "benign") %>% factor(levels = c("cancer", "benign"))) %>%
  dplyr::select(-votes)
#performance 
cm_ens <- confusionMatrix(pred_ens_i$pred_ens, test_set$y)
cm_ens

#### results ####
#model's names 
model_names = c("GLM", "Bayesian GLM", "LDA", "QDA", "Random Forest",
                "KNN", "SVM Linear Kernel", "Ensemble:Consensus")
#final results ##
model_i_perf <- data.frame(Model = character(), Accuracy_i = numeric(), 
                          Sensitivity_i = numeric(), 
                          Specificity_i = numeric())

#fill table 
for(i in 1:length(model_names)) {
  model_i_perf[i,] <- data.frame(Model = model_names[i], 
                                 Accuracy_i = confusionMatrix(pred_ens_i[,i], reference = test_set$y)$overall[1] %>% as.numeric(), 
                                 Sensitivity_i = confusionMatrix(pred_ens_i[,i], reference = test_set$y)$byClass[1] %>% as.numeric(), 
                                 Specificity_i = confusionMatrix(pred_ens_i[,i], reference = test_set$y)$byClass[2] %>% as.numeric())
  }
#Arrange results 
model_i_perf <- model_i_perf %>%
  arrange(desc(Accuracy_i))
#print
model_i_perf

#save model for PDF report 
save(model_i_perf, file = "rda/model_i_perf.rda")
