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

#adjust R options 
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
tidy_data <- raw_data %>% dplyr::select(-stage, -patient_cohort, -sample_origin, -stage, -benign_sample_diagnosis) %>% 
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

# statistical analysis #
sex_diag <- table(raw_data$sex, raw_data$diagnosis)
chisq.test(sex_diag)

# it seems that male are more prone to get cancer 

## distribution of diagnosis by age ##
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



### MODELS ###
#create a data frame with center data  
data <-  tidy_data %>% 
  group_by(biomarker) %>%
  mutate(norm_lvs = log(levels),
         norm_lvs = ifelse(norm_lvs == "-Inf", -15, norm_lvs), #set -inf values to -15 (large change)
         avg = mean(norm_lvs, na.rm = TRUE), #average of each biomarker to center data
         norm_lvs = norm_lvs/avg) %>% #center observations %>% 
  dplyr:: select(-levels, -avg) %>% 
  pivot_wider(names_from = biomarker, values_from = norm_lvs) #final shape

## extra modification ##
#make sex factor 
data$sex <-factor(data$sex)

# Add 'y' column with the diagnosis: benign or cancer
data <- data %>%
  mutate(y = ifelse(diagnosis == "cancer", "cancer", "benign") %>% factor(levels = c("cancer", "benign"))) %>% 
  dplyr::select(-diagnosis)

# create train and test set 
set.seed(14, sample.kind = "Rounding")
ind_test <- createDataPartition(data$y, times = 1, p = 0.20, list = FALSE) #20% of the data for testing 
test_set <- data[ind_test,]
train_set <- data[-ind_test,]


### Models not using REG1A and plasma_A19_9 ### 

#set x's and y in train set 
train_x <- train_set %>%dplyr::select (-REG1A, -plasma_CA19_9, -sample_id, -y)
train_y <- train_set$y

#### GLM model: was used by authors #### 
set.seed(112, sample.kind = "Rounding")
glm_model <- train(x = train_x[,-2], y = train_y, # in the original paper sex was not used as predictor 
                   method = "glm")

#glm prediction
pred_glm1 <- predict(glm_model, test_set)

#glm accuracy 
glm1_acc <- mean(pred_glm1 == test_set$y)
glm1_acc

#other performance parameters 
cm_glm <- confusionMatrix(pred_glm1, reference = test_set$y)
cm_glm

#save results 
res_glm <- tibble(Method = "GLM", Accuracy = glm1_acc, 
                   Sensitivity = cm_glm$byClass[1], 
                   Specificity = cm_glm$byClass[2])

#### GLM including sex as predictor ### 
set.seed(133, sample.kind = "Rounding")
glm_model2 <- train(x = train_x, y = train_y, 
                    method = "glm")
#glm 2 prediction 
pred_glm2 <- predict(glm_model2, test_set)
#glm2 predictions 
glm2_acc <- mean(pred_glm2 == test_set$y)
glm2_acc
#other performance parameters 
confusionMatrix(pred_glm2, reference = test_set$y)

#### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### 

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
rf_acc <- mean(pred_rf == test_set$y)

#other performance parameters 
confusionMatrix(pred_rf, reference = test_set$y)

# predictors importance 
rf_model$finalModel$importance %>% data.frame() %>%
  arrange(desc(MeanDecreaseAccuracy))
#this indicates sex is a bad predictor, so we will remove it 

#### repeat model without sex #### 

### Random forest ####
set.seed(107, sample.kind = "Rounding")
rf_model <- train(x = train_x[,-2], y = train_y, #remove sex
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
pred_rf2 <- predict(rf_model, test_set)
#accuracy 
rf_acc2 <- mean(pred_rf2 == test_set$y)

#other performance parameters 
cm_rf <- confusionMatrix(pred_rf2, reference = test_set$y)
cm_rf

#results 
res_rf <- tibble(Method = "Random Forest", Accuracy = rf_acc2, 
                 Sensitivity = cm_rf$byClass[1], 
                 Specificity = cm_rf$byClass[2])

#### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### 
#### Knn model #### 
set.seed(7, sample.kind = "Rounding")
knn_model <- train(x = train_x[,-2], y = train_y,#remove sex
                   method = "knn", tuneGrid = data.frame(k = seq(3,50,1)))

#tuning 
ggplot(knn_model, highlight = TRUE) + 
  scale_x_continuous(breaks = seq(0, 50, 5 )) + 
  theme_bw()

knn_model$bestTune

#knn predictions
pred_knn <- predict(knn_model, test_set)
#accuracy 
knn_acc <- mean(pred_knn == test_set$y)
knn_acc
# other performance parameters 
cm_knn <- confusionMatrix(pred_knn, reference = test_set$y)
cm_knn

#save results 
res_knn <- tibble(Method = "KNN", Accuracy = knn_acc, 
                  Sensitivity = cm_knn$byClass[1], 
                  Specificity = cm_knn$byClass[2])

#### ### ### ### #### ### ### ### ### #### ### ### ### #### ### ### ### ### 
### LDA model ####
if(!require(MASS)) install.packages("MASS", repos = "http://cran.us.r-project.org")
library(MASS)
set.seed(145, sample.kind = "Rounding")
lda_model <- train(x = train_x[,-2], y = train_y, method = "lda") #remove sex

#predictions
pred_lda <- predict(lda_model, test_set)
# accuracy 
lda_acc <- mean(pred_lda == test_set$y)
lda_acc

#other performance parameters 
cm_lda <- confusionMatrix(pred_lda, reference = test_set$y)
cm_lda

#save results 
res_lda <-  tibble(Method = "LDA", Accuracy = lda_acc, 
                   Sensitivity = cm_lda$byClass[1], 
                   Specificity = cm_lda$byClass[2]) 

#### ### ### ### #### ### ### ### ### #### ### ### ### #### ### ### ### ### 
#### QDA model #### 
set.seed(150, sample.kind = "Rounding")
qda_model <- train(x = train_x[,-2], y = train_y, method = "qda") #remove sex

#predictions
pred_qda <- predict(qda_model, test_set)

#accuracy 
qda_acc <- mean(pred_qda == test_set$y)
qda_acc

#other perforamnce parameters 
cm_qda <- confusionMatrix(pred_qda, reference = test_set$y)
cm_qda

#save results 
res_qda <- tibble(Method = "QDA", Accuracy = qda_acc, 
                  Sensitivity = cm_qda$byClass[1], 
                  Specificity = cm_qda$byClass[2])

#### ### ### ### #### ### ### ### ### #### ### ### ### #### ### ### ### ###  

#### Bayesian Generalized Linear Model  #### 
if(!require(arm)) install.packages("arm", repos = "http://cran.us.r-project.org")
library(arm)

set.seed(277, sample.kind = "Rounding")
bayesglm_model <- train(x = train_x[,-2], y = train_y, method = "bayesglm")

#predictions 
pred_bglm <- predict(bayesglm_model, test_set)
#accuracy 
bglm_acc <- mean(pred_bglm == test_set$y)

#other performance parameters 
cm_bglm <- confusionMatrix(pred_bglm, reference = test_set$y)
cm_bglm

#save results 
res_bglm <- tibble(Method = "Bayesian GLM", Accuracy = bglm_acc, 
                   Sensitivity = cm_bglm$byClass[1], 
                   Specificity = cm_bglm$byClass[2])

#### ### ### ### #### ### ### ### ### #### ### ### ### #### ### ### ### ###
#### Support Vector Machines with Linear Kernel Model #### 
if(!require(kernlab)) install.packages("kernlab", repos = "http://cran.us.r-project.org")
library(kernlab)
temp_train <- train_x %>% dplyr::select(-sex) %>%
  mutate(y = train_y)

set.seed(290, sample.kind = "Rounding")
svm_model <- train(y~., data =temp_train, 
                   method = "svmLinear", 
                   tuneGrid = data.frame(C = seq(1,10,1)))
#tuning 
ggplot(svm_model) +
  theme_bw()

#predictions 
temp_test <- test_set %>% dplyr:: select(age, age, creatinine, LYVE1, REG1B, TFF1)
pred_svm <- predict(svm_model, test_set)

#accuracy
svm_acc <- mean(pred_svm == test_set$y)

#other performance parameters 
cm_svm <- confusionMatrix(pred_svm, reference = test_set$y)
cm_svm

#save results 
res_svm <- tibble(Method = "SVM Linear Kernel", Accuracy = svm_acc, 
                  Sensitivity = cm_svm$byClass[1], 
                  Specificity = cm_svm$byClass[2])


#### ENSEMBLE #### 
pred_ens <- tibble(pred_glm1, pred_rf2, pred_knn, pred_lda, pred_qda, 
                   pred_bglm, pred_svm) %>% 
  mutate(votes = rowMeans(.=="cancer"), 
         pred_ens = ifelse(votes >= 0.5, "cancer", "benign") %>% factor(levels = c("cancer", "benign")))

#accuracy 
ens_acc <- mean(pred_ens$pred_ens == test_set$y)
ens_acc

#other performanc parameters 
cm_ens <- confusionMatrix(pred_ens$pred_ens, reference = test_set$y)
cm_ens

#save 
res_ens <- tibble(Method = "Ensemble:Consensus", Accuracy = ens_acc, 
                  Sensitivity = cm_ens$byClass[1], 
                  Specificity = cm_ens$byClass[2])

### Models performance ####
mod_perd <- res_bglm %>%
  rbind(res_glm) %>%
  rbind(res_knn) %>% 
  rbind(res_lda) %>% 
  rbind(res_qda) %>%
  rbind(res_rf) %>% 
  rbind(res_svm) %>%
  arrange(desc(Accuracy)) %>%
  rbind(res_ens)

#### MODEL USING other biomarkers  ########
#filter out entries with NA values 
data2 <- data %>% 
  filter(!is.na(REG1A) & !is.na(plasma_CA19_9))

### split data and create test and train data ### 
ind_test <- createDataPartition(data2$y, times = 1, p = 0.30, list = FALSE)
test_set <- data2[ind_test,] %>% dplyr::select(-sex, -sample_id)
train_x <- data2[-ind_test,] %>% dplyr::select(-sex, -sample_id, -y)
train_y <- data2[-ind_test,]$y

#### GLM MODEL #### 
set.seed(377, sample.kind = "Rounding")
model_glm <- train(x = train_x, y =train_y,
                   method = "glm") 
#predictions
pred_glm <- predict(model_glm, test_set)

#performance 
cm_glm <- confusionMatrix(pred_glm, reference = test_set$y)
cm_glm


#### Bayesian GLM #### 
set.seed(289, sample.kind = "Rounding")
model_bglm <- train(x = train_x, y = train_y, method = "bayesglm")
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

#### LDA model #### 
set.seed(414, sample.kind = "Rounding")
model_lda <- train(x = train_x, y = train_y, 
                   method = "lda")
#prediction 
pred_lda <- predict(model_lda, test_set)

#performance 
cm_lda <- confusionMatrix(pred_lda, reference = test_set$y)
cm_lda

#### QDA MODEL #### 
set.seed(425, sample.kind = "Rounding")
model_qda <- train(x =train_x, y = train_y, 
                   method = "qda")
#predictions 
pred_qda <- predict(model_qda, test_set)

#performance 
cm_qda <- confusionMatrix(pred_qda, test_set$y)
cm_qda

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

#### ENSEMBLE ### 
#predictions
pred_ens_i <- data.frame(pred_glm, pred_bglm, pred_lda, pred_qda, 
                       pred_rf, pred_knn, pred_svm) %>%
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



### analysis of errors ### 
test_sex <- data2[ind_test,]$sex
test_y <- test_set$y

errors <- pred_ens_i %>% 
  mutate(sex = test_sex, 
         y = test_y)
