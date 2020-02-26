# ---------------------Dat Mining Group Project -----------------------------
# Time series, 2007-2015 lending club loan dataset
#
# Jianhong Chen
# 10-22-2019
# ----------------------------------------------------------------------------
# intial set up
ls() # list all objects 
gc() # force R to release memory of deleted variables 
rm(list = ls()) # delete all variables in the environment
# ----- load require libraries
# reading datafiles
library(readxl) # read excel files 
library(foreign) # read Weka arff files

# data wrangling 
library(tidyverse) # inlcude "tidyr", "readr", "dplyr", "purrr"
library(reshape2) 
library(lubridate) # handling date variable 

# "big" data 
library(bigmemory)
library(biganalytics)
library(bigtabulate)

# classifier packages
library(caret)
library(dbscan) # kNN
library(factoextra) # multivariable analysis in R [hkmean]
library(e1071) # naiveBayes 
library(rpart) # decision tree
library(ROCR) # better evalulize the result of the classifier models
library(pROC)
library(kernlab) # kernal method 

# plotting packages
library(RColorBrewer)
library(ggplot2)
library(GGally)

# database packages for SQL
library(dplyr)
library(dbplyr)
library(RSQLite)
# -----------------------------------------------------------------------------------------------------

dir = "/Users/jianhongchen/Desktop/UCSC_AMS_PhD/Courses_Works/CSE243(DataMining)_fall2019/group_project/lending_club_loan_data"
setwd(dir)
#LC_raw = read_csv("loan.csv")

# conncet to the database
LC_db= DBI::dbConnect(RSQLite::SQLite(), "database.sqlite")

# make a quick small query to study the all the possible features of the data set
LC_quick = tbl(LC_db, "loan") %>%
  head(10) %>%
  show_query() %>%
  collect()
#SQL implementation in R
n = 50000 # total pulled data point
LC_df = tbl(LC_db, "loan") %>%
  #filter(issue_d %in% c("Dec-2015","Aug-2015")) %>%
  dplyr::select( loan_amnt, term, int_rate, installment, grade, sub_grade,
                 emp_title, emp_length, home_ownership, annual_inc, issue_d, 
                loan_status, purpose, addr_state, dti, inq_last_12m, total_pymnt,
                total_rec_late_fee, last_pymnt_amnt, tot_hi_cred_lim) %>%
  head(n,n) %>% 
  show_query() %>%
  collect() %>%
  transform( # if want string type, no change is needed
            loan_amnt = as.numeric(loan_amnt),
            term = as.factor(term),
            int_rate = as.numeric(int_rate),
            installment = as.numeric(installment),
            grade = as.factor(grade),
            sub_grade = as.factor(sub_grade),
            home_ownership = as.factor(home_ownership),
            annual_inc = as.numeric(annual_inc),
            loan_status = as.factor(loan_status),
            dti = as.numeric(dti),
            inq_last_12m = as.factor(inq_last_12m),
            total_pymnt = as.numeric(total_pymnt),
            total_rec_late_fee = as.numeric(total_rec_late_fee),
            last_pymnt_amnt = as.numeric(last_pymnt_amnt),
            tot_hi_cred_lim = as.numeric(tot_hi_cred_lim)
            )

#missmap(LC_df, col = c("blue", "red"), legend = T) # check missing values
# --- convert issue_date into date
date = LC_df$issue_d %>%
  str_split("-")
month = sapply(date, function(x){x[1]})
year = sapply(date, function(x){x[2]})

LC_df = dplyr::select(LC_df, -issue_d) %>%
  mutate("month" = factor(month), "year" = factor(year))

# ---------------------------- build anomaly detection ----------------------------
LC_num = select_if(LC_df, is.numeric) %>%
  scale(center = T, scale = T) %>%
  as.data.frame() %>%  
  mutate("class" = factor(ifelse(LC_df$loan_status %in% 
                                   c("Current", "Fully Paid"), 0, 1))) %>%
  na.omit() # remove missing values from any features

LC_num_long = LC_num %>%
  gather(key = "finance_terms", value= "value", -class)


# ---check the distribution of the financial terms
ggplot(LC_num_long, aes(x = sqrt(abs(value)), fill = class)) +
  stat_bin(bins = 100) + 
  facet_wrap(~finance_terms, nrow = 2) + 
  theme(plot.title = element_text(hjust = 0.5),
                text = element_text(size = 20))

LC_num_tran = sqrt(abs(select(LC_num, -class))) %>%
  mutate("class" = LC_num$class)
  


# ---- split the data into 3 sets: train, Cross Validation(CV), test
SplitTrainCVTest = function(x_df, split_pcent){
  # split the data into 3 random subsets: training, CV, testing
  # input: x_df = input dataframe;
  #        split_pcent = the splitting percentages between train vs (CV+test)
  # output: return a list of dataframe containts: 1. training; 2.CV; 3.testing
  smp_size1 = floor(split_pcent*nrow(x_df))
  train_ind = sample(seq_len(nrow(x_df)), size = smp_size1)
  x_train = x_df[train_ind,]
  x_CVtest = x_df[-train_ind, ]
  
  smp_size2 = floor(0.5*nrow(x_CVtest))
  CVtest_ind = sample(seq_len(nrow(x_CVtest)), size = smp_size2)
  x_CV = x_CVtest[CVtest_ind, ]
  x_test = x_CVtest[-CVtest_ind, ]
  return(list(x_train, x_CV,  x_test)) 
}

set.seed(78)
anom_train = SplitTrainCVTest(LC_num_tran, 0.6)[[1]] 
anom_CV = SplitTrainCVTest(LC_num_tran, 0.6)[[2]]
anom_test = SplitTrainCVTest(LC_num_tran, 0.6)[[3]] 

# ---compute the Gaussian parameters for training features
GaussianFun = function(x, mu, var){
  p = 1/sqrt(2*pi*var) * exp(-(x - mu)^2/(2*var))
  return(p)
}

train_input = select(anom_train, -class)
feature_mu = colMeans(train_input)
feature_var = colSums((train_input - feature_mu)^2) / dim(train_input)[1]

# ---- compute probability of Gaussian distibution for the CV set features
CV_input = select(anom_CV, -class)
p_CV = GaussianFun(x = CV_input, mu = feature_mu, var = feature_var )
M = dim(CV_input)[1]
p_CV_prod = matrix(0, nrow = M)
for(m in 1:M){
  p_CV_prod[m] = prod(p_CV[m, ])
}

# --- find the optimized threshold epsilon by using the CV set
ComputeFMeasure = function(prediction, ground_T){
  # compute the F-measure of a binary classifier
  # precision, recall, accuracy, f1-score
  #
  # inputs: prediction = the predicted class label by the model
  #         ground_T = ground true label from the actual file
  #
  # outputs: A list contains the four F-measure metrics. 
  tab = table(prediction, ground_T)
  
  TP = tab[1,1]
  FN = tab[2,1]
  FP = tab[1,2]
  TN = tab[2,2]
  # precision = true pos. / (total predict pos.)
  prec = TP / (TP + FP)
  # recall = true pos./(total actual pos.)
  recall = TP / (TP+FN) 
  
  accuracy = (TP + TN) /(TP + FN + FP + TN)
  # f1 = 2*(prec. * recall) / (prec. + recall)
  f1 = 2*prec*recall/(prec + recall)
  
  result = list(prec, recall, accuracy, f1)
  names(result) = c("precision", "recall", "accuracy", "f1_score")
  return(result)
}
SelectThreshold = function(p_val, y_lab){
  F1 = 0 
  threshold = 0
  stepsize = (max(p_val) - min(p_val)) / 1000
  for (epsilon in seq(min(p_val)+stepsize, max(p_val)-stepsize, stepsize)){
    anam_pred = ifelse(p_val < epsilon, 1, 0)
    anam_f1 = ComputeFMeasure(anam_pred, y_lab)$f1_score
    if (anam_f1 > F1) {
      F1 = anam_f1
      threshold = epsilon
    }
  }
  return(c(threshold, F1))
}

bestThreshold = SelectThreshold(p_CV_prod, anom_CV$class)[1]
bestF1 = SelectThreshold(p_CV_prod, anom_CV$class)[2]

# --- now test the model with test set
test_input = select(anom_test, -class)
p_test = GaussianFun(x = test_input, mu = feature_mu, var = feature_var)
p_test_prod = matrix(0, nrow = nrow(test_input))
for(m in 1:nrow(test_input)){
  p_test_prod[m] = prod(p_test[m, ])
}

anom_result = ifelse(p_test_prod < bestThreshold, 1, 0)
anom_tab = table(anom_result, anom_test$class)
anom_f1 = ComputeFMeasure(anom_result, anom_test$class)

# --------------------------- test with svm model ----------------------
SplitTrainTest = function(x_df, split_pcent){
  smp_size = floor(split_pcent*nrow(x_df))
  train_ind <- sample(seq_len(nrow(x_df)), size = smp_size)
  x_train = x_df[train_ind,]
  x_test = x_df[-train_ind, ]
  return(list(x_train, x_test)) 
}
svm_train = SplitTrainTest(LC_num, 0.75)[[1]]
svm_test = SplitTrainTest(LC_num, 0.75)[[2]]

svm_model = svm(class ~. , svm_train,
                method="C-classification", kernal="rbfdot",
                gamma=0.1, cost=10)

train_pred = predict(svm_model, svm_train)
train_tab = table(train_pred, svm_train$class)
train_f1 = ComputeFMeasure(train_pred, svm_train$class)

test_pred = predict(svm_model, svm_test)
test_tab = table(test_pred, svm_test$class)
test_f1 = ComputeFMeasure(test_pred, svm_test$class)


# -------------------------------------  clustering -----------------------------------
LC_num_PCA = prcomp(LC_num, scale = T)
PCA_eigen_per = percentVar = round(100*LC_num_PCA$sdev^2 / sum(LC_num_PCA$sdev^2), 1)
ggplot(as.data.frame(LC_num_PCA$x), aes(x = PC1, y = PC2)) + 
  geom_point(shape = 20)

# perform a clustering based on the financial terms 

LC_km = kmeans(LC_num, centers = 5)
#LC_kkm = kkmeans(as.matrix(LC_num), centers = 5, kernel = "rbfdot")

LC_GG = data.frame(LC_num_PCA$x) %>%
  mutate("cluster" = factor(LC_km$cluster))

ggplot(LC_GG, aes(PC1, PC2)) +
  geom_point(alpha = 0.6, shape =20, aes(colour = cluster)) +
  ggtitle("Visualization of k-means") +
  xlab(paste0("PC1, VarExp: ", PCA_eigen_per[1], "%")) +
  ylab(paste0("PC2, VarExp: ", PCA_eigen_per[2], "%")) +
  # ylim(-6.5, 6.5)+
  # xlim(-4,12) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = brewer.pal(8, "Set2"))




dbDisconnect(LC_db) # disconnect the database





