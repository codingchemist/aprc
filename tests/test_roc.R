#############################################################
# Script for testing implementation of performance measures
#############################################################

# load required library
library(boot)
library(glmnet)
library(ROCR)

# set path
setwd("C:/Users/Bian/Documents/code/r/machine_learning/performance_measures/")

# load data set
kcnq_df <- read.csv("dataset.csv", header = TRUE, stringsAsFactors = FALSE)

# removed the first three non-feature columns
kcnq_df <- kcnq_df[-c(1,2,3)]

number_replicates <- 10
bootstrap_results <- data.frame(numeric(nrow(kcnq_df)))
for (i in 1:number_replicates) {
  # bootstrapping
  bootstrap_indices <- sample(1:nrow(kcnq_df), nrow(kcnq_df), replace = TRUE)
  bootstrap_set <- kcnq_df[bootstrap_indices, ]
  # create feature matrix and response vector, the responses are assumed to be last column
  feature_names <- names(bootstrap_set)
  feature_names <- feature_names[-length(feature_names)]
  formula <- as.formula(paste(names(bootstrap_set)[ncol(bootstrap_set)], "~", paste(feature_names, collapse = "+")))
  # boostrap sample
  features <- model.matrix(formula, data = bootstrap_set)[, -ncol(bootstrap_set)]
  responses <- bootstrap_set[, ncol(bootstrap_set)]
  # fit logistic ridge regression model via cross-validation
  cv_model <- cv.glmnet(features, responses, alpha = 0, type.measure = "mse", 
                                 family = "binomial", nfolds = 10)
  # optimal lambda
  optimal_lambda <- cv_model$lambda.min
  # refit the ridge model to the whole data
  refit_model <- glmnet(x = features, y = responses, alpha = 0, family = "binomial")
  # compute predictions
  predictions <- predict(object = refit_model, type = "response", s = optimal_lambda, newx = features)
  # current replicate
  cur_results<- data.frame(predictions, responses)
  cur_results <- cur_results[order(cur_results$X1, decreasing = TRUE), ]
  bootstrap_results <- cbind(bootstrap_results, cur_results)
}