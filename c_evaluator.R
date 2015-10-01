##########################################################################
# @brief Implementation of performance measures of classifiers
# @author Bian Li
# @e-mail bian.li@vanderbilt.edu
# @date 9/30/2015
##########################################################################

c_evaluator <- function(predictions, labels) {
  # Constructor of c_evaluator object
  # Parameters
  # ----------
  # predictions: numeric vector
  #    Predictions (scores) produced by a classifier
  # labels: factor with 2 levels
  #    True class labels
  # Returns
  # -------
  # A c_evaluator object
  if (length(predictions) != length(labels)) {
    stop("predictions and responses must have equal lengths!")
  } else if (!is.vector(predictions, "numeric")) {
    stop("predictions must be a numeric vector!")
  } else if (length(levels(as.factor(labels))) != 2) {
    stop("labels must either be a factor with 2 levels!")
  }
  
  # number of positive instances and number of negative instances
  labels <- as.factor(labels)
  pos_label <- levels(labels)[2]
  neg_label <- levels(labels)[1]
  p <- sum(labels == pos_label)
  n <- sum(labels == neg_label)
  
  
  data_members <- list(
    m_predictions = predictions,
    m_labels = labels,
    m_pos = p,
    m_neg = n,
    m_roc = compute_roc(predictions, labels)
  )
  
  # set the number of the class
  class(data_members) <- c("c_evaluator")
  return(data_members)
}


compute_roc <- function(predictions, labels) {
  # Compute #tp, #tn, #fp, #fn at every possible cutoff.
  # Algorithm
  # ---------
  #
  # Parameters
  # ----------
  # predictions: numeric vector
  #  Predictions (scores) produced by a classifier
  # labels: factor with 2 levels
  #  True class labels
  pred_order <- order(predictions, decreasing = TRUE)
  pred_sorted <- predictions[pred_order]
  labels_sorted <- labels[pred_order]
  
  # positive label and negative label
  pos_label <- levels(labels)[2]
  neg_label <- levels(labels)[1]
  
  tps <- fps <- cutoffs <- c()
  tp <- fp <- 0
  score_prev <- Inf
  for (i in 1:length(pred_sorted)) {
    if (labels_sorted[i] == pos_label) {
      tp <- tp + 1
    } else {
      fp <- fp + 1
    }
    if (pred_sorted[i] != score_prev) {
      tps <- c(tps, tp)
      fps <- c(fps, fp)
      cutoffs <- c(cutoffs, pred_sorted[i])
      score_prev <- pred_sorted[i]
    }
  }
  return(data.frame(cutoffs, tps, fps))
}

##########################################################
# Implemeting Performance Measures Discussed in
# Pierre Baldi, et al., Bioinformatics, 2000, 412-424
##########################################################

get_sensitivity <- function(c_evaluator.obj) {
  UseMethod("get_sensitivity", c_evaluator.obj)
}

get_sensitivity.c_evaluator <- function(c_evaluator.obj) {
  sensitivity <- c_evaluator.obj$m_roc$tps / c_evaluator.obj$m_pos
  cutoff <- c_evaluator.obj$m_roc$cutoffs
  return(data.frame(cutoff, sensitivity))
}

get_specificity <- function(c_evaluator.obj) {
  UseMethod("get_specificity", c_evaluator.obj)
}

get_specificity.c_evaluator <- function(c_evaluator.obj) {
  specificity <- (c_evaluator.obj$m_neg - c_evaluator.obj$m_roc$fps) / c_evaluator.obj$m_neg
  cutoff <- c_evaluator.obj$m_roc$cutoffs
  return(data.frame(cutoff, specificity))
}

get_ppv <- function(c_evaluator.obj) {
  UseMethod("get_ppv", c_evaluator.obj)
}

get_ppv.c_evaluator <- function(c_evaluator.obj) {
  ppv <- c_evaluator.obj$m_roc$tps / (c_evaluator.obj$m_roc$tps + c_evaluator.obj$m_roc$fps)
  cutoff <- c_evaluator.obj$m_roc$cutoffs
  return(data.frame(cutoff, ppv))
}

get_npv <- function(c_evaluator.obj) {
  UseMethod("get_npv", c_evaluator.obj)
}

get_npv.c_evaluator <- function(c_evaluator.obj) {
  tns <- c_evaluator.obj$m_neg - c_evaluator.obj$m_roc$fps
  fns <- c_evaluator.obj$m_pos - c_evaluator.obj$m_roc$tpr
  npv <- tns / (tns + fns)
  cutoff <- c_evaluator.obj$m_roc$cutoffs
  return(data.frame(cutoff, npv))
}

compute_hamm_dist <- function(predictions, labels, cutoff = 0.5) {
  # An alternative implementaton would be count the total number of errors
  # which is equal to fp + fn, this approach is not implemented here though.
  if (length(predictions) != length(labels)) {
    stop("predictions and labels must have the same number of entries!")
  } else if (!is.factor(labels)) {
    labels <- as.factor(labels)
  }
  pos_label <- levels(labels)[2]
  neg_label <- levels(labels)[1]
  labels <- ifelse(labels == pos_label, 1, 0)
  # convert predictions into binary
  if (is.numeric(predictions)) {
    predictions <- ifelse(predictions >= cutoff, 1, 0)
  }  else if (!is.factor(predictions)) {
    predictions <- as.factor(predictions)
    predictions <- ifelse(predictions == pos_label, 1, 0)
  }
  # compute hamming distance
  return(sum(abs(predictions - labels)))
}

get_hamm_dist <- function(c_evaluator.obj) {
  UseMethod("get_hamm_dist", c_evaluator.obj)
}

get_hamm_dist.c_evaluator <- function(c_evaluator.obj) {
  predictions <- c_evaluator.obj$m_predictions
  labels <- c_evaluator.obj$m_labels
  cutoff <- c_evaluator.obj$m_roc$cutoffs
  hamm_dist <- c()
  for (c in cutoff) {
    hamm_dist <- c(hamm_dist, compute_hamm_dist(predictions, labels, c))
  }
  return(data.frame(cutoff, hamm_dist))
}

get_quad_dist <- function(c_evaluator.obj) {
  UseMethod("get_quad_dist", c_evaluator.obj)
}

get_quad_dist.c_evaluator <- function(c_evaluator.obj) {
  # this implementation of quadratic distance assumes predictions to be numeric
  # and no cutoffs is needed
  predictions <- c_evaluator.obj$m_predictions
  labels <- c_evaluator.obj$m_labels
  pos_label <- levels(labels)[2]
  neg_label <- levels(labels)[1]
  labels <- ifelse(labels == pos_label, 1, 0)
  return(sum((predictions - labels)^2))
}

get_mcc <- function(c_evaluator.obj) {
  UseMethod("get_mcc", c_evaluator.obj)
}

get_mcc.c_evaluator <- function(c_evaluator.obj) {
  # See Pierre Baldi, et al., Bioinformatics, 2000, 412-424
  # for a derivation of the formula for computing mcc
  tp <- c_evaluator.obj$m_roc$tps
  fp <- c_evaluator.obj$m_roc$fps
  tn <- c_evaluator.obj$m_neg - fp
  fn <- c_evaluator.obj$m_pos - tp
  mcc <- (tp * tn - fp * fn) / sqrt((tp + fn) * (tp + fp) * (tn + fp) * (tn + fn))
  cutoff <- c_evaluator.obj$m_roc$cutoffs
  return(data.frame(cutoff, mcc))
}

get_cross_entropy <- function(c_evaluator.obj) {
  UseMethod("get_cross_entropy", c_evaluator.obj)
}

entropy_log <- function(x) {
  # logarithm for entropy calculation
  ifelse(x == 0, 0, log(x))
}

get_cross_entropy.c_evaluator <- function(c_evaluator.obj) {
  #   References
  #   ----------
  #   C.M. Bishop (2006). Pattern Recognition and Machine Learning. Springer,
  #   p. 206.
  predictions <- c_evaluator.obj$m_predictions
  labels <- c_evaluator.obj$m_labels
  pos_label <- levels(labels)[2]
  neg_label <- levels(labels)[1]
  labels <- ifelse(labels == pos_label, 1, 0)
  cross_entropy <- sum(-(labels * entropy_log(predictions) + (1 - labels) * entropy_log(1 - predictions)))
  return(cross_entropy)
}