#' @brief Implementation of performance measures of classifiers
#' @author Bian Li
#' @e-mail bian.li@vanderbilt.edu
#' @date 9/30/2015

#' @TODO add support for the case when predictions is a data frame
predeval <- function(predictions, labels, positive.label, prediction.is.numeric = FALSE) {
  # Constructor of predeval object
  # Parameters
  # ----------
  # predictions: numeric vector
  #    Predictions (scores) produced by a classifier
  # labels: factor with 2 levels
  #    True class labels
  # Returns
  # -------
  # A predeval object
	if (length(predictions) != length(labels)) {
		stop("predictions and responses must have equal lengths!")
	} else if (length(levels(as.factor(labels))) != 2) {
		stop("labels must be a factor with 2 levels!")
	}
	
	# number of positive instances and number of negative instances
	labels <- as.factor(labels)
	negative.label <- ifelse(levels(labels)[1] == positive.label, 
			levels(labels)[2], levels(labels)[1])  
	p <- sum(labels == positive.label)
	n <- sum(labels == negative.label)
  
	# compute tp, tn, fp, fn
	if (!prediction.is.numeric) {
		# in this scenario, tp, tn, etc. are scalars
		tp <- sum(predictions == positive.label & predictions == labels)
		tn <- sum(predictions == negative.label & predictions == labels)
		roc <- NULL
	} else {
		# in this scenario, tp, tn, etc. are vectors
		roc <- compute_roc(predictions, labels)
		tp <- roc$tps
		tn <- n - roc$fps
	}
	
	fp <- n - tn
	fn <- p - tp
  
	predeval <- list(
			predictions = predictions,
			labels = labels,
			positive.label = positive.label,
			negative.label = negative.label,
			p = p,
			n = n,
			tp = tp,
			tn = tn,
			fp = fp,
			fn = fn,
			roc = roc
		)

  
	# set the name of the class
	class(predeval) <- c("predeval")
	return(predeval)
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

#get_sensitivity <- function(c.eval.obj) {
#	UseMethod("get_sensitivity", c.eval.obj)
#}

#' @note 
#' @param an object of predeval
#' @return specificity
get_tpr <- function(predeval.obj) {
	tpr <- predeval.obj$tp / predeval.obj$p
	if (is.null(predeval.obj$roc)) {
		return(tpr)
	} else {
		cutoff <- predeval.obj$roc$cutoffs
		return(data.frame(cutoff, tpr))
	}
}

#get_specificity <- function(predeval.obj) {
#	UseMethod("get_specificity", predeval.obj)
#}

#' @note 
#' @param an object of predeval
#' @return specificity
get_tnr <- function(predeval.obj) {
	tnr <- predeval.obj$tn / predeval.obj$n
	if (is.null(predeval.obj$roc)) {
		return(tnr)
	} else {
		cutoff <- predeval.obj$roc$cutoffs
		return(data.frame(cutoff, tnr))
	}
}

#' @note 
#' @param an object of predeval
#' @return specificity
get_ppv <- function(predeval.obj) {
	ppv <- predeval.obj$tp / (predeval.obj$tp + predeval.obj$fp)
	if (is.null(predeval.obj$roc)) {
		return(ppv)
	} else {
		cutoff <- predeval.obj$roc$cutoffs
		return(data.frame(cutoff, ppv))
	}
}

#'
#' 
get_npv <- function(predeval.obj) {
	npv <- predeval.obj$tn / (predeval.obj$tn + predeval.obj$fn)
	if (is.null(predeval.obj$roc)) {
		return(npv)
	} else {
		cutoff <- predeval.obj$roc$cutoffs
		return(data.frame(cutoff, npv))
	}
}

get_accuracy <- function(predeval.obj) {
	accuracy <- (predeval.obj$tp + predeval.obj$tn) / (predeval.obj$p + predeval.obj$n)
	if (is.null(predeval.obj$roc)) {
		return(accuracy)
	} else {
		cutoff <- predeval.obj$roc$cutoffs
		return(data.frame(cutoff, accuracy))
	}
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

get_hamm_dist <- function(predeval.obj) {
	UseMethod("get_hamm_dist", predeval.obj)
}

#'
#' 
get_hamm_dist.predeval <- function(predeval.obj) {
  predictions <- predeval.obj$m_predictions
  labels <- predeval.obj$m_labels
  cutoff <- predeval.obj$m_roc$cutoffs
  hamm_dist <- c()
  for (c in cutoff) {
    hamm_dist <- c(hamm_dist, compute_hamm_dist(predictions, labels, c))
  }
  return(data.frame(cutoff, hamm_dist))
}

get_quad_dist <- function(predeval.obj) {
  UseMethod("get_quad_dist", predeval.obj)
}

#' this implementation of quadratic distance assumes predictions to be numeric
#' and no cutoffs is needed
get_quad_dist.predeval <- function(predeval.obj) {
  predictions <- predeval.obj$m_predictions
  labels <- predeval.obj$m_labels
  pos_label <- levels(labels)[2]
  neg_label <- levels(labels)[1]
  labels <- ifelse(labels == pos_label, 1, 0)
  return(sum((predictions - labels)^2))
}

#' @note See Pierre Baldi, et al., Bioinformatics, 2000, 412-424
#' for a derivation of the formula for computing mcc
get_mcc <- function(predeval.obj) {
	tp <- predeval.obj$tp
	tn <- predeval.obj$tn
	fp <- predeval.obj$fp
	fn <- predeval.obj$fn
	mcc <- (tp * tn - fp * fn) / sqrt((tp + fn) * (tp + fp) * (tn + fp) * (tn + fn))
	if(is.null(predeval.obj$roc)) {
		return(mcc)
	} else {
		cutoff <- predeval.obj$m_roc$cutoffs
		return(data.frame(cutoff, mcc))		
	}
}

entropy_log <- function(x) {
	# logarithm for entropy calculation
	return(ifelse(x == 0, 0, log(x)))
}

get_cross_entropy<- function(predeval.obj) {
  #   References
  #   ----------
  #   C.M. Bishop (2006). Pattern Recognition and Machine Learning. Springer,
  #   p. 206.
	labels <- ifelse(predeval.obj$labels == predeval.obj$positive.label, 1, 0)
	if (is.null(predeval.obj$roc)) {
		predictions <- ifelse(predeval.obj$predictions == predeval.obj$positive.label, 1, 0)
	} else {
		predictions <- predeval.obj$predictions
	}
	cross_entropy <- sum(-(labels * entropy_log(predictions) + (1 - labels) * entropy_log(1 - predictions)))
	return(cross_entropy)
}

#' @brief Prints a predeval object
#' 
print.predeval <- function(x, ...) {
	if(is.null(x$roc)) {
		measure.names <- c("p", "n", "tp", "fp", "tn", "fn")
		measure.values <- c(x$p, x$n, x$tp, x$fp, x$tn, x$fn)
		print(data.frame(measure.names, measure.values))
	} else {
		print(x$roc)
	}
}