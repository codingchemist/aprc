##########################################################################
# @brief Implementation of the Receiver Operating Characteristic Curve 
# @author Bian Li
# @e-mail bian.li@vanderbilt.edu
# @date 9/29/2015
##########################################################################


roc <- function(fpr, tpr) {
  # Contructor
  # fpr: a numeric vector of false positive rates
  # tpr: a numeric vector of true positive rates
  stopifnot(
    length(fpr) == length(tpr) 
    && is.vector(fpr, "numeric")
    && is.vector(tpr, "numeric")
  )
  
  data_members <- list(
    m_fpr = fpr,
    m_tpr = tpr
  )
  
  # set the name of the class
  class(data_members) <- append(class(data_members), "roc")
  return(data_members)
}

set_fpr <- function(roc, new_fpr) {
  print("Calling the base set_fpr function")
  UseMethod("set_fpr", roc)
}

set_fpr.roc <- function(roc, new_fpr) {
  # Set the m_fpr member of a roc object.
  # roc: a roc object
  # new_fpr: a numeric vector
  stopifnot(
    length(new_fpr) == length(roc$m_fpr)
    && is.vector(new_fpr, "numeric")
  )
  roc$m_fpr <- new_fpr
  return(roc)
}

get_fpr <- function(roc) {
  UseMethod("get_fpr", roc)
}

get_fpr.roc <- function(roc) {
  # Return the m_fpr member of a roc object
  # roc: a roc object
  return(roc$m_fpr)
}

set_tpr <- function(roc, new_tpr) {
  UseMethod("set_tpr", roc)
}

set_tpr.roc <- function(roc, new_tpr) {
  # Set the m_tpr member of a roc object.
  # roc: a roc object
  # new_fpr: a numeric vector
  stopifnot(
    length(new_tpr) == length(roc$m_tpr)
    && is.vector(new_tpr, "numeric")
  )
  roc$m_tpr <- new_tpr
  return(roc)
}

get_tpr <- function(roc) {
  UseMethod("get_tpr", roc)
}

get_tpr.roc <- function(roc) {
  # Return the m_fpr member of a roc object
  # roc: a roc object
  return(roc$m_tpr)
}

compute_auc <- function(roc) {
  UseMethod("compute_auc", roc)
}

trapezoid_area <- function(x1, x2, y1, y2) {
  # Computes the area of the trapezoid specified by x1, x2, y1, y2.
  # 
  # Parameters
  # ----------
  # predictions: numeric vector
  #    Predictions stored as a numeric vector
  # Returns
  # -------
  # The area of the trapezoid.  
  base <- abs(x1 - x2)
  height_avg <- (y1 + y2) / 2
  return(base * height_avg)
}

compute_auc <- function(roc) {
  UseMethod("compute_auc", roc)
}

compute_auc.roc <- function(roc) {
  # Computes the area under the receiver operating characteristic curve.
  # Algorithm
  # ---------
  # See "An introduction to ROC analysis, Tom Fawcett, PRL, 2006"
  # Parameters
  # ----------
  # roc: numeric vector
  #    A roc object
  # Returns
  # -------
  # The area under the receiver operating characteristic curve.
  auc <- 0
  # summing the areas of trapezoids
  i <- 1
  while (i <= (length(roc$m_fpr) - 1)) {
    fpr_cur <- roc$m_fpr[i]
    fpr_next <- roc$m_fpr[i + 1]
    tpr_cur <- roc$m_tpr[i]
    tpr_next <- roc$m_tpr[i + 1]
    while (TRUE) {
      i <- i + 1
      if (fpr_next == fpr_cur) {
        fpr_cur <- fpr_next
        fpr_next <- roc$m_fpr[i + 1]
        tpr_cur <- tpr_next
        tpr_next <- roc$m_tpr[i + 1]
      } else {
        break
      }
    }
    auc <- auc + trapezoid_area(fpr_cur, fpr_next, tpr_cur, tpr_next)
  }
  return(auc)
}

interpolate <- function(rocp1, rocp2, fpr) {
  # Interpolate the true positive rate from two roc points.
  # Parameters
  # ----------
  # rocp1: numeric
  #    ROC point 1
  # rocp2: numeric
  #    ROC point 2
  # fpr: numeric
  #    The false positive rate for which the true positive rate is to be interpolated
  # Returns
  # -------
  # The true positive rate.
  slope <- (rocp2[2] - rocp1[2]) / (rocp2[1] - rocp1[1])
  return(rocp1[2] + slope * (fpr - rocp1[1]))
}

tpr_for_fpr <- function(fpr, roc) {
  # Compute the true positive rate given a false positive rate.
  # Algorithm
  # ----------
  # See "An introduction to ROC analysis, Tom Fawcett, PRL, 2006"
  # Parameters
  # ----------
  # fpr: numeric
  #    False positive rate
  # roc: numeric
  #    A roc object
  # Returns
  # -------
  # The true positive rate. 
  i <- 1
  npts <- length(roc$m_fpr) # number of roc points
  while (i < npts && roc$m_fpr[i] <= fpr) {
    i <- i + 1
  }
  if (roc$m_fpr[i - 1] == fpr) {
    return(roc$m_tpr[i - 1])
  } else {
    rocp1 <- c(roc$m_fpr[i - 1], roc$m_tpr[i - 1])
    rocp2 <- c(roc$m_fpr[i], roc$m_tpr[i])
    return(interpolate(rocp1, rocp2, fpr))
  }
}

average_roc.roc <- function(rocs, fp_samples) {
  # Averages the receiver operating characteristic curve.
  # Algorithm
  # ----------
  # See "An introduction to ROC analysis, Tom Fawcett, PRL, 2006"
  # Parameters
  # ----------
  # rocs: roc
  #    A vector of roc objects
  # fp_samples:
  #    The number of false positive samples
  # Returns
  # -------
  # The averaged roc object.
  s <- 1
  avg_tprs <- vector("numeric")
  for (fpr in seq(0, 1, 1 / fp_samples)) {
    tpr_sum <- 0
    nrocs <- length(rocs)
    for (i in 1:nrocs) {
      tpr_sum <- tpr_sum + tpr_for_fpr(fpr, rocs[i])
    }
    avg_tprs[s] <- tpr_sum / nrocs
    s <- s + 1
  }
  avg_roc <- roc(seq(0, 1, 1 / fp_samples), avg_tprs)
  return(avg_roc)
}