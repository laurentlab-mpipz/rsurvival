#' A function to calculate probabilities of neutral model and positive
#' selection model from a genotype matrix whit a randomly shuffled survival
#' vector
#'
#' @param gt A genotype matrix
#' @param survival A survival vector to shuffle
#' @param odds Optional. Not used yet.
#' @param verbose Logical. If TRUE (default), report status of the process
#' along the execution.
#'
#' @return 
#' A data frame of frequencies and probabilities, one snp per row.
#'
#' @export
#'
#' @seealso For more detail, see \code{link{AnalyseExpt}} which this function
#' bind.
#'
#' @examples
#' PermutSurv(gt, c(TRUE, TRUE, FALSE, ...))
#' PermutSurv(gt, survival, verbose = FALSE)

PermutSurv <- function(gt, survival, odds = NULL, verbose = TRUE){

  survival 	  <- sample(survival) # pseudo-random shuffle
  splitted.gt <- SplitGt(gt$gt, survival, verbose = verbose)
  probs       <- AnalyseExpt(splitted.gt$alive, splitted.gt$dead)

  result <- probs

  return(result)

}


#' A function to count snps with a low probability of being due to a neutral
#' model
#'
#' @param probs A probability data frame
#' @param max.p.neutral A threshold for \code{p.neutral} under which the snp
#' will be counted (default is 0.01)
#'
#' @return 
#' An integer. The number of snps with low \code{p.neutral} value in
#' \code{probs} data frame.
#'
#' @export
#'
#' @seealso For more detail, see \code{link{AnalyseExpt}} or
#' \code{link{CalcProbsSelection}} which can provide you a probs data frame.
#'
#' @examples
#' CountSignSnps(probs)
#' CountSignSnps(probs, max.p.neutral = 0.05)

CountSignSnps <- function(probs, max.p.neutral = 0.01){

  # bukd a list of p.neutral from probs data frame, exlcuding NAs
  p.neutral.list <- probs[, 'p.neutral'][!is.na(probs[, 'p.neutral'])]
  # count significant p.neutral (i.e. under max.p.neutral threshold)
  result <- sum(as.numeric(p.neutral.list) <= max.p.neutral)
	
  return(result)

}


#' A function to iterate significant snps counting in a genotype matrx with a
#' random permutation of the survival vector
#'
#' @param gt A genotype matrix
#' @param survival A survival vector to shuffle
#' @param max.p.neutral A threshold for \code{p.neutral} under which the snp
#' will be counted (default is 0.01)
#' @param iter Number of iteration (default is 10)
#' @param odds Optional. Not used yet.
#' @param verbose Logical. If TRUE (default), report status of the process
#' along the execution.
#'
#' @return 
#' An integer vector. Each item is the number of snps with low \code{p.neutral}
#' value in \code{gt} data frame after a random permutation of the survival
#' vector. Countains \code{iter} items.
#'
#' @export
#'
#' @seealso For more detail, see \code{link{PermutSurv}} and
#' \code{link{CountSignSnps}} which this function bind.
#'
#' @examples
#' IterRandomPick(gt, survival)
#' IterRandomPick(gt, survival, max.p.neutral = 0.05)
#' IterRandomPick(gt, survival, iter = 1000, verbose = FALSE)

IterRandomPick <- function(gt, survival, max.p.neutral = 0.1, iter = 10,
                            odds = NULL, verbose = TRUE){

  result <- c()

  time.stamp <- proc.time()[3] # ellapsed time since session launched
  time.left  <- NA # estimated remaining time

  for (i in 1:iter) {

    if (verbose) {
      # print the iteration number and the estimated remaining time
      if (i == 1) {
        message(paste("Iterating : ", i, " out of ", iter, " (Estimating time)",
                      sep = ""))
      } else {
        message(paste("Iterating : ", i, " out of ", iter, " (", time.left,
                      " remaining)" ,sep = ""))
      }
    }

    # calculate probs with a random survival vector
    probs.rand <- PermutSurv(gt, survival, verbose = FALSE)
    # count significant p.neutral and append the count to result
    result     <- c(result, CountSignSnps(probs.rand, max.p.neutral))

    # estimating remaining time from the last iteration time
    time.left  <- round((proc.time()[3] - time.stamp) * (iter - i))
    time.left  <- lubridate::seconds_to_period(time.left) # format to time
    time.stamp <- proc.time()[3]

  }

  if (verbose) {
    message("Iterating : DONE.")
  }
  
  return(result)

}