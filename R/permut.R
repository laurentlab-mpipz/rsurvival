#' A function to calculate probabilities of neutral model and positive
#' selection model from a genotype matrix whit a randomly shuffled survival
#' vector
#'
#' @param gt A genotype matrix
#' @param survival A survival vector to shuffle
#' @param odded.lives Optional. A vector containing a number of lives for each
#' categories of \code{odded.samples} (e.g. c("HETERO" = 32, "HOMOREF" = 18, 
#' ... ))
#' @param odded.samples Optional. A list of vectors containing index for each
#' odded sample categories. (e.g. list$HETERO could be c(1,3,7,45,141))
#' @param verbose Logical. If TRUE (default), report status of the process
#' along the execution.
#'
#' @return 
#' A data frame of frequencies and probabilities, one snp per row.
#'
#' @export
#'
#' @seealso For more detail, see \code{link{AnalyseSplittedExpt}} which this function
#' bind.
#'
#' @examples
#' PermutSurv(gt, c(TRUE, TRUE, FALSE, ...))
#' PermutSurv(gt, survival, verbose = FALSE)
#' PermutSurv(gt, survival, odded.lives = c("HETERO" = 32, "HOMOREF" = 18, ...
#' ), odded.samples = list("HETERO" = c(1,3,7,45,141), ...)))

PermutRandSurv <- function(gt, survival, odded.lives = NULL,
                           odded.samples = NULL, verbose = TRUE){

  result <- NULL

  if (!is.null(odded.lives) && !is.null(odded.samples)) {

    survival <- rep(FALSE, length(survival))
    survival[sample(odded.samples$HOMOREF)[1:odded.lives['HOMOREF']]] <- TRUE
    survival[sample(odded.samples$HETERO)[1:odded.lives['HETERO']]]   <- TRUE
    survival[sample(odded.samples$HOMOALT)[1:odded.lives['HOMOALT']]] <- TRUE

  } else {

    survival <- sample(survival) # pseudo-random shuffle

  } 

  probs <- AnalyseExpt(gt, survival, min.freq.al = 0.1)
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
#' @seealso For more detail, see \code{link{AnalyseSplittedExpt}} or
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
#' @param odded.pos Optional. Reference to a variant in \code{gt} as a string.
#' Iterate random pick according to the observed odds of genotype at this
#' position in the genome.
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

IterRandPick <- function(gt, survival, max.p.neutral = 0.1, iter = 10,
                         odded.pos = NULL, verbose = TRUE){

  result <- c()

  # Prepare number of deaths if needed ----------------------------------------

  if (is.character(odded.pos)) {

    odded.variant <- gt[odded.pos, ]
    odded.pick    <- CalcOddedPick(odded.variant, survival, n = iter,
                                   verbose = verbose)
    odded.lives   <- odded.pick$odded.lives
    odded.samples <- odded.pick$odded.samples

  }

  # Start iterating -----------------------------------------------------------

  time.stamp <- proc.time()[3] # ellapsed time since session launched
  time.left  <- NA # estimated remaining time

  for (i in 1:iter) {

    # print the iteration number and the estimated remaining time
    if (verbose) {
      message.start = paste("Iterating : ", i, " out of ", iter, sep = "")
      if (i == 1) {
        message(paste(message.start , " (Estimating time)", sep = ""))
      } else {
        message(paste(message.start, " (", time.left, " remaining)" ,sep = ""))
      }
    }

    # calculate probs with a random survival vector
    if (is.character(odded.pos)){
      probs.rand <- PermutRandSurv(gt, survival,
                                   odded.lives = odded.lives[, i],
                                   odded.samples = odded.samples,
                                   verbose = FALSE)
    } else {
      probs.rand <- PermutRandSurv(gt, survival, verbose = FALSE)
    }

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


#' A function to calculate the number of picked lives per categories of samples
#' according to MWNCHypergeo model. 
#'
#' @param variant The odded variant vector
#' @param survival A survival vector
#' @param n Number of life vectors to be generated (default is 10)
#' @param verbose Logical. If TRUE (default), report status of the process
#' along the execution.
#'
#' @return 
#' A list. list$odded.lives is a matrix containing a number of survivors per 
#' categories (e.g. list$odded.lives row names could be "HOMOREF", "HOMOALT",
#' ...) with \code{n} column, each time generated with the same input values.
#' list$odded.samples are the index of the samples belonging to each categories
#' (e.g. list$odded.samples$HETERO could be c(1, 3, 8, 14))
#' @export
#'
#' @seealso For more detail, see \code{link{IterRandPick}} which binds this
#' function
#'
#' @examples
#' CalcOddedPick(variant, survival)
#' CalcOddedPick(variant, survival, n = 1000)
#' CalcOddedPick(variant, survival, verbose = FALSE)

CalcOddedPick <- function(variant, survival, n = 10, verbose = TRUE){

  if (verbose) {
      message(paste("Preparing : ", n * sum(survival),
                    " odded survivals, it can take a while...", sep = ""))
    }

    # split variant in which positive selection will be applied
    split.variant <- SplitVect(variant, survival)
    variant.alive <- t(data.frame("variant" = split.variant$alive,
                                  stringsAsFactors = FALSE))
    variant.dead  <- t(data.frame("variant" = split.variant$dead,
                                  stringsAsFactors = FALSE))
    rownames(variant.alive) <- c("odded.variant")
    rownames(variant.dead)  <- c("odded.variant")

    # ids of samples concerned by each genotype
    gt.ids <- list("HOMOREF" = grep("0.*0", variant),
                   "HETERO" = grep("0.*1|1.*0", variant),
                   "HOMOALT" = grep("1.*1", variant))

    # stats for observed selection during the real experiment
    probs <- AnalyseSplittedExpt(variant.alive, variant.dead, genotypic = TRUE)
    freqs.all <- probs[, c("ALL.count.gt.HOMOREF", "ALL.count.gt.HETERO",
                           "ALL.count.gt.HOMOALT")]
    odds <- probs[, c("weight.gt.HOMOREF", "weight.gt.HETERO",
                      "weight.gt.HOMOALT")]

    # simulate number of deaths per genotypes with the observed odds
    nb.lives <- BiasedUrn::rMWNCHypergeo(nran = n, m = freqs.all,
                                         n = sum(survival), odds = odds)
    rownames(nb.lives) <- c("HOMOREF", "HETERO", "HOMOALT")
    colnames(nb.lives) <- 1:n

    result <- list("odded.lives" = nb.lives, "odded.samples" = gt.ids)
    
    if (verbose) {
      message("Preparing : DONE.")
    }

    return (result)

}