#' A function to analyse a selection experiment, including observed frequencies
#' of genotype and probabilities that they are due to random picking and
#' positive selection.
#'
#' @param gt.alive A genotype data frame of the alive population
#' @param gt.dead A genotype data frame of the dead population
#' @param p.values If TRUE (default), some probabilities will be included in
#' the result, as well as frequencies
#' @param backup.path Optionnal. A path where backup files can be stored
#'
#' @return 
#' A data frame of frequencies and probabilities from your experiment.
#'
#' @seealso For more information, see \code{\link{CalcProbsSelection}} and
#' \code{\link{CalcFreqGt}} which this function bind. 
#'
#' @export
#'
#' @examples
#' AnalyseExpt(gt.alive, gt.dead)
#' AnalyseExpt(gt.alive, gt.dead, p.values = FALSE)
#' AnalyseExpt(gt.alive, gt.dead, backup.path = "example.csv")

AnalyseExpt <- function(gt.alive, gt.dead, p.values = TRUE,
                              backup.path = NULL){

  # calculate frequencies
  freq.alive <- CalcFreqGt(gt.alive, genotypic = TRUE, allelic = TRUE,
                        absolute = TRUE)
  freq.dead  <- CalcFreqGt(gt.dead, genotypic = TRUE, allelic = TRUE,
                        absolute = TRUE)

  # only keep rows in common
  freq.alive <- freq.alive[rownames(freq.alive) %in% rownames(freq.dead), ]
  freq.dead  <- freq.dead[rownames(freq.dead) %in% rownames(freq.alive), ]
  freq.all   <- freq.alive + freq.dead
  freqs      <- cbind(freq.alive, freq.all)

  # modify frequencies column names
  colnames.surv   <- paste("SURVIVERS", colnames(freq.alive), sep = "_")
  colnames.all    <- paste("ALL", colnames(freq.alive), sep = "_")
  colnames(freqs) <- c(colnames.surv, colnames.all)

  if (p.values) {

    col.dead <- ncol(freq.alive) + 1
    col.all  <- ncol(freq.alive) + ncol(freq.dead) + 1
    col.end  <- ncol(freq.alive) + ncol(freq.dead) + ncol(freq.all)

    probs <- apply(cbind(freq.alive, freq.dead, freq.all), MARGIN = 1,
                      FUN = function(x) {
                        CalcProbsSelection (x[1:(col.dead - 1)],
                                            x[col.dead:(col.all - 1)], 
                                            x[col.all:col.end])
                      }
                )

    freqs <- cbind(freqs, t(probs))

  }

  result <- freqs

  if (is.character(backup.path)) {
    utils::write.csv(result, backup.path) # write to file
  }

  return(result)

}


#' A function to calculate probabilities that observed frequencies are due to
#' positive selection using MWNCHHypergeo model.
#'
#' @param freq.alive An vector of observed absolute frequencies for the three
#' genotypes in the alive population
#' @param freq.dead An vector of observed absolute frequencies for the three
#' genotypes in the dead population
#' @param freq.all An vector of observed absolute frequencies for the three
#' genotypes in the population before selection
#'
#' @return 
#' A vector of probabilities : \code{p_neutral} is the probability to obtain
#' this distribution with a random picking process, \code{p_select} is the
#' probability to obtain this distribution with a picking in a biased urn,
#' \code{weights} are the odds used to calculate \code{p_select}, \code{lrt}
#' is the result of the likelyhood test of \code{p_select} uppon
#' \code{p_neutral}, \code{p_value} is the probability that the distribution
#' is due to positive selection.
#'
#' @seealso For more information, see \code{\link{CalcWeightsSelection}} wich
#' this function bind.
#'
#' @export
#'
#' @examples
#' f.alive <- c(17, 3, 15, 1, 36)
#' f.all   <- c(24, 5, 22, 2, 53)
#' f.dead  <- f.all - f.dead
#' CalcProbsSelection(f.alive, f.dead, f.all)

CalcProbsSelection <- function(freq.alive, freq.dead, freq.all){

  sel.weights <- CalcWeightsSelection(freq.dead, freq.all)

  if ((sum(is.na(sel.weights)) == 0) && (sum(sel.weights <= 0) == 0) 
      && (sum(is.infinite(sel.weights)) == 0)) {

    x <- freq.alive[1:3] # distribution in alive population
    m <- freq.all[1:3] # distribution in aliveoriginal population
    n <- freq.alive[5] - freq.alive[4] # nb of alive samples with valid data

    prob.sel <- BiasedUrn::dMWNCHypergeo(x = x, m = m, n = n,
                                          odds = (1 / sel.weights))
    prob.neutral <- BiasedUrn::dMWNCHypergeo(x = x, m = m, n = n,
                                              odds = c(1,1,1))

    lrt <- 2 * log(prob.sel / prob.neutral)
    p.value <- 1 - stats::pchisq(q = lrt, df = 2)

    result <- c(round(prob.neutral, 4), round(prob.sel, 4),
                toString(round(sel.weights, 4)), round(lrt,4),
                round(p.value,4))

  } else {
    result <- rep(NA, 5)

  }

  names(result) <- c("p.neutral", "p.select", "weigths",  "lrt", "p.value")
  return(result)

}


#' A function to calculate odds for positive selection from observed selection
#' frequencies using MWNCHHypergeo model. 
#'
#' @param freq.dead An vector of observed absolute frequencies for the three
#' genotypes in the dead population
#' @param freq.all An vector of observed absolute frequencies for the three
#' genotypes in the population before selection
#'
#' @return 
#' A vector of odds for the 3 genotypes, according to the MWNCHHypergeo model.
#'
#' @export
#'
#' @examples
#' f.dead <- c(17, 3, 15, 1, 36)
#' f.all  <- c(24, 5, 22, 2, 53)
#' CalcWeightsSelection(f.dead, f.all)

CalcWeightsSelection <- function(freq.dead, freq.all){
    
    mu <- freq.dead[1:3] # distribution of dead samples
    m  <- freq.all[1:3] # distribution of original population
    n  <- freq.dead[5] - freq.dead[4] # number of dead samples with valid data

    # calculate odds
    result <- suppressWarnings(BiasedUrn::oddsMWNCHypergeo(mu = mu,
                                                            m = m, n = n))

    return(result)

}