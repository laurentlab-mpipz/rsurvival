#' A function to analyse a selection experiment, including observed frequencies
#' of genotype and probabilities that they are due to random picking and
#' positive selection.
#'
#' @param gt.alive A genotype data frame of the alive population
#' @param gt.dead A genotype data frame of the dead population
#' @param location.cols If TRUE, adds two columns to the result, which are the
#' scaffold name and the locus name as factors (deafult is FALSE)
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
#' AnalyseSplittedExpt(gt.alive, gt.dead)
#' AnalyseSplittedExpt(gt.alive, gt.dead, p.values = FALSE)
#' AnalyseSplittedExpt(gt.alive, gt.dead, backup.path = "example.csv")

AnalyseSplittedExpt <- function(gt.alive, gt.dead, location.cols = TRUE,
                        p.values = TRUE, genotypic = FALSE,
                        backup.path = NULL){

  # calculate frequencies
  freq.alive <- CalcFreqGt(gt.alive, genotypic = TRUE, allelic = TRUE,
                            absolute = TRUE)
  freq.dead  <- CalcFreqGt(gt.dead, genotypic = TRUE, allelic = TRUE,
                            absolute = TRUE)
  dim(freq.alive)
  # only keep rows in common
  filter.alive <- rownames(freq.alive) %in% rownames(freq.dead)
  if (sum(!filter.alive) > 0) {
    freq.alive <- freq.alive[filter.alive, ]
  }
  filter.dead <- rownames(freq.dead) %in% rownames(freq.alive)
  if (sum(!filter.dead) > 0) {
    freq.dead <- freq.dead[filter.dead, ]
  }

  freq.all   <- freq.alive + freq.dead
  freqs      <- cbind(freq.alive, freq.all)

  # modify frequencies column names
  colnames.alive  <- paste("SURVIVERS", colnames(freq.alive), sep = ".")
  colnames.all    <- paste("ALL", colnames(freq.alive), sep = ".")
  colnames(freqs) <- c(colnames.alive, colnames.all)

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

  if (!genotypic) {

    filter <- !grepl("count\\.gt\\.", colnames(freqs))
    freqs <- freqs[, filter]

  }

  if (location.cols) {
    temp.names <- colnames(freqs)
    lim <- stringr::str_locate(rownames(freqs), "_")[, 1]
    freqs <- cbind(data.frame(stringr::str_sub(rownames(freqs), lim + 1, -1)),
                    freqs)
    freqs <- cbind(data.frame(stringr::str_sub(rownames(freqs), 1, lim - 1)),
                    freqs)
    colnames(freqs) <- c("CHROM", "LOCUS", temp.names)
  }

  result <- freqs

  if (is.character(backup.path)) {
    utils::write.csv(result, backup.path) # write to file
  }

  return(result)

}


#' A function to analyse a selection experiment, including observed frequencies
#' of genotype and probabilities that they are due to random picking and
#' positive selection.
#'
#' @param gt A genotype data frame of your original population
#' @param survival A logical vector. TRUE survived the experiment, FALSE died
#' @param location.cols If TRUE, adds two columns to the result, which are the
#' scaffold name and the locus name as factors (deafult is FALSE)
#' @param p.values If TRUE (default), some probabilities will be included in
#' the result, as well as frequencies
#' @param backup.path Optionnal. A path where backup files can be stored
#'
#' @return 
#' A data frame of frequencies and probabilities from your experiment.
#'
#' @seealso For more information, see \code{\link{AnalyseSplittedExpt}} which
#' this function bind. 
#'
#' @export
#'
#' @examples
#' AnalyseExpt(gt, survival)
#' AnalyseExpt(gt, survival, p.values = FALSE)
#' AnalyseExpt(gt, survival, backup.path = "example.csv")

AnalyseExpt <- function(gt, survival, location.cols = TRUE, 
                        p.values = TRUE, genotypic = FALSE,
                        backup.path = NULL){

  splitted.gt <- SplitGt(gt, survival)
  gt.alive    <- splitted.gt$alive
  gt.dead     <- splitted.gt$dead

  result <- AnalyseSplittedExpt(gt.alive, gt.dead,
                                location.cols = location.cols,
                                p.values = p.values, genotypic = genotypic,
                                backup.path = backup.path)
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
#' @seealso For more information, see \code{\link{CalcWeightsSurvival}} wich
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

  sel.weights <- CalcWeightsSurvival(freq.alive, freq.all)

  if ((sum(is.na(sel.weights)) == 0) && (sum(sel.weights <= 0) == 0) 
      && (sum(is.infinite(sel.weights)) == 0)) {

    id.alive <- FindIdsGtCounts(freq.alive)
    id.all   <- FindIdsGtCounts(freq.all)

    # distribution in alive population
    x <- freq.alive[c(id.alive$homo.ref, id.alive$hetero, id.alive$homo.alt)]
    # distribution in aliveoriginal population
    m <- freq.all[c(id.all$homo.ref, id.all$hetero, id.all$homo.alt)]
    # nb of alive samples with valid data
    n <- freq.alive['count.gt.TOTAL'] - freq.alive['count.gt.MISSVAL']

    prob.sel     <- BiasedUrn::dMWNCHypergeo(x = x, m = m, n = n,
                                          odds = sel.weights)
    prob.neutral <- BiasedUrn::dMWNCHypergeo(x = x, m = m, n = n,
                                              odds = c(1,1,1))

    lrt     <- 2 * log(prob.sel / prob.neutral)
    p.value <- 1 - stats::pchisq(q = lrt, df = 2)  
    result  <- c(sel.weights, prob.neutral, prob.sel, lrt, p.value)

  } else {
    result <- rep(NA, 7) # to ensure that the result always has the same length
  }

  names(result) <- c("weight.gt.HOMOREF", "weight.gt.HETERO",
                      "weight.gt.HOMOALT", "p.neutral", "p.select",  "lrt",
                      "p.value")
  return(result)

}


#' A function to calculate odds for positive selection from observed selection
#' frequencies using MWNCHHypergeo model. 
#'
#' @param freq.alive A vector of observed absolute frequencies for the three
#' genotypes in the alive population
#' @param freq.all An vector of observed absolute frequencies for the three
#' genotypes in the population before selection
#'
#' @return 
#' A vector of odds for the 3 genotypes, according to the MWNCHHypergeo model.
#'
#' @export
#'
#' @examples
#' f.alive <- c(17, 3, 15, 1, 36)
#' f.all  <- c(24, 5, 22, 2, 53)
#' CalcWeightsSurvival(f.alive, f.all)

CalcWeightsSurvival <- function(freq.alive, freq.all){
    
    id.alive <- FindIdsGtCounts(freq.alive)
    id.all   <- FindIdsGtCounts(freq.all)

    # distribution of alive samples
    mu <- freq.alive[c(id.alive$homo.ref, id.alive$hetero, id.alive$homo.alt)]
    # distribution of original population
    m  <- freq.all[c(id.all$homo.ref, id.all$hetero, id.all$homo.alt)] 
    # number of alive samples with valid data
    n  <- freq.alive[id.alive$total] - freq.alive[id.alive$missval]
    


    # calculate odds
    result <- suppressWarnings(BiasedUrn::oddsMWNCHypergeo(mu = mu,
                                                            m = m, n = n))
    return(result)

}