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

AnalyseSplittedExpt <- function(gt.alive, gt.dead, min.freq.al = NULL,
                                location.cols = TRUE, p.values = TRUE, 
                                genotypic = FALSE, backup.path = NULL){

  # Calculate frequencies -----------------------------------------------------
  freq.alive <- CalcFreqGt(gt.alive, genotypic = TRUE, allelic = TRUE,
                            absolute = c(TRUE, FALSE))
  freq.dead  <- CalcFreqGt(gt.dead, genotypic = TRUE, allelic = TRUE,
                            absolute = c(TRUE, FALSE))
  freq.alive <- data.frame(freq.alive)
  freq.dead  <- data.frame(freq.dead)
  # Only keep rows in common --------------------------------------------------

  filter.alive <- rownames(freq.alive) %in% rownames(freq.dead)
  if (sum(!filter.alive) > 0) {
    freq.alive <- freq.alive[filter.alive, ]
  }
  filter.dead <- rownames(freq.dead) %in% rownames(freq.alive)
  if (sum(!filter.dead) > 0) {
    freq.dead <- freq.dead[filter.dead, ]
  }

  # Create frequency matrix ---------------------------------------------------

  first.freq.alive <- SliceDfRows(freq.alive, 1) # Freqs for the 1st variant
  map.gt.alive <- FindIdsGtCounts(first.freq.alive) # Ids of intersting freqs
  map.al.alive <- FindIdsAlFreqs(first.freq.alive)
  # how many samples survived (ratio)
  ratio.alive  <- ncol(gt.alive) / (ncol(gt.alive) + ncol(gt.dead))

  freq.all   <- cbind(SliceDfColumns(freq.alive, unlist(map.gt.alive)) 
                      + SliceDfColumns(freq.dead, unlist(map.gt.alive)),
                      SliceDfColumns(freq.alive, unlist(map.al.alive))
                      * ratio.alive 
                      + SliceDfColumns(freq.dead, unlist(map.al.alive)) 
                      * (1 - ratio.alive),
                      SliceDfColumns(freq.alive, ncol(freq.alive)) # <-- bad things here
                      + SliceDfColumns(freq.alive, ncol(freq.alive)))

  first.freq.all <- SliceDfRows(freq.all, 1)
  map.gt.all     <- FindIdsGtCounts(first.freq.all)
  map.al.all     <- FindIdsAlFreqs(first.freq.all)

  freqs <- cbind(freq.alive, freq.all)

  # modify frequency matrix column names
  colnames.alive  <- paste("SURVIVERS", colnames(freq.alive), sep = ".")
  colnames.all    <- paste("ALL", colnames(freq.alive), sep = ".")
  colnames(freqs) <- c(colnames.alive, colnames.all)

  # Remove variants with low allelic frequencies ------------------------------

  if (is.numeric(min.freq.al)) {

      freq.al.ids <- c(map.al.all$ref, map.al.all$alt)
      filter <- !as.logical(rowSums(freq.all[, freq.al.ids] < min.freq.al,
                                    na.rm = TRUE))
      freqs <- freqs[filter, ]
      freq.alive <- freq.alive[filter, ]
      freq.all <- freq.all [filter, ]

      if (is.null(dim(freqs))) {
        stop("You removed all rows from gt using parameter min.freq.al")
      }

  }

  # Calculate probabilities ---------------------------------------------------

  if (p.values) {

    col.all  <- ncol(freq.alive) + 1
    col.end  <- ncol(freq.alive) + ncol(freq.all)
    probs <- apply(cbind(freq.alive, freq.all), MARGIN = 1,
                      FUN = function(x) {
                        CalcProbsSelection (x[1:(col.all - 1)],
                                            x[col.all:col.end],
                                            map.alive = map.gt.alive,
                                            map.all = map.gt.all)
                      }
                )
    freqs <- cbind(freqs, t(probs))

  }

  # Remove gt stats columns ---------------------------------------------------

  if (!genotypic) {

    filter <- !grepl("count\\.gt\\.", colnames(freqs))
    freqs <- freqs[, filter]

  }

  # Add location columns ------------------------------------------------------

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

  # Save to a file ------------------------------------------------------------

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

AnalyseExpt <- function(gt, survival, location.cols = TRUE, min.freq.al = NULL,
                        p.values = TRUE, genotypic = FALSE,
                        backup.path = NULL){

  splitted.gt <- SplitGt(gt, survival)
  gt.alive    <- splitted.gt$alive
  gt.dead     <- splitted.gt$dead

  result <- AnalyseSplittedExpt(gt.alive, gt.dead,
                                location.cols = location.cols, 
                                min.freq.al = min.freq.al,
                                p.values = p.values, genotypic = genotypic,
                                backup.path = backup.path)
  return(result)

}


#' A function to calculate probabilities that observed frequencies are due to
#' positive selection using MWNCHHypergeo model.
#'
#' @param freq.alive An vector of observed absolute frequencies for the three
#' genotypes in the alive population
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
#' CalcProbsSelection(f.alive, f.all)

CalcProbsSelection <- function(freq.alive, freq.all, map.alive = NULL,
                               map.all = NULL){

  sel.weights <- CalcWeightsSurvival(freq.alive, freq.all)

  if ((sum(is.na(sel.weights)) == 0) && (sum(sel.weights <= 0) == 0) 
      && (sum(is.infinite(sel.weights)) == 0)) {

    if(!is.null(map.alive)) {
      id.alive <- map.alive
    } else {
      id.alive <- FindIdsGtCounts(freq.alive)
    }

    if(!is.null(map.all)) {
      id.all <- map.all
    } else {
      id.all <- FindIdsGtCounts(freq.all)
    }

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

CalcWeightsSurvival <- function(freq.alive, freq.all, map.alive = NULL,
                                map.all = NULL){

    if(!is.null(map.alive)) {
      id.alive <- map.alive
    } else {
      id.alive <- FindIdsGtCounts(freq.alive)
    }

    if(!is.null(map.all)) {
      id.all <- map.all
    } else {
      id.all <- FindIdsGtCounts(freq.all)
    }

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