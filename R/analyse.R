# -----------------------------------------------------------------------------
# Copyright © 2018 Samuel Badion, Stefan Laurent. 
#
# This file is part of Rsurvival.
#
# Rsurvival is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or any later version.
#
# Rsurvival is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# Rsurvival. If not, see <https://www.gnu.org/licenses/>
# -----------------------------------------------------------------------------


#' A function to analyse a selection experiment, including observed frequencies
#' of genotype and probabilities that they are due to random picking and
#' positive selection.
#'
#' @param gt.alive A genotype data frame of the alive population
#' @param gt.dead A genotype data frame of the dead population
#' @param min.freq.al Optional. Numeric from 0 to 1. Variants with MAF under
#' this threshold will be ommited.
#' @param location.cols If TRUE, adds two columns to the result, which are the
#' scaffold name and the locus name as factors (deafult is FALSE)
#' @param deltas If TRUE (default), include columns of differencies between
#' relative allelic frequencies after selection and before selection in the
#' result. Positive value means that the allele is more common in the survivers
#' population.
#' @param p.values If TRUE (default), some probabilities will be included in
#' the result, as well as frequencies
#' @param genotypic If TRUE, will include counts of genotype in the result.
#' @param gt.pop1 Optional. A genotype matrix. If both \code{gt.pop1} and \code{gt.pop2}
#' are set, will include the W&C FST estimation to the result.
#' @param gt.pop2 Optional. A genotype matrix. If both \code{gt.pop1} and \code{gt.pop2}
#' are set, will include the W&C FST estimation to the result.
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
#' \dontrun{
#' AnalyseSplittedExpt(gt.alive, gt.dead)
#' AnalyseSplittedExpt(gt.alive, gt.dead, min.freq.al = 0.1)
#' AnalyseSplittedExpt(gt.alive, gt.dead, location.cols = FALSE)
#' AnalyseSplittedExpt(gt.alive, gt.dead, deltas = FALSE, p.values = FALSE)
#' AnalyseSplittedExpt(gt.alive, gt.dead, backup.path = "example.csv")
#' }

AnalyseSplittedExpt <- function(gt.alive, gt.dead, min.freq.al = NULL,
                                location.cols = TRUE, deltas = TRUE,
                                p.values = TRUE, genotypic = FALSE,
                                gt.pop1 = NULL, gt.pop2 = NULL,
                                backup.path = NULL){

  # Calculate frequencies -----------------------------------------------------
  freq.alive <- CalcFreqGt(gt.alive, genotypic = TRUE, allelic = TRUE,
                            absolute = c(TRUE, FALSE))
  freq.all   <- CalcFreqGt(cbind(gt.alive, gt.dead), genotypic = TRUE,
                          allelic = TRUE, absolute = c(TRUE, FALSE))
  freq.alive <- data.frame(freq.alive)
  freq.all   <- data.frame(freq.all)

  # Only keep rows in common --------------------------------------------------
  filter.alive <- rownames(freq.alive) %in% rownames(freq.all)
  if (sum(!filter.alive) > 0) {
    freq.alive <- freq.alive[filter.alive, ]
  }
  filter.all <- rownames(freq.all) %in% rownames(freq.alive)
  if (sum(!filter.all) > 0) {
    freq.all <- freq.all[filter.all, ]
  }

  # Create frequency matrix ---------------------------------------------------

  first.freq.alive <- SliceDfRows(freq.alive, 1) # Freqs for the 1st variant
  map.gt.alive <- FindIdsGtCounts(first.freq.alive) # Ids of intersting freqs
  map.al.alive <- FindIdsAlFreqs(first.freq.alive)

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
      freqs      <- freqs[filter, ]
      freq.alive <- freq.alive[filter, ]
      freq.all   <- freq.all[filter, ]

      if (is.null(dim(freqs))) {
        stop("You removed all rows from gt using parameter min.freq.al")
      }

  }

  # Calculate deltas ----------------------------------------------------------

  if (deltas) {

    deltas <- (SliceDfColumns(freq.alive, c(map.al.alive$ref, map.al.alive$alt))
               - SliceDfColumns(freq.all, c(map.al.all$ref, map.al.all$alt)))
    colnames(deltas) <- c("DELTA.freq.al.REF", "DELTA.freq.al.ALT")
    freqs <- cbind(freqs, deltas)

  }

  # Estimate FST --------------------------------------------------------------

  if (!is.null(gt.pop1) && !is.null(gt.pop2)) {

    delimiter <- ncol(gt.pop1)
    end <- delimiter + ncol(gt.pop2) 
    fst <- apply(cbind(gt.pop1, gt.pop2), MARGIN = 1,
                 FUN = function(x){
                  x <- data.frame(t(x))
                  variant.pop1 <- SliceDfColumns(x, 1:delimiter)
                  variant.pop2 <- SliceDfColumns(x, delimiter:end)
                  result <- EstimateFST(variant.pop1, variant.pop2)
                  return(result)
                 })
    
    fst <- data.frame(fst)
    colnames(fst) <- c("FST.pop1.pop2")
    freqs <- cbind(freqs, fst)

  }

  # Calculate probabilities ---------------------------------------------------

  if (p.values) {

    col.all  <- ncol(freq.alive) + 1
    col.end  <- ncol(freq.alive) + ncol(freq.all)
    probs  <- apply(cbind(freq.alive, freq.all), MARGIN = 1,
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
#' @param population A logical vector. TRUE belong to pop1, FALSE belong to
#' pop2. A M&C FST estimation will be included in the result.
#' @param min.freq.al Optional. Numeric from 0 to 1. Variants with MAF under
#' this threshold will be ommited.
#' @param location.cols If TRUE, adds two columns to the result, which are the
#' scaffold name and the locus name as factors (deafult is FALSE)
#' @param deltas If TRUE (default), include columns of differencies between
#' relative allelic frequencies after selection and before selection in the
#' result. Positive value means that the allele is more common in the survivers
#' population.
#' @param p.values If TRUE (default), some probabilities will be included in
#' the result, as well as frequencies
#' @param genotypic If TRUE, will include counts of genotype in the result.
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
#' \dontrun{
#' AnalyseExpt(gt, survival)
#' AnalyseExpt(gt, survival, min.freq.al = 0.1, location.cols = FALSE)
#' AnalyseExpt(gt, survival, deltas = FALSE, p.values = FALSE)
#' AnalyseExpt(gt, survival, backup.path = "example.csv")
#' }

AnalyseExpt <- function(lka, survival, population = NULL, min.freq.al = NULL,
                        location.cols = TRUE, deltas = TRUE, p.values = TRUE,
                        genotypic = FALSE, backup.path = NULL){

  # splitted.gt <- SplitGt(gt, survival)
  # gt.alive    <- splitted.gt$alive
  # gt.dead     <- splitted.gt$dead
  # gt.pop1     <- NULL
  # gt.pop2     <- NULL

  # if (is.logical(population)) {
  #   gt.pop  <- SplitGt(gt, population)
  #   gt.pop1 <- gt.pop$alive
  #   gt.pop2 <- gt.pop$dead 
  # }

  # result <- AnalyseSplittedExpt(gt.alive, gt.dead, min.freq.al = min.freq.al,
  #                               location.cols = location.cols, deltas = deltas,
  #                               p.values = p.values, genotypic = genotypic,
  #                               gt.pop1 = gt.pop1, gt.pop2 = gt.pop2,
  #                               backup.path = backup.path)
  # return(result)

  return("hallo")
}