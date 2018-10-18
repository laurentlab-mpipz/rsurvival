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

#' @export

ShapeCountsR <- function(variant, genotypic = TRUE, allelic = FALSE,
                            absolute = TRUE, percentage = FALSE,
                            extrapolate.freq = TRUE, totals = TRUE,
                            min.freq.gt = NULL, min.freq.al = NULL) {

  result.al <- NULL
  result.gt <- NULL
  result    <- NULL 
  return.na <- FALSE

  if (length(absolute) == 1) {
    absolute <- rep(absolute, 2)
  }

 # actual counting using regex -----------------------------------------------

  levels.list <- list(HOMOREF = c("0/0", "0|0"), HETERO = c("1/0", "1|0", "0/1", "0|1"), HOMOALT = c("1/1", "1|1"), NA)
  fact <- factor(variant, levels=levels.list,order = TRUE, exclude=NULL)
  counts = as.vector(table(fact))
  counts[5] = sum(counts)

  # calculate frequencies of alleles ------------------------------------------

  counts.al <- c(counts[1] * 2 + counts[2], counts[3] * 2 + counts[2],
                 2 * counts[4])
  total.al <- sum(counts.al)

  if (extrapolate.freq) {
    freqs.al  <- c(counts.al[1] / sum(counts.al[1:2]),
                   counts.al[2] / sum(counts.al[1:2]),
                   counts.al[3] / total.al)
  } else {
    freqs.al <- counts.al / total.al   
  }

  # calculate frequencies of genotypes ----------------------------------------

  counts.gt <- c(counts[1], counts[2], counts[3], counts[4])
  total.gt  <- sum(counts.gt)

  if (extrapolate.freq) {
    freqs.gt  <- c(counts.gt[1] / sum(counts.gt[1:3]),
                   counts.gt[2] / sum(counts.gt[1:3]),
                   counts.gt[3] / sum(counts.gt[1:3]),
                   counts.gt[4] / total.gt)
  } else {
    freqs.gt  <- counts.gt / total.gt
  } 

  # check minimal frequencies -------------------------------------------------

  if (!is.null(min.freq.al)) {
    if (!is.numeric(min.freq.al)) {
      stop("Parameter min.freq.al must be a number.")
    }
    if ((min.freq.al < 0) || (min.freq.al > 1)) {
      warning("Parameter min.freq.al should be from 0 to 1.")
    }
    if (sum(freqs.al[1:2] < min.freq.al) > 0) {
      return.na = TRUE
    }
  } else if (!is.null(min.freq.gt)) {
    if (!is.numeric(min.freq.gt)) {
      stop("Parameter min.freq.gt must be a number.")
    }
    if ((min.freq.gt < 0) || (min.freq.gt > 1)) {
      warning("Parameter min.freq.gt should be from 0 to 1.")
    }
    if (sum(freqs.gt[1:3] < min.freq.gt) > 0) {
      return.na = TRUE
    }
  }

  if (return.na) {
      counts.gt <- rep(NA, length(counts.gt))
      freqs.gt  <- rep(NA, length(freqs.gt))
      total.gt  <- rep(NA, length(total.gt))
      counts.al <- rep(NA, length(counts.al))
      freqs.al  <- rep(NA, length(freqs.al))
      total.al  <- rep(NA, length(total.al))
  }

  # build resulting frequencies -----------------------------------------------

  if (!absolute[1]) { # for genotype
    prefix.gt    <- "freq"
    result.gt <- freqs.gt
    if (percentage) {
      prefix.gt <- "perc"
      result.gt <- result.gt * 100
    }
  } else {
    prefix.gt <- "count"
    result.gt <- counts.gt
  }

  if (!absolute[2]) { # for alleles
    prefix.al <- "freq"
    result.al <- freqs.al
    if (percentage) {
      prefix.al <- "perc"
      result.al <- result.al * 100
    }
  } else {
    prefix.al  <- "count"
    result.al <- counts.al
  }

  # build titles for each frequencies returned --------------------------------

  names(result.al) <- paste(prefix.al,
                            c("al.REF", "al.ALT", "al.MISSVAL"), sep=".")
  names(result.gt) <- paste(prefix.gt,
                            c("gt.HOMOREF", "gt.HETERO", "gt.HOMOALT",
                              "gt.MISSVAL"),
                            sep=".")
  names(total.al)  <- c("count.al.TOTAL")
  names(total.gt)  <- c("count.gt.TOTAL")

  # concatenate results if needed ---------------------------------------------

  if (totals) {
    result.al <- c(result.al, total.al)
    result.gt <- c(result.gt, total.gt)
  }

  if (genotypic && allelic) {
    result <- c(result.gt, result.al)
  } else if (genotypic) {
    result <- result.gt
  } else if (allelic) {
    result <- result.al
  } else {
    result <- NULL
  }

  if (include.ids) {
    result <- list("freq" = result, "ids" = ids)
  }

  return(result)

}