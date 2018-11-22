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

#' A function to count occurences of genotypes in multiple combined variant
#' vectors. 
#'
#' @param variants A vector of variants, subset of a genotype matrix.
#' @param group.nas If FALSE (default), each combinations with missing values
#' will be considered on its own. Else, combinations with missing values will
#' be grouped as one category.
#'
#' @return 
#' A vector of frequencies.
#' Each Item is a calculated frequency for this combination of variants.
#'
#' @seealso \code{\link{CalcFreqGt}} which calculate frequencies for a 
#' genotype matrix.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' CalcFreqMultiVariant(variants, group.nas = TRUE)
#' }

CalcFreqMultiVariant <- function(variants, group.nas = FALSE){

  sorted <- SortSamplesVariant(variants, group.nas = group.nas)
  result <- unlist(lapply(sorted, FUN = length))
  result[["TOTAL"]] <- sum(result)
  names(result) <- paste("count.gt.", names(result), sep = "")

  return(result)

}


#' A function to sort occurences of genotypes in multiple combined variant
#' vectors.
#'
#' @param variants A vector of variants, subset of a genotype matrix.
#' @param group.nas If FALSE (default), each combinations with missing values
#' will be considered on its own. Else, combinations with missing values will
#' be grouped as one category.
#'
#' @return 
#' A list of vectors which contains the ids of categorized samples.
#'
#' @seealso \code{\link{CalcFreqMultiVariant}} which binds this function
#'
#' @export
#'
#' @examples
#' \dontrun{
#' SortSamplesMultiVariant(variants, group.nas = TRUE)
#' }

SortSamplesVariant <- function(variants, group.nas = FALSE){

  result     <- NULL
  n.variants <- nrow(variants)
  names      <- names(IdentifyGt(SliceDfRows(variants, 1)))

  if (n.variants == 1) {

    ids   <- IdentifyGt(SliceDfRows(variants, 1))
    names <- names(ids)
    gt    <- lapply(ids,
                    FUN = function(x){
                      SliceDfColumns(variants, x)
                    })

    result <- gt

  } else {

    for (i in 1:n.variants) {

      gt <- SortSamplesVariant(SliceDfRows(variants, i))
      names(gt) <- paste("VAR", i, "_", names(gt), sep = "")

      if (length(result)){
        gt <- lapply(result,
                     FUN = function(x){
                       res.x <- lapply(gt,
                                       FUN = function(y){
                                         filter <- colnames(x) %in% colnames(y)
                                         res.y  <- SliceDfColumns(x, filter)
                                       })
                     })
        gt <- unlist(gt, recursive = FALSE)
      }

      result <- c(gt)

    }

  }

  if (group.nas) {

    filter <- !grepl("MISSVAL", names(result))
    temp <- result
    result <- result[filter]
    result[["MISSVAL"]] <- data.frame(t(unlist(temp[!filter])))
  
  }

  return(result)

}



#' @export

CalcProbsMultiVariant <- function(freq.alive, freq.all, map.alive, map.all){

  freq.dead <- freq.all - freq.alive

  # distribution of dead samples
  mu <- freq.dead
  # distribution of original population
  m  <- freq.all 
  # number of dead samples with valid data
  n  <- sum(freq.dead)

  # replace 0s with non-zero values
  mu[mu == 0] <- 10^(-12)
  m[m == 0]   <- 10^(-12) # .Machine$double.xmin is to low

  pred.odds <- suppressWarnings(BiasedUrn::oddsMWNCHypergeo(mu = mu,
                                                       		m = m,
                                                       		n = n))  

  if (sum(is.na(pred.odds)) != 0) {
    stop("Enable to calculate probabilities for these frequencies.")
  }
  
  sel.odds <- 1 / pred.odds

  # distribution in survivors population
  x <- freq.alive # [c(map.alive[c(-length(map.alive), -length(map.alive) + 1)])]
  # distribution in original population
  m <- freq.all # [c(map.all[c(-length(map.all), -length(map.all) + 1)])]
  # nb of survivivor samples with valid data
  n <- sum(freq.alive) # [map.alive[length(map.alive)]] - freq.alive[map.alive[length(map.alive)]]

  prob.sel     <- BiasedUrn::dMWNCHypergeo(x = x, m = m, n = n,
                                           odds = 1 / pred.odds)
  prob.neutral <- BiasedUrn::dMWNCHypergeo(x = x, m = m, n = n,
                                           odds = rep(1, length(pred.odds)))

  lrt     <- 2 * log(prob.sel / prob.neutral)
  p.value <- 1 - stats::pchisq(q = lrt, df = 2)
  result  <- c(prob.neutral, prob.sel, lrt, p.value)

  names(result) <- c("p.neutral", "p.select",  "lrt", "p.value")
  return(result)

}
