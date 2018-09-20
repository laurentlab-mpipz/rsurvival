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


#' A function to calculate probabilities that observed frequencies are due to
#' positive selection using MWNCHHypergeo model.
#'
#' @param freq.alive An vector of observed absolute frequencies for the three
#' genotypes in the alive population
#' @param freq.all An vector of observed absolute frequencies for the three
#' genotypes in the population before selection
#' @param map.alive Optional. Integer list containing the ids of specific
#' columns of \code{freq.alive}. Increase speed for large scale calls. Items of
#' \code{map.alive} must be named "homo.ref", "hetero", "homo.alt", "missval",
#' "total".
#' @param map.all Optional. Integer list containing the ids of specific
#' columns of \code{freq.all}. Increase speed for large scale calls. Items of
#' \code{map.all} must be named "homo.ref", "hetero", "homo.alt", "missval",
#' "total".
#' @return 
#' A vector of probabilities : \code{p_neutral} is the probability to obtain
#' this distribution with a random picking process, \code{p_select} is the
#' probability to obtain this distribution with a picking in a biased urn,
#' \code{weights} are the odds used to calculate \code{p_select}, \code{lrt}
#' is the result of the likelyhood test of \code{p_select} uppon
#' \code{p_neutral}, \code{p_value} is the probability that the distribution
#' is due to positive selection.
#'
#' @seealso For more information, see \code{\link{CalcWeightsSurvival}} and
#' \code{\link{FindIdsGtCounts}} wich this function binds.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' f.alive <- c(17, 3, 15, 1, 36)
#' f.all   <- c(24, 5, 22, 2, 53)
#' CalcProbsSelection(f.alive, f.all)
#' }

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
#' @param map.alive Optional. Integer list containing the ids of specific
#' columns of \code{freq.alive}. Increase speed for large scale calls. Items of
#' \code{map.alive} must be named "homo.ref", "hetero", "homo.alt", "missval",
#' "total".
#' @param map.all Optional. Integer list containing the ids of specific
#' columns of \code{freq.all}. Increase speed for large scale calls. Items of
#' \code{map.all} must be named "homo.ref", "hetero", "homo.alt", "missval",
#' "total".
#'
#' @return 
#' A vector of odds for the 3 genotypes, according to the MWNCHHypergeo model.
#'
#' @seealso For more information, see \code{\link{FindIdsGtCounts}} which this
#' function binds.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' f.alive <- c(17, 3, 15, 1, 36)
#' f.all  <- c(24, 5, 22, 2, 53)
#' CalcWeightsSurvival(f.alive, f.all)
#' }

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