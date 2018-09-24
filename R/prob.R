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
#' @seealso For more information, see \code{\link{CalcWeightsPredation}} and
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

  pred.weights <- CalcWeightsPredation(freq.alive, freq.all,
                                       map.alive = map.alive,
                                       map.all = map.all)

  if ((sum(is.na(pred.weights)) == 0) && (sum(pred.weights <= 0) == 0) 
      && (sum(is.infinite(pred.weights)) == 0)) {

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

    homo.weights     <- c(pred.weights[1], pred.weights[3])
    min.homo.weights <- min(homo.weights)
    pred.weights     <- (1 / min.homo.weights) * pred.weights # set min w to 1
    homo.weights     <- c(pred.weights[1], pred.weights[3])
    max.homo.weights <- max(homo.weights)
    hetero.weight    <- pred.weights[2]

    s <- max.homo.weights - 1
    h <- hetero.weight + s - 1

    # distribution in survivors population
    x <- freq.alive[c(id.alive$homo.ref, id.alive$hetero, id.alive$homo.alt)]
    # distribution in original population
    m <- freq.all[c(id.all$homo.ref, id.all$hetero, id.all$homo.alt)]
    # nb of survivivor samples with valid data
    n <- freq.alive[id.alive$total] - freq.alive[id.alive$missval]

    prob.sel     <- BiasedUrn::dMWNCHypergeo(x = x, m = m, n = n,
                                          odds = 1 / pred.weights)
    prob.neutral <- BiasedUrn::dMWNCHypergeo(x = x, m = m, n = n,
                                              odds = c(1,1,1))

    lrt     <- 2 * log(prob.sel / prob.neutral)
    p.value <- 1 - stats::pchisq(q = lrt, df = 2)
    result  <- c(pred.weights, s, h, prob.neutral, prob.sel, lrt, p.value)

  } else {
    result <- rep(NA, 9) # to ensure that the result always has the same length
  }

  names(result) <- c("weight.gt.HOMOREF", "weight.gt.HETERO",
                      "weight.gt.HOMOALT", "s", "h", "p.neutral", "p.select",  "lrt",
                      "p.value")
  return(result)

}


#' A function to calculate odds for predation from observed selection
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
#' CalcWeightsPredation(f.alive, f.all)
#' }

CalcWeightsPredation <- function(freq.alive, freq.all, map.alive = NULL,
                                map.all = NULL){

    standard.map <- list("homo.ref" = 1 , "hetero" = 2, "homo.alt" = 3,
                         "missval" = 4, "total" = 5)

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


    freq.dead <- freq.all[unlist(id.all)] - freq.alive[unlist(id.alive)]
    id.dead   <- standard.map

    # distribution of dead samples
    mu <- freq.dead[c(id.dead$homo.ref, id.dead$hetero, id.dead$homo.alt)]
    # distribution of original population
    m  <- freq.all[c(id.all$homo.ref, id.all$hetero, id.all$homo.alt)] 
    # number of dead samples with valid data
    n  <- freq.dead[id.dead$total] - freq.dead[id.dead$missval]

    # calculate odds
    result <- suppressWarnings(BiasedUrn::oddsMWNCHypergeo(mu = mu,
                                                            m = m, n = n))
    return(result)

}