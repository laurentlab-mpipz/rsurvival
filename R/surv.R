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

GenerateSurvVector <- function(variant, s, h, nb.dead, nb.vect = 1) {

  f.variant <- CalcFreqVariant(variant, include.ids = TRUE);
  freqs     <- f.variant$freqs
  ids       <- f.variant$ids 

  nb.HMR <- freqs["count.gt.HOMOREF"]
  nb.HET <- freqs["count.gt.HETERO"]
  nb.HMA <- freqs["count.gt.HOMOALT"]

  surv.odds <- ConvertSHToOdds(s,h)
  pred.odds <- 1 / surv.odds

  deaths <- BiasedUrn::rMWNCHypergeo(nran = nb.vect,
                                     m = c(nb.HMR, nb.HET, nb.HMA),
                                     n = nb.dead, odds = pred.odds)

  if (is.null(dim(deaths))) {
  
      result <- ConvertNbDeathsToLogic(deaths, ids, length(variant))
  
  } else {

    result <- apply(deaths, MARGIN = 2,
                    FUN = function(x) {
                      ConvertNbDeathsToLogic(x, ids, length(variant))
                    })
  
  }

  return(result)

}

#' @export

ConvertOddsToSH <- function(odds) {

  homo.odds     <- c(odds[1], odds[3])
  min.homo.odds <- min(homo.odds)
  sel.odds      <- (1 / min.homo.odds) * odds # set min w to 1
  homo.odds     <- c(sel.odds[1], sel.odds[3])
  max.homo.odds <- max(homo.odds)
  hetero.odd    <- sel.odds[2]

  s <- max.homo.odds - 1
  h <- (hetero.odd - 1) / s

  result <- list("s" = s, "h" = h)

  return(result)

}

#' @export

ConvertSHToOdds <- function(s, h) {

  hetero.odd    <- 1 + s*h
  max.homo.odds <- 1 + s
  odds          <- rep(0, 3)
  odds[1]       <- 1
  odds[2]       <- hetero.odd
  odds[3]       <- max.homo.odds

  result <- odds

  return(result)

}

ConvertNbDeathsToLogic <- function(deaths, ids, length) {

  deaths.HOMOREF <- deaths[1]
  deaths.HETERO  <- deaths[2]
  deaths.HOMOALT <- deaths[3]

  surv.HOMOREF <- c(rep(FALSE, deaths.HOMOREF),
                    rep(TRUE, length(ids$HOMOREF) - deaths.HOMOREF))
  surv.HETERO  <- c(rep(FALSE, deaths.HETERO),
                    rep(TRUE, length(ids$HETERO) - deaths.HETERO))
  surv.HOMOALT <- c(rep(FALSE, deaths.HOMOALT),
                    rep(TRUE, length(ids$HOMOALT) - deaths.HOMOALT))

  ids$HOMOREF <- ids$HOMOREF[sample(surv.HOMOREF)]
  ids$HETERO  <- ids$HETERO[sample(surv.HETERO)]
  ids$HOMOALT <- ids$HOMOALT[sample(surv.HOMOALT)]

  result <- rep(FALSE, length)
  result[unlist(ids)] <- TRUE

  return(result)

}
