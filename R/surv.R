GenerateSurvVector <- function(variant, s, h, nbSurv) {

    return(NULL)

}

#' @export

ConvertOddsToSH <- function(odds) {

    homo.odds     <- c(odds[1], odds[3])
    min.homo.odds <- min(homo.odds)
    sel.odds      <- (1 / min.homo.odds) * odds # set min w to 1
    homo.odds     <- c(odds[1], odds[3])
    max.homo.odds <- max(homo.odds)
    hetero.odd    <- odds[2]

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
    odds[1]       <- max.homo.odds
    odds[2]       <- hetero.odd
    odds[3]       <- 1

    result <- odds

    return(result)

}