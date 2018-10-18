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

Benchmark <- function(gt, variant.name, s, h, nb.dead, n, verbose = TRUE) {

  message("Preparing : Generate survival vectors.")

  surv <- GenerateSurvVector(gt[variant.name, ] , s, h, nb.dead, n)
  ordered <- list(rep(NA, n))
  variant.names <- rownames(gt)

  n.cores <- parallel::detectCores()
  cluster <- parallel::makeCluster(n.cores) # prevent overloading
  doSNOW::registerDoSNOW(cluster)

  time.start <- proc.time()[3] # ellapsed time since session launched
  time.stamp <- time.start
  # time.left  <- NA # estimated remaining time

  message("Preparing : Cluster is ready. (", n.cores, " cores used)")
  message("Iterating : Iterations start.")

  progress.bar <- txtTimerBar(n=n)
  progress <- function(i) setTxtProgressBar(progress.bar, i)
  options <- list(progress = progress)

  ordered <- foreach::foreach(i = 1:n, .options.snow = c("progress" = progress)) %dopar% {

    # operator %dopar% binds foreach::"%dopar%"
    # direct use of foreach::"%dopar%" will cause a package building failure
    # binding is made in file /R/config.R 

    # if (verbose) {
    #   message.start = paste("Iterating : ", i, " out of ", n, sep = "")
    #   if (i == 1) {
    #     message(paste(message.start , " (Estimating time)", sep = ""))
    #   } else {
    #     message(paste(message.start, " (", time.left, " remaining)", sep = ""))
    #   }
    # }
    
    if (class(gt) == "matrix") {
      gt <- data.frame(gt, stringsAsFactors = FALSE)
    } else if (class(gt) == "list") {
      return(lapply(gt, FUN = function(x){ SplitGt(x, survival) }))
    }

    alive <- gt[, surv[,i]]
    dead  <- gt[, !surv[,i]]

    # force 1 row data frame outputs when outputs are 1 row vectors
    if (is.null(dim(alive))) {
      alive <- t(data.frame(alive))
      rownames(alive) <- rownames(gt)
    }
    if (is.null(dim(dead))) {
      dead <- t(data.frame(dead))
      rownames(dead) <- rownames(gt)
    }

    #gt.alive    <- splitted.gt$alive
    #gt.dead     <- splitted.gt$dead

    analysis <- AnalyseSplittedExpt(alive, dead, min.freq.al = 0.1)

    ##analysis <- AnalyseExpt(gt, surv[, i], min.freq.al = 0.1)
    ord.analysis <- analysis[order(analysis[, "p.value"]), ]
    ord.analysis <- ord.analysis[!is.na(ord.analysis[, "p.value"]), ]
    ord.analysis <- ord.analysis[ord.analysis[, "p.value"] < 0.05, ]
    ord.analysis <- rownames(ord.analysis)
    ord.analysis

    # estimating remaining time from the last iteration time
    # time.left  <- round((proc.time()[3] - time.stamp) * (n - i))
    # time.left  <- lubridate::seconds_to_period(time.left) # format to time
    # time.stamp <- proc.time()[3]

  }

  close(progress.bar)
  parallel::stopCluster(cluster)

  if (verbose) {
    # time.stamp <- proc.time()[3]
    time.ellapsed <- round((proc.time()[3] - time.start))
    time.ellapsed <- lubridate::seconds_to_period(time.ellapsed)
    message("\nIterating : DONE. (",  time.ellapsed, " ellapsed)")
  }

  message("Ending : Data shaper launch.")

  positions <- as.data.frame(lapply(ordered,
                                   FUN = function(x){
                                     match(variant.names, x)
                                   })
                            )
  pos.variant <- match(variant.name, variant.names)
  result <- apply(positions, MARGIN = 2, FUN = function(x){x[pos.variant]})
  result <- result[!is.na(result)]

  message("Ending : DONE.")

  return(result)
  pos.mod <- apply(positions, MARGIN =  1, FUN = mod, na.rm = TRUE)
  names(pos.mod) <- variant.names
  result <- pos.mod[order(pos.mod)]

  return(result)

}

BenchmarkSH <- function(variant, s, h, nb.dead, n){
  
  surv <- GenerateSurvVector(variant , s, h, nb.dead, n)

  analysis <- apply(surv, MARGIN = 2,
                    FUN = function(x){
                      AnalyseExpt(variant, x)
                    })
  
  s <- unlist(lapply(analysis,
                     FUN = function(x){
                      x$s
                     }))

  
  h <- unlist(lapply(analysis,
                   FUN = function(x){
                    x$h
                   }))

  return(list("s" = s, "h" = h))
}

FindSH <- function(benchmark){

  s.vect <- benchmark$s
  h.vect <- benchmark$h

  censore <- is.na(s.vect) || is.na(h.vect)
  s.vect <- s.vect[!censore]
  h.vect <- h.vect[!censore]

  # s.density <- density(s.vect)
  # h.density <- density(h.vect)

  # s <- s.density$x[match(max(s.density$y), s.density$y)]
  # h <- h.density$x[match(max(h.density$y), h.density$y)]

  s <- mod(s.vect)
  h <- mod(h.vect)

  return(list("s" = s, "h" = h))

}

TestPb <- function(n){
  progress.bar <- txtProgressBar(max = n, style = 3)
  progress <- function(i) setTxtProgressBar(progress.bar, i)

  for (i in 1:n) {
    progress(i)
  }
  close(progress.bar)
}

#' Text progress bar with time.
#'
#' A textual progress bar that estimates time remaining. It displays the
#' estimated time remaining and, when finished, total duration.  Lifted from the
#' plyr package.  Update with \code{setTxtProgressBar}
#'
#' @export
txtTimerBar <- function(n = 1) {
  start <- .last_update_time <- proc.time()[3]
  times <- numeric(n)
  value <- NULL

  killed <- FALSE

  width <- getOption("width") - nchar('||100%  (9999 H 99 M 99 S remaining)')

  update <- function(i) {
    if (i == 0) return()

    value <<- i
    times[i] <- proc.time()[3] - start

    avg <- times[i] / i
    time_left <- (n - i) * avg

    nbars <- trunc(i / n * width)

    cat_line("|", str_rep("=", nbars), str_rep(" ", width - nbars), "|",
      format(i / n * 100, nsmall=0, digits=4, width = 3), "% (",
             show_time(time_left), " remaining)")
    #Note that I added `nsmall` and `digits` above, - NR 5/10/13
  }
  getVal <- function() value
  kill <- function(){
    if (killed) return()
    killed <<- TRUE

    if (value == n) {
      cat_line("|", str_rep("=", width), "|100%")
      #cat("\nCompleted after", show_time(proc.time()[3] - start), "\n")
    } else {
      cat("\nKilled after", show_time(proc.time()[3] - start), "\n")
    }
  }

  cat_line("|", str_rep(" ", width), "|  0%")

  structure(
    list(getVal = getVal, up = update, kill = kill),
    class = "txtProgressBar")
}

show_time <- function(x) {
  lubridate::seconds_to_period(round(x))
}

cat_line <- function(...) {
  msg <- paste(..., sep = "", collapse = "")
  gap <- max(c(0, getOption("width") - nchar(msg, "width")))
  cat("\r", msg, rep.int(" ", gap), sep = "")
  flush.console()
}

str_rep <- function(x, i) {
  paste(rep.int(x, i), collapse = "")
}

#' @export

BenchmarkRMSE <- function(gt, variant.name, surv.ratio){

  n.iter <- 100
  s <- 2
  h <- 0.5
  n.surv.vect <- round(ncol(gt) * surv.ratio)
  result <- list(rep(NA, n.iter))
  i <- 1

  for (n.surv in n.surv.vect) {

    message("\n########## SITUATION ", i, " OUT OF ", length(n.surv.vect), " ##########")
    result[i] <- list(Benchmark(gt, variant.name, s, h, n.surv, n.iter))
    #result[i] <- c(temp, rep(NA, n.iter - length(temp)))
    message("########## END OF SITUATION ", i, " ##########")
    i <- i + 1

  }

  return(result)

}