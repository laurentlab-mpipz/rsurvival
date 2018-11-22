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

#' A function to perform multiple analysis using the same parameters.
#'
#' @param gt A genotype matrix.
#' @param variant.name Name of odded variant (e.g. scaffolf123.987654).
#' @param s Fitness parameter s for \code{variant.name}.
#' @param h Fitness parameter h for \code{variant.name}.
#' @param n.dead Number of dead samples for each iteration.
#' @param n Number of iterations. \code{n} analysis will be performed with the
#' same parameters.
#' @param max.p.value Threshold to identify significant SNPs. Default is 0.05.
#' @param min.freq.al OPTIONAL. Numeric from 0 to 1. Variants with MAF under
#' this threshold will be ommited.
#' @param verbose If TRUE (default), info on progression will be printed.
#' @param shape If TRUE (default), variants ranking will be returned. Else,
#' variants will be ordered and returned.
#' @param include.all If FALSE (default), only names of the variants will be
#' returned. Else, all the calculated frequencies and probs will be returned.
#' @param max.cores OPTIONAL. Number of cores from the CPU that R can use for
#' this job. If NULL, all detected cores will be used. 
#' 
#' @return 
#' Ordered variants names or a vector of variants rank. Will return this 
#' \code{n} times.
#'
#' @seealso \code{\link{BenchmarkScan}} which binds this function.
#'
#' @export

Benchmark <- function(gt, variant.name, s, h, n.dead, n, min.freq.al = 0,
                      verbose = TRUE, shape = TRUE, include.all = FALSE,
                      max.cores = NULL) {

  message("Preparing : Generate survival vectors.")

  surv <- GenerateSurvVector(gt[variant.name, ] , s, h, n.dead, n)
  ordered <- list(rep(NA, n))
  variant.names <- rownames(gt)

  # Handle multiple cores processing ------------------------------------------

  if (is.null(max.cores)) {
    n.cores <- max(1, parallel::detectCores() - 1)
  } else {
    n.cores <- min(max.cores, parallel::detectCores())
  } 

  cluster <- parallel::makeCluster(n.cores) # prevent overloading with -1
  doSNOW::registerDoSNOW(cluster)

  message("Preparing : Cluster is ready. (", n.cores, " cores used)")

  # Set timer and progress bar ------------------------------------------------

  time.start <- proc.time()[3] # ellapsed time since session launched
  time.stamp <- time.start
  
  message("Iterating : Iterations start.")

  progress.bar <- TxtTimerBar(n = n)
  progress <- function(i) setTxtProgressBar(progress.bar, i)
  options <- c("progress" = progress)

  # Iterations ----------------------------------------------------------------

  ordered.variants <- foreach::foreach(i = 1:n,
                                       .options.snow = options) %dopar% {

    # NOTE : Operator [%dopar%] binds [foreach::"%dopar%"].
    # Direct use of [foreach::"%dopar%"] will cause a package building failure.
    # Binding is made in file [/R/config.R].
    
    # Split data --------------------------------------------------------------
    if (class(gt) == "matrix") {
      gt <- data.frame(gt, stringsAsFactors = FALSE)
    } else if (class(gt) == "list") {
      return(lapply(gt, FUN = function(x){ SplitGt(x, survival) }))
    }

    alive <- gt[, surv[,i]]
    dead  <- gt[, !surv[,i]]

    # NOTE : This is a replicate of rsurvival::AnalyseExpt. 
    # It looks like rsurvival::SplitGt do not work when called in a multiple
    # cores loop. This function has been replaced here by brackets.

    # Shape data as dataframe -------------------------------------------------

    if (is.null(dim(alive))) {
      alive <- t(data.frame(alive))
      rownames(alive) <- rownames(gt)
    }
    if (is.null(dim(dead))) {
      dead <- t(data.frame(dead))
      rownames(dead) <- rownames(gt)
    }

    # Analyse and order variants ----------------------------------------------

    analysis <- AnalyseSplittedExpt(alive, dead, min.freq.al = min.freq.al,
                                    location.cols = FALSE, deltas = FALSE)
    
    ord.analysis <- analysis[order(analysis[, "p.value"]), ]
    ord.analysis <- ord.analysis[!is.na(ord.analysis[, "p.value"]), ]
    ord.analysis <- ord.analysis[ord.analysis[, "p.value"] < max.p.value, ]

    if (include.all) {
      ord.analysis
    } else {
      rownames(ord.analysis)
    }

  }

  # End multiple cores processing and progress bar ----------------------------

  close(progress.bar)     
  parallel::stopCluster(cluster)

  if (verbose) {
    time.ellapsed <- round((proc.time()[3] - time.start))
    time.ellapsed <- lubridate::seconds_to_period(time.ellapsed)
    message("\nIterating : DONE. (",  time.ellapsed, " ellapsed)")
  }

  # Shaping data --------------------------------------------------------------

  message("Ending : Data shaper launch.")

  if (shape) {

    if (sum(unlist(lapply(ordered.variants, # if items from list have rownames
                          FUN = function(x){
                            !is.null(rownames(x))
                          })))) {
      ordered.variants <- rownames(ordered.variants) # only keep rownames
    }

    print(length(ordered.variants))
    result <- ShapeAsVariantPos(ordered.variants, variant.name)

  } else {

    result <- ordered.variants

  }

  message("Ending : DONE.")

  return(result)
  
  # pos.mod <- apply(positions, MARGIN =  1, FUN = mod, na.rm = TRUE)
  # names(pos.mod) <- variant.names
  # result <- pos.mod[order(pos.mod)]

  # return(result)

}



# BenchmarkSHL <- function(variant, s, h, n.dead, n){
  
#   surv <- GenerateSurvVector(variant , s, h, n.dead, n)

#   analysis <- apply(surv, MARGIN = 2,
#                     FUN = function(x){
#                       AnalyseExpt(variant, x)
#                     })
  
#   s <- unlist(lapply(analysis,
#                      FUN = function(x){
#                       x$s
#                      }))

  
#   h <- unlist(lapply(analysis,
#                      FUN = function(x){
#                        x$h
#                      }))

#   l <- unlist(lapply(analysis,
#                      FUN = function(x){
#                       x$LOCUS
#                      }))

#   return(list("s" = s, "h" = h, "l" = l))
#}


#' A function to retrieve s and h fitness parameters from benchmark results.
#'
#' @param benchmark Results of \code{link{Benchmark}}
#' 
#' @return 
#' Return the mod of estimated s and h fitness parameters from benchmark
#' results.
#'
#' @seealso \code{\link{Benchmark}} needed for \code{benchmark} parameter.

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

#' A function to create a text progress bar.
#'
#' @param n Maximum value of progression.
#' 
#' @return 
#' Print a progress bar, use .update(progress) to update and .kill() to kill.

TxtTimerBar <- function(n = 1) {
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




#' A function to perform statistical bencharking for a range of parameters.
#' Launch multiple \code{link{Benchmark}}.
#'
#' @param gt A genotype matrix.
#' @param variant.name Name of odded variant (e.g. scaffolf123.987654).
#' @param surv.ratio Survival ratio from 0 to 1. Could be a vector.
#' @param s Fitness parameter s for \code{variant.name}. Could be a vector.
#' @param h Fitness parameter h for \code{variant.name}. Could be a vector.
#' @param n.dead Number of dead samples for each iteration.
#' @param n.iter Number of iterations. \code{n} analysis will be performed for
#' each combinations of parameters.
#' @param include.all If FALSE (default), only names of the variants will be
#' returned. Else, all the calculated frequencies and probs will be returned.
#' @param max.cores OPTIONAL. Number of cores from the CPU that R can use for
#' this job. If NULL, all detected cores will be used. 
#' 
#' @return 
#' Ordered variants names, calculated \code{n.iter} times for each combinations of
#' parameters. Only one parameter can variate during a scan. For multiple
#' values of multiple parameters, please make a loop on top of
#' \code{link{BenchmarkScan}}.
#'
#' @seealso \code{\link{BenchmarkScan}} which this function binds.
#'
#' @export

BenchmarkScan <- function(gt, variant.name, surv.ratio = 0.5, s = 2, h = 0.5,
                          n.iter = 100, shape = TRUE, include.all = FALSE,
                          max.cores = NULL){

  len.s <- length(s)
  len.h <- length(h)
  len.r <- length(surv.ratio)

  if (sum(c(len.r > 1, len.s > 1, len.h > 1)) > 1) {
    stop("Only one parameter can variate during a scan.")
  }

  n.dead <- round(ncol(gt) * (1 - surv.ratio))
  result <- list(rep(NA, n.iter))

  if (len.r > 1) {
    variations <- n.dead
  } else if (len.s > 1) {
    variations <- s
  } else if (len.h > 1) {
    variations <- h
  } else {
    warning("No variations has been detected. Aborting.")
    return(NA)
  }

  i <- 1

  for (variation in variations) {

    message("\n########## SITUATION ", i, " OUT OF ", length(variations),
            " ##########")

    if (len.r > 1) {
      result[i] <- list(Benchmark(gt, variant.name, s = s, h = h,
                                  n.dead = variation, n = n.iter,
                                  shape = shape, max.cores = max.cores,
                                  include.all = include.all))
    } else if (len.s > 1) {
      result[i] <- list(Benchmark(gt, variant.name, s = variation, h = h,
                                  n.dead = n.dead, n = n.iter, shape = shape,
                                  max.cores = max.cores,
                                  include.all = include.all))
    } else if (len.h > 1) {
      result[i] <- list(Benchmark(gt, variant.name, s = s, h = variation,
                                  n.dead = n.dead, n = n.iter, shape = shape,
                                  max.cores = max.cores,
                                  include.all = include.all))
    } else {
      warning("No variations has been detected. Aborting.")
      return(NA)
    }

    message("############ END OF SITUATION ", i, " ############")

    i <- i + 1

  }

  return(result)

}



#' A function to find the ranking of a variant in a list of ordered variants.
#'
#' @param ordered.variants A list of ordered variants name.
#' @param variant.name The name of a variant to look for.
#' 
#' @return 
#' The rank of \code{variant.name} in \code{ordered.variant}.
#'
#' @seealso \code{\link{Benchmark}} which binds this function.

ShapeAsVariantPos <- function(ordered.variants, variant.name){

  result <- as.data.frame(lapply(ordered.variants,
                                   FUN = function(x){
                                     match(variant.name, x) # variant.names
                                   })
                            )
  # pos.variant <- match(variant.name, variant.names)
  # result <- apply(positions, MARGIN = 2, FUN = function(x){x[pos.variant]})
  result <- result[!is.na(result)]

}