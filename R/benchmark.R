
#' @export

Benchmark <- function(gt, variant.name, s, h, n.dead, n, verbose = TRUE,
                      shape = TRUE, max.cores = NULL) {

  message("Preparing : Generate survival vectors.")

  surv <- GenerateSurvVector(gt[variant.name, ] , s, h, n.dead, n)
  ordered <- list(rep(NA, n))
  variant.names <- rownames(gt)

  # Handle multiple cores processing ------------------------------------------

  if (is.null(max.cores)) {
    n.cores <- parallel::detectCores()
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

  progress.bar <- txtTimerBar(n = n)
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

    analysis <- AnalyseSplittedExpt(alive, dead, min.freq.al = 0.1)

    ord.analysis <- analysis[order(analysis[, "p.value"]), ]
    ord.analysis <- ord.analysis[!is.na(ord.analysis[, "p.value"]), ]
    ord.analysis <- ord.analysis[ord.analysis[, "p.value"] < 0.05, ]

    rownames(ord.analysis)

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




#' @export

BenchmarkSH <- function(variant, s, h, n.dead, n){
  
  surv <- GenerateSurvVector(variant , s, h, n.dead, n)

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




#' @export

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


#' @export

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

BenchmarkScan <- function(gt, variant.name, surv.ratio = 0.5, s = 2, h = 0.5, n.iter = 100, shape = TRUE, max.cores = NULL){

  len.s <- length(s)
  len.h <- length(h)
  len.r <- length(surv.ratio)

  if (sum(c(len.r > 1, len.s > 1, len.h > 1)) > 1) {
    stop("Only one parameter can variate during a scan.")
  }

  n.dead <- round(ncol(gt) * (1 / surv.ratio))
  result <- list(rep(NA, n.iter))

  if (len.r > 1) {
    variations <- n.surv
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
                                  surv.ratio = variation, n = n.iter,
                                  shape = shape, max.cores = max.cores))
    } else if (len.s > 1) {
      result[i] <- list(Benchmark(gt, variant.name, s = variation, h = h,
                                  n.dead = n.dead, n = n.iter, shape = shape,
                                  max.cores = max.cores))
    } else if (len.h > 1) {
      result[i] <- list(Benchmark(gt, variant.name, s = s, h = variation,
                                  n.dead = n.dead, n = n.iter, shape = shape,
                                  max.cores = max.cores))
    } else {
      warning("No variations has been detected. Aborting.")
      return(NA)
    }

    message("############ END OF SITUATION ", i, " ############")

    i <- i + 1

  }

  return(result)

}



#' @export

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