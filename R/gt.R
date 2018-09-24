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


#' A function to extract a genotype from a vcfR::VCF object.
#'
#' @param vcf A vcfR::VCF object
#' @param min.depth Optional. A threshold depth under which the data will be
#' censored (possibly a positive Integer)
#' @param min.sample.qual Optional. A threshold ratio under which the samples
#' with too much missing data will be removed (from 0 to 1)
#' @param min.variant.qual Optional. A threshold ratio under which the variants
#' with too much missing data will be removed (from 0 to 1)
#' @param include.depth If TRUE, will return a list of matrix including the
#' genotype as \code{list$gt} and the depth as \code{list$dp}
#' @param verbose Logical. If TRUE (default), report status of the process
#' along the execution
#'
#' @return 
#' A genotype matrix from the vcfR object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ExtractGt(vcf)
#' ExtractGt(vcf, min.sample.qual = 0.65, min.variant.qual = 0.4)
#' ExtractGt(vcf, min.depth = 4, verbose = FALSE)
#' }

ExtractGt <- function(vcf, min.depth = NULL, min.sample.qual = NULL,
                      min.variant.qual = NULL, include.depth = FALSE,
                      verbose = TRUE) {

  # A matrix of strings (e.g. "1/0") or NA.
  gt <- vcfR::extract.gt(vcf, element = "GT")

  # A matrix of positive integers or NA.
  dp <- vcfR::extract.gt(vcf, element = "DP", as.numeric = TRUE)

  if (is.numeric(min.depth)) {
    if (min.depth %% 1 != 0 || min.depth < 0) {
      warning("Parameter min.depth should be a positive integer.")
    }
    censored <- CensorLowDepth(gt, dp, min.depth, verbose = verbose)
    gt <- censored$gt
    if (include.depth) {
      dp <- censored$dp
    }
  } else if (!is.null(min.depth)) {
    stop("Parameter min.depth must be a number or NULL.")
  }

  if (is.numeric(min.variant.qual)) {
    if (min.variant.qual > 1 || min.variant.qual < 0) {
      warning("Parameter min.variant.qual should be a float between 0 and 1")
    }
    omitted <- OmitPoorVariants(gt, min.variant.qual, dp = dp,
                                verbose = verbose)
    gt <- omitted$gt
    if (include.depth) {
      dp <- omitted$dp
    }
  } else if (!is.null(min.variant.qual)) {
    stop("Parameter min.variant.qual must be a number or NULL.")
  }

  if (is.numeric(min.sample.qual)) {
    if (min.sample.qual > 1 || min.sample.qual < 0) {
      warning("Parameter min.sample.qual should be a float between 0 and 1")
    }
    omitted <- OmitPoorSamples(gt, min.sample.qual, dp = dp, verbose = verbose)
    gt <- omitted$gt
    if (include.depth) {
      dp <- omitted$dp
    }
  } else if (!is.null(min.sample.qual)) {
    stop("Parameter min.sample.qual must be a number or NULL.")
  }

  gt <- data.frame(gt, stringsAsFactors = FALSE)
  dp <- data.frame(gt, stringsAsFactors = FALSE)
  
  if (include.depth) {
    result <- list("gt" = gt, "dp" = dp)
  } else {
    result <- gt
  }

  return(result)

}


#' A function to replace genotype data with low depth by NA
#'
#' @param gt A dataframe to process (e.g. a genotype dataframe)
#' @param dp A dataframe from which the depth of \code{gt} will be extracted
#' (e.g. a depth dataframe)
#' @param min.depth Integer. The minimum depth of each \code{gt} cells, under wich the
#' data will be replaced by NA
#' @param verbose Logical. If TRUE (default), report status of the process
#' along the execution
#'
#' @return 
#' A list of dataframes. 
#' \code{list$gt} is \code{gt} with NA instead of the data with low depth.
#' \code{list$dp} is \code{dp} with NA instead of the low values of depth.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' censorLowDepth(genotype, depth, 4)
#' censorLowDepth(genotype, depth, 4, verbose = FALSE)
#' }

CensorLowDepth <- function(gt, dp, min.depth, verbose = TRUE){

  if (verbose) {
    message("Censoring shallow data...")
  }

  filter <- (dp < min.depth)

  gt[filter] <- NA
  dp[filter] <- NA

  if (verbose) {
    message(paste(sum(filter, na.rm = TRUE), " genotype cells censored."))
  }

  result <- list("gt" = gt, "dp" = dp)

  return(result)

}


#' A function to remove columns with a high ratio of missing data from a
#' dataframe. (e.g. samples from genotype)
#'
#' @param gt A dataframe to process (e.g. a genotype dataframe)
#' @param min.qual The minimum valid data ratio under which the column
#' will be deleted from the dataframe (from 0 to 1)
#' @param dp An optionnal dataframe on which the process of \code{gt} will be
#' executed too (e.g. a depth dataframe)
#' @param verbose Logical. If TRUE (default), report status of the process
#' along the execution.
#'
#' @return 
#' A dataframe or a list of dataframes. 
#' If the \code{dp} is not set, return \code{gt} without its columns where the
#' missing data ratio was too high. 
#' If the \code{dp} is set, return the modified \code{gt} as \code{list$gt};
#' \code{list$dp} is the modified \code{dp} dataframe with the same missing 
#' columns as \code{list$gt}.
#' 
#' @seealso \code{\link{OmitPoorVariants}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' OmitPoorSamples(genotype, 0.15)
#' OmitPoorSamples(genotype, 0.2, dp = depth)
#' OmitPoorSamples(genotype, 0.2, verbose = FALSE)
#' }

OmitPoorSamples <- function(gt, min.qual, dp = NULL, verbose = TRUE) {

  if (verbose) {
    message("Omitting low quality samples...")
  }

  # FALSE if too much NAs per columns
  filter <- apply(gt, MARGIN = 2,
                  FUN = function(x){!IsAtLeastNARatio(x, 1 - min.qual)})

  gt <- gt[, filter]

  if (!is.null(dp)) { #if dp is set, process and return dp too
    dp <- dp[, filter]
    result <- list("gt" = gt, "dp" = dp)
  } else { #else just return gt
    result <- gt
  }

  if (verbose) {
    message(paste(sum(!filter, na.rm = TRUE), " samples omitted."))
  }

  return(result)

}


#' A function to remove rows with a high ratio of missing data from a
#' dataframe. (e.g. variants from genotype)
#'
#' @param gt A dataframe to process (e.g. a genotype dataframe)
#' @param min.qual The minimum present data ratio under which the row will be
#' deleted from the dataframe (from 0 to 1)
#' @param dp An optionnal dataframe on which the process of \code{gt} will be
#' executed too (e.g. a depth dataframe)
#' @param verbose Logical. If TRUE (default), report status of the process
#' along the execution.
#'
#' @return 
#' A dataframe or a list of dataframes. 
#' If the parameter \code{dp} is not set, return the \code{gt} dataframe
#' without its rows where the missing data ratio was too high. 
#' If the \code{dp} is set, return the modified \code{gt} as \code{list$gt};
#' \code{list$dp} is the modified \code{dp} dataframe with the same missing 
#' rows as \code{list$gt}.
#' 
#' @seealso \code{\link{OmitPoorSamples}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' OmitPoorVariants(genotype, 0.15)
#' OmitPoorVariants(genotype, 0.2, dp = depth)
#' OmitPoorVariants(genotype, 0.2, verbose = FALSE)
#' }

OmitPoorVariants <- function(gt, min.qual, dp = NULL, verbose = TRUE) {

  if (verbose) {
    message("Omitting low quality variants...")
  }

  #FALSE if too much NAs per rows
  filter <- apply(gt, MARGIN = 1,
                  FUN = function(x) {!IsAtLeastNARatio(x, 1 - min.qual)})

  gt <- gt[filter, ]

  if (!is.null(dp)) { #if dp is set, process and return dp too
    dp <- dp[filter, ]
    result <- list("gt" = gt, "dp" = dp)
  } else { #else just return gt
    result <- gt
  }

  if (verbose) {
    message(paste(sum(!filter, na.rm = TRUE), " variants omitted."))
  }

  return(result)

}


#' A function to split a genotype matrix two genotype matrix according to a
#' logical vector
#'
#' @param gt The original genotype matrix to split
#' @param survival A logical vector for splitting the data into two gt matrix.
#' The length of the vector must match the number of samples in the gt matrix.
#' @param verbose Logical. If TRUE (default), report status of the process
#' along the execution.
#'
#' @return 
#' A list of genotype matrix.
#' \code{list$alive} is the gt matrix containing the samples which position
#' corresponds with a TRUE in the \code{survival} vector.
#' The other samples are contained in \code{list$dead} as an other gt matrix
#'
#' @export
#'
#' @examples
#' \dontrun{
#' SplitGt(gt, c(TRUE, TRUE, FALSE, ...))
#' SplitGt(gt, survival.vector, verbose = FALSE)
#' }

SplitGt <- function(gt, survival, verbose = TRUE){

  alive <- NULL
  dead  <- NULL

  if (!is.logical(survival)) {
    stop("Parameter survival must be a logical vector")
  }
  if (class(gt) == "list") {
    return(lapply(gt, FUN = function(x){ SplitGt(x, survival) }))
  }
  if (class(gt) != "matrix" ) {
    stop("Parameter gt must be a matrix")
  }

  # number of alive samples
  len.survival <- length(survival)
  # number of samples in the gt object
  len.gt <- ncol(gt)

  if (len.survival < len.gt){
    if (verbose) {
      warning(paste("Parameter surivals is shorter than parameter gt (",
                    len.survival, " versus ", len.gt, ")"))
    }
  } else if (len.survival > len.gt) {
    stop(paste("Parameter surivals is longer than parameter gt (",
                    len.survival, " versus ", len.gt, ")"))
  }

  alive <- gt[, survival]
  dead  <- gt[, !survival]

  # force 1 row data frame outputs when outputs are 1 row vectors
  if (is.null(dim(alive))) {
    alive <- t(data.frame(alive))
    rownames(alive) <- rownames(gt)
  }
  if (is.null(dim(dead))) {
    dead <- t(data.frame(dead))
    rownames(dead) <- rownames(gt)
  }

  # $alive items correspond to a TRUE in the survival vector
  result <- list("alive" = alive, "dead" = dead) 

  return(result)

}