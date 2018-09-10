#' A function to extract a genotype from a vcfR::VCF object.
#'
#' @param vcf A vcfR::VCF object
#' @param min.depth Optional. A threshold depth under which the data will be
#' censored (possibly a positive Integer)
#' @param min.sample.qual Optional. A threshold ratio under which the samples
#' with too much missing data will be removed (from 0 to 1)
#' @param min.variant.qual Optional. A threshold ratio under which the variants
#' with too much missing data will be removed (from 0 to 1)
#' @param verbose Logical. If TRUE (default), report status of the process
#' along the execution
#'
#' @return 
#' A list of dataframes. \code{list$gt} is the genotype dataframe, \code{list$dp} is the
#' depth dataframe.
#'
#' @export
#'
#' @examples
#' ExtractGt(vcf)
#' ExtractGt(vcf, min.sample.qual = 0.65, min.variant.qual = 0.4)
#' ExtractGt(vcf, min.depth = 4, verbose = FALSE)

ExtractGt <- function(vcf, min.depth=NULL, min.sample.qual=NULL,
                      min.variant.qual=NULL, verbose=TRUE) {

  # A data.frame of strings (e.g. "1/0") or NA.
  gt <- vcfR::extract.gt(vcf, element = "GT")

  # A data.frame of positive integers or NA.
  dp <- vcfR::extract.gt(vcf, element = "DP", as.numeric = TRUE)

  if (is.numeric(min.depth)) {

    censored <- CensorLowDepth(gt, dp, min.depth, verbose = verbose)
    gt <- censored$gt
    dp <- censored$dp

  }

  if (is.numeric(min.variant.qual)) {

    omitted <- OmitPoorVariants(gt, min.variant.qual, dp = dp,
                                verbose = verbose)
    gt <- omitted$gt
    dp <- omitted$dp

  }

  if (is.numeric(min.sample.qual)) {

    omitted <- OmitPoorSamples(gt, min.sample.qual, dp = dp, verbose = verbose)
    gt <- omitted$gt
    dp <- omitted$dp

  }

  result <- list("gt" = gt, "dp" = dp)

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
#' censorLowDepth(genotype, depth, 4)
#' censorLowDepth(genotype, depth, 4, verbose = FALSE)

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
#' OmitPoorSamples(genotype, 0.15)
#' OmitPoorSamples(genotype, 0.2, dp = depth)
#' OmitPoorSamples(genotype, 0.2, verbose = FALSE)

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
#' OmitPoorVariants(genotype, 0.15)
#' OmitPoorVariants(genotype, 0.2, dp = depth)
#' OmitPoorVariants(genotype, 0.2, verbose = FALSE)

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