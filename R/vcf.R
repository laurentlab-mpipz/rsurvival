#' A function to get a filtered (and possibly splitted) vcfR::VCF object from
#' a VCF file.
#'
#' @param file.path The path to the VCF file you want to load
#' @param survivals Optional. A logical vector for splitting the data into two
#' vcfR::VCF objects. The length of the vector must match the number of samples
#' in the VCF file.
#' @param only.biallelic Logical. If TRUE (default), only biallelic variants
#' will be loaded.
#' @param only.snp Logical. If TRUE (default), the indel variants will be
#' omited.
#' @param verbose Logical. If TRUE (default), report status of the process
#' along the execution.
#'
#' @return 
#' A vcR:VCF object or a list of vcfR:VCF objects if a survivals vector was
#' provided. 
#' If a list is returned, list$alive is the vcfR::VCF object containing the
#' samples which the position corresponded with a TRUE in the survivals vector.
#' The other samples are contained in list$dead as an other vcfR::VCF object. 
#'
#' @export
#'
#' @seealso For more information, see \code{\link{SplitVcf}},
#' \code{\link{RmNonBiallelics}} and \code{\link{RmIndels}} which this
#' function bind.
#'
#' @examples
#' LoadVCF("example.vcf")
#' LoadVCF("example.vcf", survivals = c(TRUE, FALSE, FALSE, TRUE, ...))
#' LoadVCF("example.vcf", survivals = mySurvVector, verbose = FALSE)
#' LoadVCF("example.vcf", only.biallelic = FALSE, only.snp = FALSE)

LoadVcf <- function(file.path, survivals = NULL, only.biallelic = TRUE,
                    only.snp = TRUE, verbose = TRUE) {

  vcf <- NULL

  if (!is.character(file.path)){
    stop("Parameter file.path must be a character vector")
  } else if (!grepl("^.*\\.vcf$", file.path)){
    stop("Parameter file.path leads to a non VCF file")
  } else if (!is.logical(only.biallelic)){
    stop("Parameter only.biallelic must be a logical")
  } else if (!is.logical(only.snp)){
    stop("Parameter only.snp must be a logical")
  }

  if (file.exists(file.path)) {
    
    vcf <- vcfR::read.vcfR(file.path, verbose = verbose)

    if (is.logical(survivals) && !is.na(survivals)) {
      vcf <- SplitVcf(vcf, survivals, verbose = verbose)
    }
    if (only.biallelic) {
      vcf <- RmNonBiallelics(vcf, verbose = verbose)
    }
    if (only.snp) {
      vcf <- RmIndels(vcf, verbose = verbose)
    }

  } else {

    stop("File not found following parameter file.path")

  }

  return(vcf)

}


#' A function to split a vcfR::VCF object in two vcfR::VCF object according to
#' a logical vector
#'
#' @param vcf The original vcfR::VCF object to split
#' @param survivals A logical vector for splitting the data into two vcfR::VCF
#' objects.
#' The length of the vector must match the number of samples in the VCF file.
#' @param verbose Logical. If TRUE (default), report status of the process
#' along the execution.
#'
#' @return 
#' A list of vcfR:VCF objects.
#' \code{list$alive} is the vcfR::VCF object contains the samples which the position
#' corresponded with a TRUE in the \code{survivals} vector.
#' The other samples are contained in \code{list$dead} as an other vcfR::VCF object.
#'
#' @export
#'
#' @examples
#' SplitVcf(vcf, c(TRUE, TRUE, FALSE, ...))
#' SplitVcf(vcf, survivals.vector, verbose = FALSE)

SplitVcf <- function(vcf, survivals, verbose=TRUE){

  alive <- NULL
  dead  <- NULL

  if (!is.logical(survivals)) {
    stop("Parameter survivals must be a logical vector")
  }
  if (class(vcf) != "vcfR" ) {
    stop("Parameter vcf must be a vcfR object")
  }

  # number of alive samples
  len.survival <- length(survivals)
  # number of samples in th VCF file, "-1" uncount the FORMAT column
  len.vcf <- length(vcf[1, ]@gt) - 1 

  if (len.survival < len.vcf){
    if (verbose) {
      warning(paste("Parameter surivals is shorter than parameter vcf  (",
                    len.survival, " versus ", len.vcf, ")"))
    }
  } else if (len.survival > len.vcf) {
    stop(paste("Parameter surivals is longer than parameter vcf  (",
                    len.survival, " versus ", len.vcf, ")"))
  }

  # add TRUE at the beggining to include the FORMAT column as well
  alive <- vcf[, c(TRUE, survivals)]
  dead <- vcf[, c(TRUE, !survivals)]

  # $alive items correspond to a TRUE in the survivals vector
  result <- list("alive" = alive, "dead" = dead) 

  return(result)

}


#' A function to remove non-biallelic variants from a vcfR::VCF object
#'
#' @param vcf The original vcfR::VCF object to process
#' @param verbose Logical. If TRUE (default), report status of the process
#' along the execution.
#'
#' @return 
#' A vcfR::VCF object rid of its non-biallelic variants
#'
#' @export
#'
#' @examples
#' RmNonBiallelics(vcf)
#' RmNonBiallelics(vcf, verbose = FALSE)

RmNonBiallelics <- function(vcf, verbose = TRUE) {

  if (class(vcf) != "vcfR" && class(vcf) != "list") {
    stop("Parameter vcf must be a vcfR object")
  }

  if (verbose) {
    message("Omitting non-biallelic variants...")
  }

  if (typeof(vcf) == "list") {
    vcf <- lapply(vcf,
                  FUN = function(x){
                    return(x[vcfR::is.biallelic(x),])
                    }
            ) 
  } else {
    vcf <- vcf[vcfR::is.biallelic(vcf), ] #remove non-biallelics
  }

  return(vcf)

}


#' A function to remove indel variants from a vcfR::VCF object
#'
#' @param vcf The original vcfR::VCF object to process
#' @param verbose Logical. If TRUE (default), report status of the process
#' along the execution.
#'
#' @return 
#' A vcfR::VCF object rid of its indel variants
#'
#' @export
#'
#' @examples
#' RmIndels(vcf)
#' RmIndels(vcf, verbose = FALSE)

RmIndels <- function(vcf, verbose = TRUE) {

  result <- NULL

  if (class(vcf) != "vcfR" && class(vcf) != "list") {
    stop("Parameter vcf must be a vcfR object")
  }

  if (verbose) {
    message("Omitting indel variants...")
  }

  if (typeof(vcf) == "list") {
    result <- lapply(vcf, 
                  FUN = function(x){
                    x <- RmIndels(x, verbose=verbose)
                    return(x)
                  }
            )
  } else {
    result <- vcfR::extract.indels(vcf) #remove indels
  }

  return(result)

}