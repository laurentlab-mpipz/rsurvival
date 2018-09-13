#' A function to count occurences of genotypes in a genotype matrix.
#'
#' @param gt The genotype matrix to analyse
#' @param genotypic If TRUE (default), include frequencies of genotypes in the
#' result
#' @param allelic If TRUE, include frequencies of alles in the result,
#' default set to FALSE
#' @param absolute If TRUE (default), occurences will be shown instead of
#' relative frequencies
#' @param percentage If TRUE and absolute is FALSE, show the relative
#' frequencies in percentage. Default is FALSE.
#' @param backup.path Optional. Path to the CSV backup file the result must be
#' saved to.
#' @param totals If TRUE (default), returned vector will include totals.
#' @param min.freq.gt Optional. Minimal value of relative frequency for each
#' genotype in a variant of \code{gt}. Variants which does not fulfill this
#' condition will be omitted. Default is NULL (disabled).
#' @param min.freq.al Optional. Minimal value of relative frequency for each
#' allele in a variant of \code{gt}. Variants which does not fulfill this
#' condition will be omitted. Default is 0.1, set it to NULL to disable this
#' feature.
#'
#' @return 
#' A data frame of frequencies.
#' Each row is a variant and each column is a calculated frequency. 
#' Column names depend on provided parameters.
#' Typical column names would be "count_gt_HomoREF" (absolute frequency of homo
#' REF genotype) or "freq_al_ALT" (relative frequency of ALT alleles).
#' Column names ending with "MissVal" concern samples with missing values for
#' this variant.
#' Column names ending with "Total" concern all samples for this variant,
#' including the ones with missing data.
#'
#' @seealso \code{\link{CalcFreqVariant}} which this function bind
#'
#' @export
#'
#' @examples
#' CalcFreqGt(genotype)
#' CalcFreqGt(genotype, genotypic = FALSE, allelic = TRUE)
#' CalcFreqGt(genotype, absolute = FALSE, percentage = TRUE)
#' CalcFreqGt(genotype, backup.path = "example.csv")
#' CalcFreqGt(genotype, totals= FALSE, min.freq.al = 0.15)

CalcFreqGt <- function(gt, genotypic = TRUE, allelic = FALSE, absolute = TRUE,
                        percentage = FALSE, totals = TRUE, backup.path = NULL,
                        min.freq.gt = NULL, min.freq.al = NULL) {

  result <- apply(gt, MARGIN = 1, 
                  FUN = function(x) {
                    CalcFreqVariant(x, genotypic = genotypic,
                                    allelic = allelic, absolute = absolute,
                                    percentage = percentage, totals = totals,
                                    min.freq.gt = min.freq.gt,
                                    min.freq.al = min.freq.al)
                  }
            )

  result <- t(result)

  if (!is.null(min.freq.gt) || !is.null(min.freq.al)) {
    # filter the variants without NA frequencies
    # this appends when a freq is lower than min.freq.al or min.freq.gt
    filter <- apply(result, MARGIN = 1,
                    FUN = function(x) {
                      return(!as.logical(sum(is.na(x))))
                    }
              )
    if (sum(!filter) != 0){
      result <- result[filter, ]
    }
    #result <- result[!as.logical(sum(is.na(result[, 1])))]
  }

  if(is.character(backup.path)){
    utils::write.csv(result, backup.path) # write to file
  }

  return(result)

}


#' A function to count occurences of genotypes in a variant vector.
#'
#' @param variant The genotype matrix to analyse
#' @param genotypic If TRUE (default), include frequencies of genotypes in the
#' result
#' @param allelic If TRUE, include frequencies of alles in the result,
#' default set to FALSE
#' @param absolute If TRUE (default), occurences will be shown instead of
#' relative frequencies
#' @param percentage If TRUE and absolute is FALSE, show the relative
#' frequencies in percentage
#' @param totals If TRUE (default), returned vector will include totals.
#' @param min.freq.gt Optional. Minimal value of relative frequency for each
#' genotype in a \code{variant}. If \code{variant} does not fulfill this
#' condition, return NULL. Default is NULL (disabled).
#' @param min.freq.al Optional. Minimal value of relative frequency for each
#' allele in a \code{variant}. If \code{variant} does not fulfill this
#' condition, return NULL. Default is 0.1, set it to NULL to disable this
#' feature.
#'
#' @return 
#' A vector of frequencies.
#' Each Item is a calculated frequency for this variant. 
#' Item names depend on provided parameters.
#' Typical item names would be "count_gt_HomoREF" (absolute frequency of homo
#' REF genotype) or "freq_al_ALT" (relative frequency of ALT alleles).
#' Item names ending with "MissVal" concern samples with missing values for
#' this variant.
#' Item names ending with "Total" concern all samples for this variant,
#' including the ones with missing data.
#'
#' @seealso \code{\link{CalcFreqGt}} which bind this function
#'
#' @export
#'
#' @examples
#' CalcFreqVariant(variant)
#' CalcFreqVariant(variant, genotypic = FALSE, allelic = TRUE)
#' CalcFreqVariant(variant, absolute = FALSE, percentage = TRUE)
#' CalcFreqVariant(variant, totals= FALSE, min.freq.al = 0.15)

CalcFreqVariant <- function(variant, genotypic = TRUE, allelic = FALSE,
                            absolute = TRUE, percentage = FALSE,
                            totals = TRUE, min.freq.gt = NULL,
                            min.freq.al = NULL) {

  result.al <- NULL
  result.gt <- NULL
  result    <- NULL 
  return.na <- FALSE

  # handle lists --------------------------------------------------------------

  if (typeof(variant) == "list") {
    result <- lapply(variant,
                     FUN = function(x){
                       CalcFreqVariant(x, genotypic = genotypic,
                                       allelic = allelic, absolute = absolute,
                                       percentage = percentage, totals = totals,
                                       min.freq.gt = min.freq.gt,
                                       min.freq.al = min.freq.al)
                     }
                   )
    return(result)
  }

  # actual counting using regex -----------------------------------------------

  counts <- c(sum(grepl("0.*0",variant)), sum(grepl("0.*1|1.*0",variant)),
               sum(grepl("1.*1",variant)), sum(is.na(variant)))

  # calculate frequencies of alleles ------------------------------------------

  counts.al <- c(counts[1] * 2 + counts[2], counts[3] * 2 + counts[2],
                2 * counts[4])
  total.al <- sum(counts.al)
  freqs.al <- counts.al / total.al

  # calculate frequencies of genotypes ----------------------------------------

  counts.gt <- c(counts[1], counts[2], counts[3], counts[4])
  total.gt  <- sum(counts.gt)
  freqs.gt  <- counts.gt / total.gt

  # check minimal frequencies -------------------------------------------------

  if (!is.null(min.freq.al)) {
    if (sum(freqs.al < min.freq.al) > 0) {
      return.na = TRUE
    }
  } else if (!is.null(min.freq.gt)) {
    if (sum(freqs.gt < min.freq.gt) > 0) {
      return.na = TRUE
    }
  }

  if (return.na) {
      counts.gt <- rep(NA, length(counts.gt))
      freqs.gt  <- rep(NA, length(freqs.gt))
      total.gt  <- rep(NA, length(total.gt))
      counts.al <- rep(NA, length(counts.al))
      freqs.al  <- rep(NA, length(freqs.al))
      total.al  <- rep(NA, length(total.al))
  }

  # build resulting frequencies -----------------------------------------------

  if (!absolute) {
    prefix    <- "freq"
    result.al <- freqs.al
    result.gt <- freqs.gt
    if (percentage) {
      prefix    <- "perc"
      result.al <- result.al * 100
      result.gt <- result.gt * 100
    }
  } else {
    prefix    <- "count"
    result.al <- counts.al
    result.gt <- counts.gt
  }

  # build titles for each frequencies returned --------------------------------

  names(result.al) <- paste(prefix,
                            c("al.REF", "al.ALT", "al.MISSVAL"), sep=".")
  names(result.gt) <- paste(prefix,
                            c("gt.HOMOREF", "gt.HETERO", "gt.HOMOALT",
                              "gt.MISSVAL"),
                            sep=".")
  names(total.al)  <- c("count.al.TOTAL")
  names(total.gt)  <- c("count.gt.TOTAL")

  # concatenate results if needed ---------------------------------------------

  if (totals) {
    result.al <- c(result.al, total.al)
    result.gt <- c(result.gt, total.gt)
  }

  if (genotypic && allelic) {
    result <- c(result.gt, result.al)
  } else if (genotypic) {
    result <- result.gt
  } else {
    result <- result.al
  }

  return(result)

}