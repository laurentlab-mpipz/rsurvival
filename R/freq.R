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
#' @param extrapolate.freq If TRUE (default), the relative frequencies of the 
#' signinficant data will be extrapolated. (i.e. freq.al.REF + freq.al.ALT = 1,
#' even if there are missing datas
#' @param totals If TRUE (default), returned vector will include totals.
#' @param min.freq.gt Optional. Minimal value of relative frequency for each
#' genotype in a variant of \code{gt}. Variants which does not fulfill this
#' condition will be omitted. Default is NULL (disabled).
#' @param min.freq.al Optional. Minimal value of relative frequency for each
#' allele in a variant of \code{gt}. Variants which does not fulfill this
#' condition will be omitted. Default is NULL (disabled).
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
#' \dontrun{
#' CalcFreqGt(genotype)
#' CalcFreqGt(genotype, genotypic = FALSE, allelic = TRUE)
#' CalcFreqGt(genotype, absolute = FALSE, percentage = TRUE)
#' CalcFreqGt(genotype, backup.path = "example.csv")
#' CalcFreqGt(genotype, totals= FALSE, min.freq.al = 0.15)
#' }

CalcFreqGt <- function(gt, genotypic = TRUE, allelic = FALSE, absolute = TRUE,
                        percentage = FALSE, extrapolate.freq = TRUE,
                        totals = TRUE, backup.path = NULL, min.freq.gt = NULL,
                        min.freq.al = NULL) {

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

  # if (!is.null(min.freq.gt) || !is.null(min.freq.al)) {
  #   # filter the variants without NA frequencies
  #   # this appends when a freq is lower than min.freq.al or min.freq.gt
  #   filter <- apply(result, MARGIN = 1,
  #                   FUN = function(x) {
  #                     return(!as.logical(sum(is.na(x))))
  #                   }
  #             )
  #   if (sum(!filter) != 0){
  #     result <- result[filter, ]
  #   }
  #   #result <- result[!as.logical(sum(is.na(result[, 1])))]
  # }

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
#' @param extrapolate.freq If TRUE (default), the relative frequencies of the 
#' signinficant data will be extrapolated. (i.e. freq.al.REF + freq.al.ALT = 1,
#' even if there are missing datas)
#' @param totals If TRUE (default), returned vector will include totals.
#' @param min.freq.gt Optional. Minimal value of relative frequency for each
#' genotype in a \code{variant}. If \code{variant} does not fulfill this
#' condition, return NULL. Default is NULL (disabled).
#' @param min.freq.al Optional. Minimal value of relative frequency for each
#' allele in a \code{variant}. If \code{variant} does not fulfill this
#' condition, return NULL. Default is NULL (disabled).
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
#' \dontrun{
#' CalcFreqVariant(variant)
#' CalcFreqVariant(variant, genotypic = FALSE, allelic = TRUE)
#' CalcFreqVariant(variant, absolute = FALSE, percentage = TRUE)
#' CalcFreqVariant(variant, totals= FALSE, min.freq.al = 0.15)
#' }

CalcFreqVariant <- function(variant, genotypic = TRUE, allelic = FALSE,
                            absolute = TRUE, percentage = FALSE,
                            extrapolate.freq = TRUE, totals = TRUE,
                            min.freq.gt = 0, min.freq.al = 0,
                            include.ids = FALSE) {

  # fact <- ConvertVariantToFactors(variant)
  
  ids  <- NA

  if (include.ids) {

    ids <- list("HOMOREF" = which(variant == "HOMOREF"),
                "HETERO"  = which(variant == "HETERO"),
                "HOMOALT" = which(variant == "HOMOALT"),
                "NA"      = which(variant == "NA"))

    if (!genotypic && !allelic) {
      return(ids)
    }

  }

  counts    <- as.vector(table(unlist(variant)))
  counts[5] <- sum(counts)

  if (genotypic || allelic) {

    if (!is.null(min.freq.al)) {
      if (!is.numeric(min.freq.al)) {
        stop("Parameter min.freq.al must be a number.")
      }
    } else {
      min.freq.al <- 0
    }

    if (!is.null(min.freq.gt)) {
      if (!is.numeric(min.freq.gt)) {
        stop("Parameter min.freq.gt must be a number.")
      }
    } else {
      min.freq.gt <- 0
    }

    freqs <- ShapeCountsCpp(counts, absolute, genotypic = genotypic[1],
                             allelic = allelic[1], percentage = percentage[1],
                             extrapolateFreq = extrapolate.freq[1],
                             totals = totals[1], minFreqGt = min.freq.gt[1],
                             minFreqAl = min.freq.al[1])
    freqs[freqs < 0] <- 0

  } else {

    freqs <- c()

  }

  if (include.ids) {
    result <- list("freqs" = freqs, "ids" = ids)
  } else {
    result <- freqs
  }

  return(result)

}

#' A function to find the ids of counts.gt columns in a frequency matrix.
#'
#' @param freq the absolute frequency matrix
#'
#' @return 
#' A list of integers containing ids of specific columns. Items names are :
#' "homo.ref", "hetero", "homo.alt", "missval", "total". The value are the 
#' corresponding column ids.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' FindIdsGtCounts(freq)
#' }

FindIdsGtCounts <- function(freq){

  # index of the frequencies of the genotypes
  homo.ref <- grep("^.*count\\.gt\\.HOMOREF$", names(freq))
  hetero   <- grep("^.*count\\.gt\\.HETERO$", names(freq))
  homo.alt <- grep("^.*count\\.gt\\.HOMOALT$", names(freq))
  missval  <- grep("^.*count\\.gt\\.MISSVAL$", names(freq))
  total    <- grep("^.*count\\.gt\\.TOTAL$", names(freq))

  result <- list("homo.ref" = homo.ref, "hetero" = hetero,
                 "homo.alt" = homo.alt, "missval" = missval, "total"= total)

  return(result)

}

#' A function to find the ids of freq.al columns in a frequency matrix.
#'
#' @param freq the realtive frequency matrix
#'
#' @return 
#' A list of integers containing ids of specific columns. Items names are :
#' "ref", "alt", "missval". The value are the corresponding column ids.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' FindIdsAlFeqs(freq)
#' }

FindIdsAlFreqs <- function(freq){

  # index of the frequencies of the genotypes
  ref <- grep("^.*freq\\.al\\.REF$", names(freq))
  alt <- grep("^.*freq\\.al\\.ALT$", names(freq))
  missval  <- grep("^.*freq\\.al\\.MISSVAL$", names(freq))

  result <- list("ref" = ref, "alt" = alt, "missval" = missval)

  return(result)

}



#' A function to find the ids of counts.gt columns in a frequency matrix.
#'
#' @param freq the absolute frequency matrix
#'
#' @return 
#' A list of integers containing ids of specific columns. Items names are :
#' "homo.ref", "hetero", "homo.alt", "missval", "total". The value are the 
#' corresponding column ids.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' FindIdsGtCounts(freq)
#' }

FindIdsGtCounts <- function(freq){

  # index of the frequencies of the genotypes
  homo.ref <- grep("^.*count\\.gt\\.HOMOREF$", names(freq))
  hetero   <- grep("^.*count\\.gt\\.HETERO$", names(freq))
  homo.alt <- grep("^.*count\\.gt\\.HOMOALT$", names(freq))
  missval  <- grep("^.*count\\.gt\\.MISSVAL$", names(freq))
  total    <- grep("^.*count\\.gt\\.TOTAL$", names(freq))

  result <- list("homo.ref" = homo.ref, "hetero" = hetero,
                 "homo.alt" = homo.alt, "missval" = missval, "total"= total)

  return(result)

}

#' A function to find the ids of freq.al columns in a frequency matrix.
#'
#' @param freq the realtive frequency matrix
#'
#' @return 
#' A list of integers containing ids of specific columns. Items names are :
#' "ref", "alt", "missval". The value are the corresponding column ids.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' FindIdsAlFeqs(freq)
#' }

FindIdsAlFreqs <- function(freq){

  # index of the frequencies of the genotypes
  ref <- grep("^.*freq\\.al\\.REF$", names(freq))
  alt <- grep("^.*freq\\.al\\.ALT$", names(freq))
  missval  <- grep("^.*freq\\.al\\.MISSVAL$", names(freq))

  result <- list("ref" = ref, "alt" = alt, "missval" = missval)

  return(result)

}

#' A function to classify samples by genotypes for a single variant, use regex.
#' NOT USED ANYMORE.
#'
#' @param variant a variant character data frame.
#'
#' @return 
#' A list of vectors which contain the ids of samples for each genotype.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  IdentifyGt(variant)
#' }

IdentifyGt <- function(variant){

  result <- list("HOMOREF" = grep("0.*0", variant),
                 "HETERO" = grep("0.*1|1.*0", variant),
                 "HOMOALT" = grep("1.*1",variant),
                 "MISSVAL" = which(is.na(variant)))

  return(result)

}
