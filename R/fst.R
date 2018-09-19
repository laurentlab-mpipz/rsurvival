#' A function to estimate the FST of a variant using the Weir & Cockerham
#' (1984) method
#' 
#' @param variant.pop1 A variant vector from a genotype matrix
#' @param vect.pop2 A variant vector from a genotype matrix
#'
#' @return
#' A numeric value. The FST estimation for the allelic records of this  two
#' variants.
#' 
#' @export

EstimateFST <- function(variant.pop1, variant.pop2){

  n1 <- ncol(variant.pop1)
  n2 <- ncol(variant.pop2)

  p1 <- data.frame(CalcFreqGt(variant.pop1, genotypic = FALSE,
                        allelic = TRUE, absolute = FALSE))$freq.al.REF
  p2 <- data.frame(CalcFreqGt(variant.pop2, genotypic = FALSE,
                        allelic = TRUE, absolute = FALSE))$freq.al.REF

  a <- (n1 * n2) / (n1 + n2)
  b <- 1 / (n1 + n2 - 2)
  c <- n1 * p1 * (1 - p1) + n2 * p2 * (1- p2)
  d <- (p1 - p2) ^ 2

  result <- 1 - (2 * a * b * c) / (a * d + (2 * a - 1) * b * c)

  return(result)

}