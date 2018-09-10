#' A function to turn a file into a logical vector. Useful for making a
#' survivals vector to split your data in two parts.
#'
#' @param file.path The path to the VCF file you want to analyse
#' @param type The type of the content of your file, can be "number" (default)
#' or "character"
#'
#' @return 
#' A logical vector representing the content of your file. If \code{type} is
#' set to "number" (default), for each line of your file, a "0" will result in
#' a FALSE and any other number will result in a TRUE. If \code{type} is set to
#' "character", "TRUE" or "T" or "true" or "True" will return TRUE, "FALSE" or
#' "F" or "False" or "false" will return a FALSE.
#'
#' @export
#'
#' @examples
#' ConvertFileToLogic("myFile.txt")
#' ConvertFileToLogic("myFile.txt", type = "character")


ConvertFileToLogic <- function(file.path, type = "number") {

  result <- NULL

  if (file.exists(file.path)) {
    if (type == "number") {
      result <- as.logical(scan(file.path, what = integer()))
    } else if (type == "character") {
      result <- as.logical(scan(file.path, what = character()))
    } else {
      stop("Parameter type is not recognised")
    }
  }
  else {
    stop("Parameter file.path leads to a non-existing file")
  }

  if (length(result) < 1) {
    warning('Result is empty')
  }

  if (sum(is.na(result))) {
    warning("Result contains NAs")
  }
  
  return(result)                                     

}


#' A function to check if the ratio of NA in a vector is over a threshold.
#'
#' @param vect A vector to analyse
#' @param min.ratio A threshold ratio of NA values in vect (from 0 to 1)
#'
#' @return 
#' Logical. TRUE if the ratio of NA in vect is at least minRatio, else FALSE.
#'
#' @export
#'
#' @examples
#' IsAboveNARatio(vector, 0.2)

IsAtLeastNARatio <- function(vect, min.ratio) {

  if (!(is.numeric(min.ratio))) {
    stop("Parameter min.ratio must be a number")
  } else if (length(min.ratio) != 1){
    stop("Parameter min.ratio must be of length 1")
  } else if(min.ratio < 0  || min.ratio > 1){
    warning("Parameter min.ratio should in range [0:1]")
  } 

  if (length(vect) < 1){
    result <- NA
  } else if ((sum(is.na(vect)) / length(vect)) >= min.ratio) { #if too much NAs
    result <- TRUE
  } else {
     result <- FALSE
  }

  return(result)

}