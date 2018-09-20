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
#' \dontrun{
#' ConvertFileToLogic("myFile.txt")
#' ConvertFileToLogic("myFile.txt", type = "character")
#' }


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
#' \dontrun{
#' vector.true  <- c(NA, 1:3) 
#' vector.false <- c(NA, 1:5) 
#' IsAboveNARatio(vector.true, 0.2)
#' }

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


#' A function to split a vector into two vectors according to a survival 
#' logical vector
#'
#' @param vect The original vector to split
#' @param survival A logical vector for splitting the data into two vectors.
#' The length of \code{survival} must match the length of \code{vector}
#' @param verbose Logical. If TRUE (default), report status of the process
#' along the execution.
#'
#' @return 
#' A list of vectors.
#' \code{list$alive} is the vector containing the samples which position
#' corresponds with a TRUE in the \code{survival} vector.
#' The other samples are contained in \code{list$dead} as an other vector
#'
#' @export
#'
#' @examples
#' \dontrun{
#' SplitGt(gt, c(TRUE, TRUE, FALSE, ...))
#' SplitGt(gt, survival.vector, verbose = FALSE)
#' }

SplitVect <- function(vect, survival, verbose = TRUE){

  if ((length(vect) != length(survival)) && verbose) {
    warning("Both parameter should have the same length")
  }

  alive <- vect[survival]
  dead  <- vect[!survival]

  result <- list("alive" = alive, "dead" = dead)
  return(result)

}


#' A function to slice data frames without loosing row names
#'
#' @param df A dataframe to slice
#' @param rows A integer vector. The ids of the rows to return.
#'
#' @return
#' A dataframe containing only the rows number \code{rows} from \code{df}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' SliceDfRows(df, 1)
#' SliceDfRows(df, 1:5)
#' SliceDfRows(df, c(1,3,5))
#' }

SliceDfRows <- function(df, rows){

  if (length(rows) > 1){
    result <- df[rows, ]
  } else {
    df <- data.frame(df)
    df <- tibble::rownames_to_column(df, var = "rowname")
    df <- dplyr::slice(df, rows)
    result <- tibble::column_to_rownames(df, var = "rowname")
  }

  return(result)

}


#' A function to slice data frames without loosing column names
#'
#' @param df A dataframe to slice
#' @param columns A integer vector. The ids of the columns to return.
#'
#' @return
#' A dataframe containing only the columns number \code{columns} from
#' \code{df}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' SliceDfColumns(df, 1)
#' SliceDfColumns(df, 1:5)
#' SliceDfColumns(df, c(1,3,5))
#' }

SliceDfColumns <- function(df, columns){

  if (length(columns) > 1){
    result <- df[, columns]
  } else {
    result <- t(SliceDfRows(t(df), columns))
  }

  return(result)

}