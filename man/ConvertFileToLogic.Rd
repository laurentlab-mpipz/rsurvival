% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{ConvertFileToLogic}
\alias{ConvertFileToLogic}
\title{A function to turn a file into a logical vector. Useful for making a
survivals vector to split your data in two parts.}
\usage{
ConvertFileToLogic(file.path, type = "number")
}
\arguments{
\item{file.path}{The path to the VCF file you want to analyse}

\item{type}{The type of the content of your file, can be "number" (default)
or "character"}
}
\value{
A logical vector representing the content of your file. If \code{type} is
set to "number" (default), for each line of your file, a "0" will result in
a FALSE and any other number will result in a TRUE. If \code{type} is set to
"character", "TRUE" or "T" or "true" or "True" will return TRUE, "FALSE" or
"F" or "False" or "false" will return a FALSE.
}
\description{
A function to turn a file into a logical vector. Useful for making a
survivals vector to split your data in two parts.
}
\examples{
\dontrun{
ConvertFileToLogic("myFile.txt")
ConvertFileToLogic("myFile.txt", type = "character")
}
}
