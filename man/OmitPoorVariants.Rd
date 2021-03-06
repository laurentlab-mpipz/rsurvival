% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gt.R
\name{OmitPoorVariants}
\alias{OmitPoorVariants}
\title{A function to remove rows with a high ratio of missing data from a
dataframe. (e.g. variants from genotype)}
\usage{
OmitPoorVariants(gt, min.qual, dp = NULL, verbose = TRUE)
}
\arguments{
\item{gt}{A dataframe to process (e.g. a genotype dataframe)}

\item{min.qual}{The minimum present data ratio under which the row will be
deleted from the dataframe (from 0 to 1)}

\item{dp}{An optionnal dataframe on which the process of \code{gt} will be
executed too (e.g. a depth dataframe)}

\item{verbose}{Logical. If TRUE (default), report status of the process
along the execution.}
}
\value{
A dataframe or a list of dataframes. 
If the parameter \code{dp} is not set, return the \code{gt} dataframe
without its rows where the missing data ratio was too high. 
If the \code{dp} is set, return the modified \code{gt} as \code{list$gt};
\code{list$dp} is the modified \code{dp} dataframe with the same missing 
rows as \code{list$gt}.
}
\description{
A function to remove rows with a high ratio of missing data from a
dataframe. (e.g. variants from genotype)
}
\examples{
\dontrun{
OmitPoorVariants(genotype, 0.15)
OmitPoorVariants(genotype, 0.2, dp = depth)
OmitPoorVariants(genotype, 0.2, verbose = FALSE)
}
}
\seealso{
\code{\link{OmitPoorSamples}}
}
