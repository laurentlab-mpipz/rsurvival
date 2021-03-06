% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{SliceDfColumns}
\alias{SliceDfColumns}
\title{A function to slice data frames without loosing column names}
\usage{
SliceDfColumns(df, columns)
}
\arguments{
\item{df}{A dataframe to slice}

\item{columns}{A integer vector. The ids of the columns to return.}
}
\value{
A dataframe containing only the columns number \code{columns} from
\code{df}.
}
\description{
A function to slice data frames without loosing column names
}
\examples{
\dontrun{
SliceDfColumns(df, 1)
SliceDfColumns(df, 1:5)
SliceDfColumns(df, c(1,3,5))
}
}
