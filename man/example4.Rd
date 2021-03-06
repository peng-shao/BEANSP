% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/example4.R
\docType{data}
\name{example4}
\alias{example4}
\title{National Dove Nest Survival Data - Date Format}
\format{A data frame with 217 observations on the following 7 variables.:
  \describe{
    \item{\code{id}}{Nest Id.}
    \item{\code{D1}}{the date of the first encounter.}
    \item{\code{D2}}{the date of the second-to-last visit (before the outcome).}
    \item{\code{D3}}{the date of the last visit (after the outcome)}.
    \item{\code{Fate}}{Nest fate; if "F", y=0; if "S", y=1.}
    \item{\code{Censor}}{The missing fate indicator; if 1, the nest fage is missing;
    if 0, the nest fate is not missing}
    \item{\code{x1}}{the first nest-specific covariate.}
  }}
\usage{
example4
}
\description{
The original data were collected by Geissler et al.(1987). This
example dataset only contains a fragment of whold dataset.
}
\references{
blank
}
\keyword{datasets}
