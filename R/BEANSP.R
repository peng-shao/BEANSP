#' BEANSP: A package for computating the notorious bar statistic.
#'
#' Package mBEANSP computes age specific nest survival rates based on a Bayesian
#' hierarchical model.
#'
#' The package mBEANSP can be applied on either one of the two simple data
#' formats, interval format and date format. For more details of data format,
#' see example 1 and example 2 instruction. For either data format, first use
#' the function "nestconv" to convert a data to the standard format including
#' (ul,ur,zl,zr,y); then apply the function "nestsr" to generate the nest
#' survival summary.
#'
#' @references Cao, J., He, C., Suedkamp Wells, K.M., Millspaugh, J.J., and
#' Ryan, M.R. (2009).  Modeling age and nest-specific survival using a
#' hierarchical Bayesian approach. Biometrics, 65, 1052-1062.
#'
#' @docType package
#' @name BEANSP
NULL