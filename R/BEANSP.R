#' BEANSP: A package for computating the notorious bar statistic.
#'
#' @description Package BEANSP computes age specific nest survival rates based on a Bayesian
#' hierarchical model.
#'
#' @details  The package BEANSP can be applied on the case where either nest fate is
#' known or unknown. The datasets \code{example1} and \code{example2} contains missouri
#' dickcissel nest survival data, where the nest fate is known, and the datasets
#' \code{exmaple3} and \code{exmaple4} contains national mourning dove nest survival data,
#' where some nest fates are unknown.
#'
#' The package can be also applied on either one of the two simple data formats,
#' interval format and date format. For either data format, first use
#' the function "nestconv" to convert a data to the standard format including
#' \code{(zl,zr,y)} (for the uncensored data, the return also contains \code{(ul, ur)});
#' then apply the function "nestsr" to generate the nest survival summary.
#'
#' @seealso \code{\link{nestsr}}, \code{\link{nestconv}}
#' @author
#' Chong He, Yiqun Yang, Jing Cao, Peng Shao
#' @references Cao, J., He, C., Suedkamp Wells, K.M., Millspaugh, J.J., and
#' Ryan, M.R. (2009).  Modeling age and nest-specific survival using a hierarchical
#' Bayesian approach. Biometrics, 65, 1052-1062.
#'
#' @docType package
#' @name BEANSP
NULL
