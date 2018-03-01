#' @name example2
#' @title Missouri Dickcissel Nest Survival Data - Date Format
#' @description The Missouri dickcissel dataset collected by Suedkamp Wells(2005)
#' is analyzed via BEANSP for illustration. For more detailed analysis, see Cao et.al(2009).
#' @format A data frame with 233 observations on the following 7 variables.:
#'   \describe{
#'     \item{\code{id}}{Nest Id.}
#'     \item{\code{date1}}{the date of the first encounter.}
#'     \item{\code{date2}}{the date of the second-to-last visit (before the outcome).}
#'     \item{\code{date3}}{the date of the last visit (after the outcome)}.
#'     \item{\code{fate}}{Nest fate; if "F", y=0; if "S", y=1.}
#'     \item{\code{cov1}}{the first nest-specific covariate.}
#'     \item{\code{cov2}}{the second nest-specific covariate.}
#'   }
#' @docType data
#' @keywords datasets
#' @references Cao, J., He, C., Suedkamp Wells, K.M., Millspaugh, J.J.,
#' and Ryan, M.R. (2009).  Modeling age and nest-specific survival using a
#' hierarchical Bayesian approach. Biometrics, 65, 1052-1062.
#'
"example2"
