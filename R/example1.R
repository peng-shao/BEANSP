#' example1
#'
#' @title Missouri Dickcissel Nest Survival Data - Interval Format
#'
#' @description The Missouri dickcissel dataset collected by Suedkamp Wells(2005)
#' is analyzed via BEANSP for illustration. For more detailed analysis, see Cao et.al(2009).
#'
#' @format A data frame with 217 observations on the following 8 variables.:
#'   \describe{
#'     \item{\code{id}}{Nest id}
#'     \item{\code{fate}}{Nest fate; if 1, the nest fate is success;if 0, the
#'     nest fate is failure.}
#'     \item{\code{ul}}{The youngest possible age that the nest could have been
#'     when first encounterd.}
#'     \item{\code{ur}}{The oldest possible age that the nest could have been
#'     when first encounterd.}
#'     \item{\code{zl}}{The smallest possible number of time units from the
#'     first encounter date to the outcome date.}
#'     \item{\code{zr}}{The largest possible number of time units from the first
#'      encounter date to the outcome date.}
#'     \item{\code{cov1}}{The first nest-specific covariate.}
#'     \item{\code{cov2}}{The second nest-specific covariate.}
#'   }
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @references Cao, J., He, C., Suedkamp Wells, K.M., Millspaugh, J.J.,
#' and Ryan, M.R. (2009).  Modeling age and nest-specific survival using a
#' hierarchical Bayesian approach. Biometrics, 65, 1052-1062.
#'
"example1"
