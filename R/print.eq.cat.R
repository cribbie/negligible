#' @title Print function for the eq.cat function
#'
#' @param x Data frame from eq.cat
#' @param ... extra arguments
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' }
print.eq.cat <- function (x, ...) {
  cat("\n\n")
  cat("--------------------", "\n\n")
  cat("-- Equivalence Test of the Relationship --", "\n")
  cat("-- Between Two Categorical Variables --", "\n\n")
  cat("-------------------", "\n\n")
  cat("Nominal Type I error rate (alpha):", x$alpha, "\n\n")
  cat("-------------------", "\n\n")
  cat("Cramer's V: ", x$cramv, "\n\n")
  cat(100*(1-2*x$alpha), "% CI for Cramer's V: ", "(",x$cil,", ",x$ciu,")", "\n\n", sep="")
  cat("-------------------", "\n\n")
  cat("Proportion of Shared Variability: ", x$propvar, "\n\n")
  cat("-------------------", "\n\n")
  cat("Upper Bound of the Equivalence Interval (Correlation Metric): ", x$eiU, "\n\n")
  cat("Upper Bound of the ", 100*(1-2*x$alpha), "% CI for Cramer's V: ", x$ciu, "\n\n", sep="")
  cat("NHST Decision:", "\n")
  cat(x$decis,"\n\n")
  cat("-------------------", "\n\n")
}
