#library(DescTools)

#' @title Equivalence Testing for Categorical Variables
#' @description Testing for the presence of a negligible association between two categorical variables
#'
#' @param x first categorical variable
#' @param y second categorical variable
#' @param tab contingency table for the two predictor variables
#' @param eiU upper limit of equivalence interval
#' @param data data file containing the categorical variables
#' @param alpha nominal acceptable Type I error rate level
#'
#' @return returns a \code{list} containing each analysis and their respective statistics
#'   and decision
#' @export
#'
#' @examples
#' \dontrun{
#' #Example 1
#' x<-rbinom(10,1,.5)
#' y<-rbinom(10,1,.5)
#' neg.cat(x,y)
#' #Example 2
#' sex<-rep(c("m","f"),c(12,22))
#' haircol<-rep(c("bld","brn","bld","brn"),c(9,7,11,7))
#' d <- data.frame(sex,haircol)
#' tab<-table(sex,haircol)
#' neg.cat(tab=tab, alpha=.05)
#' neg.cat(x=sex,y=haircol)
#' neg.cat(x=sex,y=haircol,data=d)
#' }
neg.cat <- function (x = NULL, y = NULL,
      tab = NULL, eiU = .2, data = NULL, alpha = .05) {

  if (is.null(x) & is.null(y) & !is.null(data)) {
    tab <- deparse(substitute(tab))
    tab <- as.matrix(tab)
  }
  if (is.null(tab) & !is.null(data)) {
    d <- data.frame(x,y)
    d <- d[stats::complete.cases(d),]
    tab <- table(d$x, d$y)
  }
  if (is.null(tab) & is.null(data)) {
    tab <- table(x,y)
  }
    cv <- DescTools::CramerV(tab,conf.level=(1-2*alpha))
  propvar = cv[1]^2
  ifelse (cv[3] <= eiU,
        decis <- "The null hypothesis that the relationship between the categorical variables is substantial can be rejected",
        decis <- "The null hypothesis that the relationship between the categorical variables is substantial CANNOT be rejected")
  ret <- data.frame(cramv = cv[1],
                  propvar = propvar,
                  cil = cv[2],
                  ciu = cv[3],
                  eiU = eiU,
                  decis = decis,
                  alpha = alpha)
  class(ret) <- "neg.cat"
  return(ret)
}
#' @rdname neg.cat
#' @param x Data frame from neg.cat
#' @param ... extra arguments
#' @return
#' @export
#'
print.neg.cat <- function (x, ...) {
  cat("\n\n")
  cat("--------------------", "\n\n")
  cat("-- Negligible Effect Test of the Relationship --", "\n")
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

