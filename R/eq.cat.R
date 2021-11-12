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
#' @param ... extra arguments
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' #Example 1
#' x<-rbinom(10,1,.5)
#' y<-rbinom(10,1,.5)
#' eq.cat(x,y)
#' #Example 2
#' sex<-rep(c("m","f"),c(12,22))
#' haircol<-rep(c("bld","brn","bld","brn"),c(9,7,11,7))
#' d <- data.frame(sex,haircol)
#' tab<-table(sex,haircol)
#' eq.cat(tab=tab, alpha=.05)
#' eq.cat(x=sex,y=haircol)
#' eq.cat(x=sex,y=haircol,data=d)
#' }
eq.cat <- function (x = NULL, y = NULL,
      tab = NULL, eiU = .2, data = NULL,
      alpha = .05, ...) {

  if (is.null(x) & is.null(y) & !is.null(data)) {
    tab <- deparse(substitute(tab))
    tab <- as.matrix(tab)
  }
  if (is.null(tab) & !is.null(data)) {
    d <- data.frame(x,y)
    d <- d[complete.cases(d),]
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
  class(ret) <- "eq.cat"
  return(ret)
}




