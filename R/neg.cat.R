#' @title Equivalence Testing for Categorical Variables
#' @description Testing for the presence of a negligible association between two categorical variables
#'
#' @param v1 first categorical variable
#' @param v2 second categorical variable
#' @param tab contingency table for the two predictor variables
#' @param eiU upper limit of equivalence interval
#' @param data optional data file containing the categorical variables
#' @param alpha nominal acceptable Type I error rate level
#' @param plot logical; should a plot be printed out with the effect and the proportional distance
#' @param save should the plot be saved to 'jpg' or 'png'
#' @param nbootpd number of bootstrap samples for calculating the CI for the proportional distance
#'
#' @return A \code{list} containing the following:
#' \itemize{
#'   \item \code{cramv} Cramer's V statistic
#'   \item \code{propvar} Proportion of variance explained (V^2)
#'   \item \code{cil} Lower bound of the confidence interval for Cramer's V
#'   \item \code{ciu} Upper bound of the confidence interval for Cramer's V
#'   \item \code{eiU} Upper bound of the negligible effect (equivalence) interval
#'   \item \code{decis} NHST decision
#'   \item \code{PD} Proportional distance
#'   \item \code{CI95L} Lower bound of the 1-alpha CI for the PD
#'   \item \code{CI95U} Upper bound of the 1-alpha CI for the PD
#'   \item \code{alpha} Nominal Type I error rate
#' }
#' @export
#' @details This function evaluates whether a negligible relationship exists among two categorical variables.
#'
#' The statistical test is based on the Cramer's V statistic; namely addressing the question of whether the upper limit of the confidence interval for Cramer's V falls below the upper bound of the negligible effect (equivalence) interval (eiU).
#'
#' If the upper bound of the CI for Cramer's V falls below eiU, we can reject Ho: The relationship is nonnegligible (V >= eiU).
#'
#' eiU is set to .2 by default, but should be set based on the context of the research. Since Cramer's V statistic is in a correlation metric, setting eiU is a matter of determining what correlation is the minimally meaningful effect size (MMES) given the context of the research.
#'
#' Users can input either the names of the categorical variables (v1, v2) or a frequency (contingency) table (tab).
#'
#' The proportional distance (V/eiU) estimates the proportional distance of the effect from 0 to eiU, and acts as an alternative effect size measure.
#'
#' The confidence interval for the proportional distance is computed via bootstrapping (percentile bootstrap).
#'
#' @examples
#' sex<-rep(c("m","f"),c(12,22))
#' haircol<-rep(c("bld","brn","bld","brn"),c(9,7,11,7))
#' d <- data.frame(sex,haircol)
#' tab<-table(sex,haircol)
#' neg.cat(tab=tab, alpha=.05, nbootpd=5)
#' neg.cat(v1=sex, v2=haircol, data=d, nbootpd=5)
neg.cat <- function (v1 = NULL, v2 = NULL,
                     tab = NULL, eiU = .2, data = NULL,
                     plot = TRUE, save = FALSE, nbootpd = 1000,
                     alpha = .05) {

  if (!is.null(tab) & is.null(data)) {
    countsToCases <- function(x, countcol = "Freq") {
      idx <- rep.int(seq_len(nrow(x)), x[[countcol]])
      x[[countcol]] <- NULL
      x[idx, ]
    }
    v1<-countsToCases(as.data.frame(tab))[,1]
    v2<-countsToCases(as.data.frame(tab))[,2]
    dat<-data.frame(v1,v2)
  }
  if (is.null(tab) & !is.null(data)) {
    v1<-deparse(substitute(v1))
    v2<-deparse(substitute(v2))
    v1<-factor(data[[v1]])
    v2<-factor(data[[v2]])
    dat <- data.frame(v1,v2)
    names(dat)<-c("v1","v2")
    dat <- stats::na.omit(dat)
    v1 <- dat[, 1]
    v2 <- dat[, 2]
    tab <- table(v1, v2)
  }
  if (is.null(tab) & is.null(data)) {
    tab <- table(v1,v2)
    dat <- data.frame(v1,v2)
    dat <- stats::na.omit(dat)
  }

  if (any(tab<=1)) {
    print("Do not trust results, frequencies in one or more categories are too low for the analysis. You might also receive an error related to the low frequencies.")
  }
  cv <- DescTools::CramerV(tab,conf.level=(1-2*alpha))
  cv <- round(cv, 3)
  propvar = round(cv[1]^2, 3)
  ifelse (cv[3] <= eiU,
          decis <- "The null hypothesis that the relationship between the categorical variables is substantial can be rejected. A negligible relationship among the variables is concluded. Be sure to interpret the magnitude (and precision) of the effect size.",
          decis <- "The null hypothesis that the relationship between the categorical variables is substantial CANNOT be rejected. There is insufficient evidence to conclude a negligible effect. Be sure to interpret the magnitude (and precision) of the effect size.")
  cv
  #### Plots ####
  # Calculate Proportional Distance
  PD <- cv[1]/eiU

  # Confidence Interval for Proportional Distance
  propd<-numeric(nbootpd)
  for (i in 1:nbootpd) {
    xx<-dplyr::sample_n(dat,size=nrow(dat),replace=TRUE)
    tabxx<-table(xx$v1,xx$v2)
    cvpd<-DescTools::CramerV(tabxx)
    propd[i]<-cvpd/eiU
  }
  CI95L<-round(stats::quantile(propd,.025,na.rm=TRUE),3)
  CI95U<-round(stats::quantile(propd,.975,na.rm=TRUE),3)
  ret <- data.frame(cramv = cv[1],
                    propvar = propvar,
                    cil = cv[2],
                    ciu = cv[3],
                    eiU = eiU,
                    decis = decis,
                    PD = PD,
                    CI95L = CI95L,
                    CI95U = CI95U,
                    pl = plot,
                    alpha = alpha,
                    oe = "Cramer's V",
                    save = save)
  class(ret) <- "neg.cat"
  return(ret)
}

#' @rdname neg.cat
#' @param x Data frame from neg.cat
#' @param ... extra arguments
#' @export
#'
print.neg.cat <- function (x, ...) {
  cat("\n\n")
  cat("********************", "\n\n")
  cat("** Negligible Effect Test of the Relationship **", "\n")
  cat("** Between Two Categorical Variables **", "\n\n")
  cat("********************", "\n\n")
  cat("Nominal Type I error rate (alpha):", x$alpha, "\n\n")
  cat("********************", "\n\n")
  cat("Cramer's V: ", x$cramv, "\n\n")
  cat(100*(1-2*x$alpha), "% CI for Cramer's V: ", "(",x$cil,", ",x$ciu,")", "\n\n", sep="")
  cat("*******************", "\n\n")
  cat("Proportion of Shared Variability: ", x$propvar, "\n\n")
  cat("*******************", "\n\n")
  cat("Upper Bound of the Equivalence Interval (Correlation Metric): ", x$eiU, "\n\n")
  cat("Upper Bound of the ", 100*(1-2*x$alpha), "% CI for Cramer's V: ", x$ciu, "\n\n", sep="")
  cat("NHST Decision:", "\n")
  cat(x$decis,"\n\n")
  cat("*******************", "\n\n")
  cat("Proportional Distance","\n\n")
  cat("Proportional Distance:", x$PD,"\n")
  cat("Confidence Interval for the Proportional Distance: (",x$CI95L, ",",x$CI95U,")","\n\n",sep="")
  cat("Note: Confidence Interval for the Proportional Distance may not be precise with small N","\n")
  cat("*******************", "\n\n")


  if (x$pl == TRUE) {
    neg.pd(effect=x$cramv, PD = x$PD, eil=x$eiU, eiu=x$eiU,
    PDcil=x$CI95L, PDciu=x$CI95U, cil=x$cil, ciu=x$ciu,
    Elevel=100*(1-2*x$alpha), Plevel=100*(1-x$alpha),
    save = x$save, oe=x$oe)
  }

}

