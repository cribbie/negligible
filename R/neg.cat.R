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
#' @param plot should a plot be printed out with the effect and the proportional distance
#' @param save should the plot be saved to 'jpg' or 'png'
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
      tab = NULL, eiU = .2, data = NULL,
      plot = TRUE, save = FALSE, alpha = .05) {

  if (!is.null(tab) & is.null(data)) {
    #tab <- deparse(substitute(tab))
    #tab <- as.matrix(tab)
    countsToCases <- function(x, countcol = "Freq") {
      idx <- rep.int(seq_len(nrow(x)), x[[countcol]])
      x[[countcol]] <- NULL
      x[idx, ]
    }
    x<-countsToCases(as.data.frame(tab))[,1]
    y<-countsToCases(as.data.frame(tab))[,2]
    dat<-data.frame(x,y)
  }
  if (is.null(tab) & !is.null(data)) {
    x<-deparse(substitute(x))
    y<-deparse(substitute(y))
    x<-factor(data[[x]])
    y<-factor(data[[y]])
    dat <- data.frame(x,y)
    names(dat)<-c("x","y")
    dat <- stats::na.omit(dat)
    x <- dat[, 1]
    y <- dat[, 2]
    tab <- table(x, y)
  }
  if (is.null(tab) & is.null(data)) {
    tab <- table(x,y)
    dat <- data.frame(x,y)
    dat <- stats::na.omit(dat)
  }

  if (any(tab<=1)) {
    print("Do not trust results, frequencies in one or more categories are too low for the analysis. You might also receive an error related to the low frequencies.")
  }
      cv <- DescTools::CramerV(tab,conf.level=(1-2*alpha))
  propvar = cv[1]^2
  ifelse (cv[3] <= eiU,
        decis <- "The null hypothesis that the relationship between the categorical variables is substantial can be rejected",
        decis <- "The null hypothesis that the relationship between the categorical variables is substantial CANNOT be rejected")

  #### Plots ####
  # Calculate Proportional Distance
  PD <- cv[1]/eiU

  # confidence interval for Proportional distance
  propd<-numeric(1000)
  for (i in 1:1000) {
    xx<-dplyr::sample_n(dat,size=nrow(dat),replace=TRUE)
    tabxx<-table(xx$x,xx$y)
    cvpd<-DescTools::CramerV(tabxx)
    propd[i]<-cvpd/eiU
  }
  CI95L<-stats::quantile(propd,.025,na.rm=TRUE)
  CI95U<-stats::quantile(propd,.975,na.rm=TRUE)
  #  statfun <- function(x, data) {
  #    tabl <- table(x,y)
  #    propdis <- DescTools::CramerV(tabl)[1]/eiU
  #    return(propdis)
  #  }
  #  npbs <- np.boot(x = 1:length(x), statistic = statfun, data = data.frame(x,y), level = c(0.95), method = c("perc"))

  #  CI95 <- npbs$perc
  #  CI95L<-npbs$perc[1]
  #  CI95U<-npbs$perc[2]

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
                  save = save)
  class(ret) <- "neg.cat"
  return(ret)
}

#' @rdname neg.cat
#' @param z Data frame from neg.cat
#' @param ... extra arguments
#' @return
#' @export
#'
print.neg.cat <- function (z, ...) {
  cat("\n\n")
  cat("********************", "\n\n")
  cat("** Negligible Effect Test of the Relationship **", "\n")
  cat("** Between Two Categorical Variables **", "\n\n")
  cat("********************", "\n\n")
  cat("Nominal Type I error rate (alpha):", z$alpha, "\n\n")
  cat("********************", "\n\n")
  cat("Cramer's V: ", z$cramv, "\n\n")
  cat(100*(1-2*z$alpha), "% CI for Cramer's V: ", "(",z$cil,", ",z$ciu,")", "\n\n", sep="")
  cat("*******************", "\n\n")
  cat("Proportion of Shared Variability: ", z$propvar, "\n\n")
  cat("*******************", "\n\n")
  cat("Upper Bound of the Equivalence Interval (Correlation Metric): ", z$eiU, "\n\n")
  cat("Upper Bound of the ", 100*(1-2*z$alpha), "% CI for Cramer's V: ", z$ciu, "\n\n", sep="")
  cat("NHST Decision:", "\n")
  cat(z$decis,"\n\n")
  cat("*******************", "\n\n")
  cat("Proportional Distance","\n\n")
  cat("Proportional Distance:", z$PD,"\n")
  cat("Confidence Interval for the Proportional Distance: (",z$CI95L, ",",z$CI95U,")","\n\n",sep="")
  cat("Note: Confidence Interval for the Proportional Distance may not be precise with small N","\n")
  cat("*******************", "\n\n")


  if (z$pl == TRUE) {
    neg.pd(effect=z$cramv, PD = z$PD, EIsign=z$eiU, PDcil=z$CI95L, PDciu=z$CI95U, cil=z$cil, ciu=z$ciu, Elevel=100*(1-2*z$alpha), Plevel=100*(1-z$alpha), save = z$save)
  }

}


