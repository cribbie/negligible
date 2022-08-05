#' @title Equivalence Tests for RMSEA
#'
#' @description Function performs one of four equivalence tests for RMSEA fit index.
#'
#' @aliases neg.rmsea
#'
#' @param alpha desired alpha level (default = .05)
#' @param mod lavaan model object
#' @param eq.bound upper end of equivalence interval for comparison; must be .01, .05, .08 or .10 if modif.eq.bound = TRUE
#' @param modif.eq.bound should the upper end of the equivalence interval be modified? (default = FALSE)
#' @param ci.method method used to calculate confidence interval; options are "not.close" or "yhy.boot"; "not.close" corresponds to 100(1-2alpha) percent CI, "yhy.boot" corresponds to 100(1-2alpha) percent boot CI (default = "not.close")
#' @param nboot number of bootstrap samples if "yhy.boot" is selected as ci.method (default = 250L)
#' @param nbootpd number of bootstrap samples by "yhy.boot" for pd function
#' @param round number of digits to round equivalence bound and confidence interval bounds (default = 3)
#' @param plot logical, plotting the results (default = TRUE)
#' @param saveplot saving plots (default = FALSE)
#'
#' @return returns a \code{list} including the following:
#' \itemize{
#'    \item \code{title1} The title of the RMSEA equivalence test. The appropriate title of the test will be displayed depending on the ci.method chosen and whether modif.eq.bound is TRUE or FALSE.
#'    \item \code{rmsea_index} The RMSEA index.
#'    \item \code{ci.method} The method for confidence interval calculation (direct computation or bootstrap).
#'    \item \code{rmsea_eq} The upper end of the 1-2*alpha confidence interval for the RMSEA index.
#'    \item \code{eq.bound} The equivalence bound.
#'    \item \code{PD} Proportional distance (PD).
#'    \item \code{cilpd} Lower bound of the 1-alpha CI for the PD.
#'    \item \code{ciupd} Upper bound of the 1-alpha CI for the PD.
#' }
#' @export
#' @details The user specifies the lavaan fitted model object, the desired equivalence bound, and method of confidence interal computation. By default, the function does not modify the equivalence bounds according to Yuan et al. (2016). The user can also choose to instead run an equivalence test using a modified equivalence bound if the equivalence bound to be modified is .01, .05, .08, or .10. Alpha level can also be modified.
#'
#' For information on modified equivalence bounds see Yuan, K. H., Chan, W., Marcoulides, G. A., & Bentler, P. M. (2016). Assessing structural equation models by equivalence testing with adjusted fit indexes. Structural Equation Modeling: A Multidisciplinary Journal, 23(3), 319-330. doi: https://doi.org/10.1080/10705511.2015.1065414.
#'
#' The proportional distance quantifies the proportional distance from 0 to the nearest negligible effect (equivalence) interval (here, eiU). As values get farther from 0 the relationship becomes more substantial, with values greater than 1 indicating that the effect falls outside of the negligible effect (equivalence) interval.
#'
#' @author Rob Cribbie \email{cribbie@@yorku.ca} and
#'   Nataly Beribisky \email{natalyb1@@yorku.ca}
#'
#' @export neg.rmsea
#' @examples
#' d <- lavaan::HolzingerSwineford1939
#' hs.mod <- 'visual =~ x1 + x2 + x3
#' textual =~ x4 + x5 + x6
#' speed =~ x7 + x8 + x9'
#' fit1 <- lavaan::cfa(hs.mod, data = d)
#' neg.rmsea(alpha = .05, mod = fit1, eq.bound = .05, ci.method = "not.close", modif.eq.bound = FALSE,
#' round = 5, nboot = 25, nbootpd = 25)

neg.rmsea <- function(mod, alpha = .05, eq.bound, modif.eq.bound = FALSE,
                     ci.method = "not.close", nbootpd = 50L,
                     nboot = 250L, round = 3, plot = TRUE, saveplot = FALSE) {

  # title
  if (ci.method == "not.close" & modif.eq.bound == FALSE) {
  title1 <- "EBF-RMSEA: Equivalence Based Fit Test for RMSEA; the Not-Close Fit Test for RMSEA by MacCallum et al. (1996)"
  }

  if (ci.method == "yhy.boot" & modif.eq.bound == FALSE) {
    title1 <- "EBFB-RMSEA: Equivalence Based Fit Test for RMSEA Using YHY Bootstrap"
  }

  if (ci.method == "not.close" & modif.eq.bound == TRUE) {
    title1 <- "EBF-RMSEA-A: Equivalence Based Fit Test for RMSEA Using Yuan et al. (2016) Modified Equivalence Bounds"
  }

  if (ci.method == "yhy.boot" & modif.eq.bound == TRUE) {
    title1 <- "EBFB-RMSEA-A: Equivalence Based Fit Test for RMSEA Using YHY Bootstrap and Yuan et al. (2016) Modified Equivalence Bounds"
  }

  # rmsea index
  rmsea_index <- lavaan::fitmeasures(mod)[23]

  if (ci.method == "not.close") {

  #rmsea equiv conf
    rmsea_eq <- lavaan::fitmeasures(mod)[25]
  }


  # rmsea equiv bound
  if (modif.eq.bound == TRUE) {
    if(eq.bound != .01 & eq.bound != .05 & eq.bound != .08 &  eq.bound != .1) stop("The specified equivalence bound must either be .01, .05, .08, or .10 when using modified equivalence bounds.")

    n <- lavaan::lavInspect(mod, what = "nobs")-1
    df <- lavaan::fitmeasures(mod)[4]

    if (eq.bound == .01) {
    eq.bound <- exp(
      1.34863-.51999*log(df)+.01925*log(df)*log(df)-.59811*log(n)+.00902*sqrt(n)+.01796*log(df)*log(n)
    )
    }


    if (eq.bound == .05) {
    eq.bound <- exp(
      2.06034-.62974*log(df)+.02512*log(df)*log(df)-.98388*log(n)
      +.05442*log(n)*log(n)-.00005188*n+.05260*log(df)*log(n)
    )
    }


    if (eq.bound == .08) {
    eq.bound <- exp(
      2.84129-.54809*log(df)+.02296*log(df)*log(df)-.76005*log(n)
      +.10229*log(n)*log(n)-1.11167*(n^.2)+.04845*log(df)*log(n)
    )
    }


    if (eq.bound == .10) {
    eq.bound <- exp(
      2.36352-.49440*log(df)+.02131*log(df)*log(df)-.64445*log(n)
      +.09043*log(n)*log(n)-1.01634*(n^.2)+.04422*log(df)*log(n)
    )
    }

  }

  if (ci.method == "yhy.boot") {
    YHY.boot <- lavaan::bootstrapLavaan(mod, R = nboot, type = "yuan", FUN = lavaan::fitMeasures,
                                        fit.measures = "rmsea")
    rmsea_eq <- stats::quantile(YHY.boot,1-alpha)
  }
  # decision
  ifelse(rmsea_eq < eq.bound,
         decision <-"REJECT HO: There is evidence to reject the hypothesis of not-close fit.",
         decision <-"FAIL TO REJECT HO: We fail to find evidence to reject the hypothesis of not-close fit.")

  # plots
  EIc <- eq.bound
  PD <- rmsea_index/abs(EIc)

  YHY.boot.for.pd <- lavaan::bootstrapLavaan(mod, R = nbootpd, type = "yuan", FUN = lavaan::fitMeasures,
                                      fit.measures = "rmsea")
  ci_upd <- stats::quantile(YHY.boot.for.pd/EIc,1-alpha/2)
  ci_pd <- c(0, ci_upd)

  ret <- data.frame(rmsea_eq = rmsea_eq,
                    eq.bound = eq.bound,
                    il = 0,
                    title1 = title1,
                    rmsea_index = rmsea_index,
                    decision = decision,
                    alpha = alpha,
                    ci.method = ci.method,
                    modif.eq.bound = modif.eq.bound,
                    round = round,
                    cilpd = ci_pd[1],
                    ciupd = ci_pd[2],
                    PD = PD,
                    pl = plot,
                    saveplot = saveplot,
                    oe="RMSEA")
  class(ret) <- "neg.rmsea"
  return(ret)

}

#' @rdname neg.rmsea
#' @param x object of class \code{neg.rmsea}
#' @param ... extra arguments
#' @export
#'

print.neg.rmsea <- function(x, ...) {

  cat("----",x$title1, "----\n\n")
  cat("RMSEA index:", x$rmsea_index,"\n")
  cat("*************************************\n")
  cat("Confidence Interval Method Selected:", x$ci.method, "\n")
  cat("Upper end of ",(1-2*x$alpha)*100,sep = "", "% CI for RMSEA: ", round(x$rmsea_eq, digits = x$round),"\n")
  cat("*************************************\n")
  cat("Modified Equivalence Bound:", ifelse(x$modif.eq.bound == TRUE, "yes","no"),"\n")
  cat("Equivalence Bound:", round(x$eq.bound, digits = x$round), "\n")
  cat("*************************************\n")
  cat("Test Decision (comparing confidence interval to equivalence bound):\n",x$decision, "\n")

  cat("**********************\n\n")
  cat("***", "Proportional Distance: The proportional distance from the effect of interest to the equivalence interval of the same sign","***","\n\n")
  cat("Proportional Distance (PD):", x$PD, "\n")
  cat(100*(1-x$alpha), "% CI for PD: ", "(",x$cilpd,", ",x$ciupd,")", "\n", sep="")
  cat("**********************\n")

  if (x$pl == TRUE) {
    neg.pd(effect=x$rmsea_index, PD = x$PD, eil=x$eq.bound, eiu=x$eq.bound, PDcil=x$cilpd, PDciu=x$ciupd, cil=x$il, ciu=x$rmsea_eq, Elevel=100*(1-2*x$alpha), Plevel=100*(1-x$alpha), save = x$saveplot, oe=x$oe)
  }
}


