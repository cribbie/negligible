#' Equivalence Tests for SRMR
#'
#' Function performs one of four equivalence tests for SRMR fit index.
#'
#' @aliases neg.srmr
#' @param alpha desired alpha level (default = .05)
#' @param mod lavaan model object
#' @param eq.bound upper bound of equivalence interval for comparison; must be .05 or .10 if modif.eq.bound = TRUE
#' @param modif.eq.bound should the upper bound of the equivalence interval be modified? (default = FALSE)
#' @param ci.method method used to calculate confidence interval; options are "MO" or "yhy.boot"; "MO" corresponds to (1-2alpha) percent CI, "yhy.boot" corresponds to (1-2alpha) percent boot CI (default = "MO")
#' @param usrmr fit index around which equivalence test should be structured (usrmr = TRUE which is the default states that usrmr from Maydeu-Olivares, 2017 will be used, otherwise srmr from fitmeasures() output in lavaan will be used)
#' @param nboot number of bootstrap samples if "yhy.boot" is selected as ci.method (default = 250L)
#' @param round number of digits to round equivalence bound and confidence interval bounds (default = 3)
#' @param ... extra arguments
#' @return returns a \code{list} including the following:
#' \itemize{
#'    \item \code{title1} The title of the SRMR equivalence test. The appropriate title of the test will be displayed depending on the ci.method chosen whether usrmr and modif.eq.bound are TRUE or FALSE.
#'    \item \code{srmr_index} The SRMR index.
#'    \item \code{ci.method} The method for confidence interval calculation (direct computation or bootstrap).
#'    \item \code{srmr_ci} The upper bound of the 1-2*alpha confidence interval for the RMSEA index.
#'    \item \code{eq.bound} The equivalence bound.
#'    \item \code{PD} Proportional distance (PD).
#'    \item \code{cilpd} Lower bound of the 1-alpha CI for the PD.
#'    \item \code{ciupd} Upper bound of the 1-alpha CI for the PD.
#'}
#' @export
#' @details
#'
#' The user specifies the lavaan fitted model object, the desired equivalence bound, the method of confidence interval computation, and whether unbiased SRMR or original SRMR should be used. By default, the function does not modify the equivalence bounds. The user can also choose to instead run an equivalence test using a modified equivalence bound of .05 or .10 multiplied by the average communality of the observed indicators. Alpha level can also be modified.
#'
#' For information on unbiased SRMR and its confidence interval computation see Maydeu-Olivares, A. (2017). Maximum likelihood estimation of structural equation models for continuous data: Standard errors and goodness of fit. Structural Equation Modeling: A Multidisciplinary Journal, 24(3), 383-394. https://doi.org/10.1080/10705511.2016.1269606
#'

#' @author Rob Cribbie \email{cribbie@@yorku.ca} and
#'   Nataly Beribisky \email{natalyb1@@yorku.ca}
#'
#' @export neg.srmr
#' @examples
#' d <- lavaan::HolzingerSwineford1939
#' hs.mod <- 'visual =~ x1 + x2 + x3
#' textual =~ x4 + x5 + x6
#' speed =~ x7 + x8 + x9'
#' fit1 <- lavaan::cfa(hs.mod, data = d)
#' neg.srmr(mod=fit1,alpha=.05,eq.bound=.08,usrmr = TRUE)

neg.srmr <- function (mod, alpha = 0.05, eq.bound, modif.eq.bound = FALSE,
                      ci.method = "MO", usrmr = TRUE, nboot = 250L, round = 3
                      )
{

  # Titles
  if (ci.method == "MO" & usrmr == TRUE & modif.eq.bound == FALSE) {
    title1 <- "Equivalence Based Fit Test for Unbiased SRMR"
  }
  if (ci.method == "yhy.boot" & usrmr == TRUE & modif.eq.bound == FALSE) {
    title1 <- "Equivalence Based Fit Test for Unbiased SRMR with Bootstrapping"
  }
  if (ci.method == "yhy.boot" & usrmr == FALSE & modif.eq.bound == FALSE) {
    title1 <- "Equivalence Based Fit Test for Biased SRMR with Bootstrapping"
  }
  if (ci.method == "MO" & usrmr == TRUE & modif.eq.bound == TRUE) {
    title1 <- "Equivalence Based Fit Test for Unbiased SRMR and Modified Equivalence Interval"
  }
  if (ci.method == "yhy.boot" & usrmr == TRUE & modif.eq.bound == TRUE) {
    title1 <- "Equivalence Based Fit Test for Unbiased SRMR with Bootstrapping and Modified Equivalence Interval"
  }
  if (ci.method == "yhy.boot" & usrmr == FALSE & modif.eq.bound == TRUE) {
    title1 <- "Equivalence Based Fit Test for Biased SRMR with Bootstrapping and Modified Equivalence Interval"
  }

  # SRMR index

  # function to obtain both uSRMR and SRMR
  get_both_srmrs <- function(out) {
    resid <- lavaan::lavResiduals(out)
    uSRMR <-  resid$summary$cov[5] # <- unbiased SRMR
    SRMR <-  resid$summary$cov[1] # <- biased SRMR (original)
    srmrs <- c(uSRMR, SRMR)
    names(srmrs) <- c("uSRMR","SRMR")
    srmrs
  }

  # obtain index
  if (usrmr == TRUE) {
    srmr_index <-  get_both_srmrs(mod)[1]
  } else {
    srmr_index <- get_both_srmrs(mod)[2]
  }

  # Equivalence Bounds

  # Modified Equivalence Bound
  if (modif.eq.bound == TRUE) {
    if (eq.bound != 0.05 & eq.bound != 0.1)
      stop("The specified equivalence bound must either be .05 or .10 when using modified equivalence bounds.")


    # obtaining R2 for modified equivalence bound test
    obj <- lavaan::summary(mod, rsquare = TRUE)
    rs <- obj$pe$est[obj$pe$op== "r2" & obj$pe$rhs %in%
                       dimnames((lavaan::lavInspect(mod, "observed"))$cov)[[1]]]
    avg_rs <- mean(rs)

    if (eq.bound == 0.05) {
      eq.bound <- 0.05 * avg_rs
    }

    if (eq.bound == 0.1) {
      eq.bound <- .10 * avg_rs
    }
  }



  # Confidence intervals

  # Bootstrap option
  if (ci.method == "yhy.boot" & usrmr == TRUE) {

    YHY.boot <- lavaan::bootstrapLavaan(mod, R = nboot,
                                        type = "yuan", FUN = get_both_srmrs)
    ## Calculate 95% CI for CFI under the YHY bootstrap
    uSRMR_CIs_boot <- stats::quantile(YHY.boot[,1], c(0.050, .950))
    srmr_ci <-uSRMR_CIs_boot[2]
    srmr_avg_b <- mean(YHY.boot[,1])
  }

  if (ci.method == "yhy.boot" & usrmr == FALSE) {
    YHY.boot <- lavaan::bootstrapLavaan(mod, R = nboot,
                                        type = "yuan", FUN = get_both_srmrs)
    ## Calculate 95% CI for CFI under the YHY bootstrap
    SRMR_CIs_boot <- stats::quantile(YHY.boot[,2], c(0.050, .950))
    srmr_ci <-SRMR_CIs_boot[2]
    srmr_avg_b <-mean(YHY.boot[,2])

  }

  # Analytic option
  if (ci.method == "MO" & usrmr == TRUE) {
    resid <- lavaan::lavResiduals(mod)
    srmr_ci <- resid$summary$cov[8]
    srmr_avg_b <- NA
  }
  ifelse(srmr_ci < eq.bound, decision <- "REJECT HO: The null hypothesis that the population SRMR exceeds the equivalence bound can be rejected. There is evidence to support satisfactory fit, given the value of the equivalence bound.",
         decision <- "FAIL TO REJECT HO: The null hypothesis that the population SRMR exceeds the equivalence bound cannot be rejected.")

  ret <- data.frame(srmr_ci = srmr_ci, eq.bound = eq.bound,
                    title1 = title1, srmr_index = srmr_index, decision = decision,
                    alpha = alpha, ci.method = ci.method, modif.eq.bound = modif.eq.bound,
                    round = round, usrmr = usrmr, srmr_avg_b = srmr_avg_b)
  class(ret) <- "neg.srmr"
  return(ret)
}


#' @rdname neg.srmr
#' @param x object of class \code{neg.srmr}
#' @export
#'

print.neg.srmr <- function(x, ...) {

  cat("----",x$title1, "----\n\n")
  cat(ifelse(x$usrmr == TRUE, "uSRMR index:", "SRMR index:"), x$srmr_index,"\n")
  cat(ifelse(x$ci.method == "yhy.boot", ifelse(x$usrmr == TRUE, "Average of bootstrapped uSRMRs: ","Average of bootstrapped SRMRs: "),""))
  cat(ifelse(x$ci.method == "yhy.boot", x$srmr_avg_b, ""),sep = "",ifelse(x$ci.method == "yhy.boot","\n",""))
  cat("*************************************\n")
  cat("Confidence Interval Method Selected:", x$ci.method, "\n")
  cat("Upper bound of ",(1-2*x$alpha)*100,sep = "", "% CI for SRMR: ", round(x$srmr_ci, digits = x$round),"\n")
  cat("*************************************\n")
  cat("Modified Equivalence Bound:", ifelse(x$modif.eq.bound == TRUE, "yes","no"),"\n")
  cat("Equivalence Bound:", round(x$eq.bound, digits = x$round), "\n")
  cat("*************************************\n")
  cat("Test Decision (comparing confidence interval to equivalence bound):\n",x$decision, "\n")


}


