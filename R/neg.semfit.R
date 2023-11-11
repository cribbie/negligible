#' Equivalence Tests for Fit Indices
#'
#' Function performs equivalence tests for RMSEA, CFI, and SRMR.
#'
#' @aliases neg.semfit
#' @param alpha desired alpha level (default = .05)
#' @param mod lavaan model object
#' @param rmsea.eq.bound upper bound of the equivalence interval for RMSEA for comparison; must be .01, .05, .08, or .10 if rmsea.modif.eq.bound = TRUE
#' @param rmsea.modif.eq.bound should the upper bound of the equivalence interval for RMSEA be modified (default = FALSE)
#' @param rmsea.ci.method method used to calculate confidence interval for RMSEA; options are "not.close" or "yhy.boot"; "not.close" corresponds to (1-2alpha) percent CI, "yhy.boot" corresponds to (1-2alpha) percent boot CI (default = "not.close")
#' @param rmsea.nboot number of bootstrap samples if "yhy.boot" is selected as rmsea.ci.method (default = 250L)
#' @param cfi.eq.bound lower bound of equivalence interval for CFI for comparison; must be .99, .95, .92 or .90 if cfi.modif.eq.bound = TRUE
#' @param cfi.modif.eq.bound should the lower bound of the equivalence interval for CFI be modified (default = FALSE)
#' @param cfi.ci.method method used to calculate confidence interval for CFI; options are "yuan", "equiv" or "yhy.boot"; "yuan" corresponds to (1-alpha) percent CI, "equiv" corresponds to (1-2alpha) percent CI, "yhy.boot" corresponds to (1-2alpha) percent boot CI (default = "equiv")
#' @param cfi.nboot number of bootstrap samples if "yhy.boot" is selected as cfi.ci.method (default = 250L)
#' @param srmr.eq.bound upper bound of equivalence interval for SRMR for comparison; must be .05 or .10 if modif.eq.bound = TRUE
#' @param srmr.modif.eq.bound should the upper bound of the equivalence interval for SRMR be modified? (default = FALSE)
#' @param srmr.ci.method method used to calculate confidence interval for SRMR; options are "MO" or "yhy.boot"; "MO" corresponds to (1-2alpha) percent CI, "yhy.boot" corresponds to (1-2alpha) percent boot CI (default = "MO")
#' @param usrmr fit index around which equivalence test should be structured (usrmr = TRUE which is the default states that usrmr from Maydeu-Olivares, 2017 will be used, otherwise srmr from fitmeasures() output in lavaan will be used)
#' @param srmr.nboot number of bootstrap samples if "yhy.boot" is selected as srmr.ci.method (default = 250L)
#' @param round number of digits to round equivalence bound and confidence interval bounds (default = 3)
#'
#' @return returns a \code{list} containing analysis and respective statistics
#'   and decision.
#' \itemize{
#'    \item \code{title1} The appropriate title of the test will be displayed depending on the ci.method chosen and whether modif.eq.bound is TRUE or FALSE.
#'    \item \code{cfi_index} The CFI index.
#'    \item \code{ci.method} The method for confidence interval calculation.
#'    \item \code{cfi_eq} The lower end of the confidence interval for the CFI index.
#'    \item \code{eq.bound} The equivalence bound.
#'}
#' @export
#' @details
#'#'
#' The user specifies the lavaan fitted model object, the desired equivalence bound, and method of confidence interval computation for RMSEA, CFI, and SRMR. By default, the function does not modify the equivalence bounds according to Yuan et al. (2016) or according to Shi et al. (2018). The user can also choose to instead run an equivalence test using a modified equivalence bound if the equivalence bound to be modified is .01, .05, .08, or .10 for RMSEA,.99, .95, .92 or .90 for CFI, .05 or .10 for SRMR.
#' Alpha level can also be modified.
#'
#' For information on modified equivalence bounds for CFI and RMSEA see Yuan, K. H., Chan, W., Marcoulides, G. A., & Bentler, P. M. (2016). Assessing structural equation models by equivalence testing with adjusted fit indexes. Structural Equation Modeling: A Multidisciplinary Journal, 23(3), 319-330. doi: https://doi.org/10.1080/10705511.2015.1065414.
#' For information on uSRMR and modified cut-offs for SRMR see:
#' Maydeu-Olivares, A. (2017). Maximum likelihood estimation of structural equation models for continuous data: Standard errors and goodness of fit. Structural Equation Modeling: A Multidisciplinary Journal, 24(3), 383-394.
#' Shi, D., Maydeu-Olivares, A., & DiStefano, C. (2018). The relationship between the standardized root mean square residual and model misspecification in factor analysis models. Multivariate Behavioral Research, 53(5), 676-694.
#'
#'
#' @author Rob Cribbie \email{cribbie@@yorku.ca} and
#'   Nataly Beribisky \email{natalyb1@@yorku.ca}
#' @export neg.cfi
#' @examples
#' d <- lavaan::HolzingerSwineford1939
#' hs.mod <- 'visual =~ x1 + x2 + x3
#' textual =~ x4 + x5 + x6
#' speed =~ x7 + x8 + x9'
#' fit1 <- lavaan::cfa(hs.mod, data = d)
#' neg.cfi(mod = fit1, alpha = .05, eq.bound = .95,  modif.eq.bound = FALSE, ci.method = "equiv",
#' round = 3, plot = TRUE)
#'


neg.semfit <- function(mod, alpha = 0.05, round = 3,
                        rmsea.eq.bound = 0.05, rmsea.modif.eq.bound = FALSE, rmsea.ci.method = "not.close", rmsea.nboot = 250L,
                        cfi.eq.bound = 0.95, cfi.modif.eq.bound = FALSE, cfi.ci.method = "yhy.boot", cfi.nboot = 250L,
                        srmr.eq.bound = 0.08, srmr.modif.eq.bound = FALSE, srmr.ci.method = "MO", usrmr = TRUE, srmr.nboot = 250L)
{
# ----------SRMR----------------
  neg.srmr <- function (mod, alpha, eq.bound, modif.eq.bound,
                      ci.method, usrmr, nboot, round,
                      nbootpd)
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

#---------- RMSEA--------
neg.rmsea <- function(mod, alpha, eq.bound, modif.eq.bound,
                      ci.method, nboot, round) {

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


  ret <- data.frame(rmsea_eq = rmsea_eq,
                    eq.bound = eq.bound,
                    il = 0,
                    title1 = title1,
                    rmsea_index = rmsea_index,
                    decision = decision,
                    alpha = alpha,
                    ci.method = ci.method,
                    modif.eq.bound = modif.eq.bound,
                    round = round)
  class(ret) <- "neg.rmsea"
  return(ret)

}

# ----------CFI----------------
neg.cfi <- function(mod, alpha, eq.bound, modif.eq.bound,
                    ci.method, nbootpd,
                    nboot, round) {

  ncp_chi2=function(alpha, T_ml,df){
    z=stats::qnorm(1-alpha);
    z2=z*z; z3=z2*z; z4=z3*z; z5=z4*z;
    sig2=2*(2*T_ml-df+2);
    sig=sqrt(sig2); sig3=sig*sig2; sig4=sig2*sig2;sig5=sig4*sig;
    sig6=sig2*sig4;

    delta=T_ml-df+2+sig*
      (
        z+(z2-1)/sig-z/sig2 + 2*(df-1)*(z2-1)/(3*sig3)
        +( -(df-1)*(4*z3-z)/6+(df-2)*z/2 )/sig4
        +4*(df-1)*(3*z4+2*z2-11)/(15*sig5)
        +(
          -(df-1)*(96*z5+164*z3-767*z)/90-4*(df-1)*(df-2)*(2*z3-5*z)/9
          +(df-2)*z/2
        )/sig6
      );
    delta=max(delta,0);
    return(delta)
  }


  # model T statistic and df
  T_ml<- lavaan::fitmeasures(mod)[3]
  df <- lavaan::fitmeasures(mod)[4]

  # baseline T statistic and df
  T_mli <- lavaan::fitmeasures(mod)[6]
  df_i <- lavaan::fitmeasures(mod)[7]

  # title
  if (modif.eq.bound == FALSE & ci.method == "equiv") {
    title1 <- "EBF-CFI: Equivalence Based Fit Test for CFI"
  }

  if (modif.eq.bound == FALSE & ci.method == "yuan") {
    title1 <- "EBF alpha-CFI: Equivalence Based Fit Test for CFI with Modified CI"
  }

  if (modif.eq.bound == FALSE & ci.method == "yhy.boot") {
    title1 <- "EBFB-CFI: Equivalence Based Fit Test for CFI using YHY Bootstrap for CI"
  }

  if (modif.eq.bound == TRUE & ci.method == "equiv") {
    title1 <- "EBF-CFI-A: Equivalence Based Fit Test for CFI using Yuan et al. (2016) Modified Equivalence Bounds"
  }

  if (modif.eq.bound == TRUE & ci.method == "yuan") {
    title1 <- "EBF alpha-CFI-A: Equivalence Based Fit Test for CFI with Modified CI using Yuan et al. (2016) Modified Equivalence Bounds"
  }

  if (modif.eq.bound == TRUE & ci.method == "yhy.boot") {
    title1 <- "EBFB-CFI-A: Equivalence Based Fit Test for CFI using YHY Bootstrap for CI and Yuan et al. (2016) Modified Equivalence Bounds"
  }

  if (modif.eq.bound == TRUE) {
    if(eq.bound != .99 & eq.bound != .95 & eq.bound != .92 &  eq.bound != .9) stop("The specified equivalence bound must either be .99, .95, .92, or .90 when using modified equivalence bounds.")

    n <- lavaan::lavInspect(mod, what = "nobs")-1

    if (eq.bound == .99) {
      eq.bound <- 1-exp(4.67603-.50827*log(df)+.87087*(df^(1/5))-.59613*((df_i)^(1/5))-1.89602*log(n)
                        + .10190*((log(n))^2)+ .03729*log(df)*log(n))
    }

    if (eq.bound == .95) {
      eq.bound <-  1-exp(4.12132-.46285*log(df)+.52478*(df^(1/5))-.31832*((df_i)^(1/5))-1.74422*log(n)
                         +.13042*((log(n))^2)-.02360*(n^(1/2))+.04215*log(df)*log(n))
    }

    if (eq.bound == .92) {
      eq.bound <-  1-exp(6.31234-.41762*log(df)+.01554*((log(df))^2)-.00563*((log(df_i))^2)-1.30229*log(n)
                         +.19999*((log(n))^2)-2.17429*(n^(1/5))+.05342*log(df)*log(n)-.01520*log(df_i)*log(n))
    }

    if (eq.bound == .90) {
      eq.bound <- 1-exp(5.96633-.40425*log(df)+.01384*((log(df))^2)-.00411*((log(df_i))^2)-1.20242*log(n)
                        +.18763*((log(n))^2)-2.06704*(n^(1/5))+.05245*log(df)*log(n)-.01533*log(df_i)*log(n))
    }

  }
  # cfi index
  cfi_index <- lavaan::fitmeasures(mod)[9]

  if (ci.method == "equiv") {
    # ncp for model
    delta_t <-ncp_chi2(alpha = alpha, T_ml = T_ml, df = df)

    # ncp for baseline
    delta_it <- ncp_chi2(alpha = 1-alpha, T_ml = T_mli, df = df_i)


    #cfi equiv conf
    cfi_eq <- 1-max(delta_t,0, na.rm = T)/max(delta_t, delta_it,0, na.rm = T)
    cfi_eq
  }

  if (ci.method == "yuan") {
    # ncp for model
    delta_t <-ncp_chi2(alpha = alpha/2, T_ml = T_ml, df = df)

    # ncp for baseline
    delta_it <- ncp_chi2(alpha = 1-alpha/2, T_ml = T_mli, df = df_i)

    #cfi equiv conf
    cfi_eq <- 1-max(delta_t,0, na.rm = T)/max(delta_t, delta_it,0, na.rm = T)
    cfi_eq
  }

  if (ci.method == "yhy.boot") {
    YHY.boot <- lavaan::bootstrapLavaan(mod, R = nboot, type = "yuan", FUN = lavaan::fitMeasures,
                                        fit.measures = "cfi")
    cfi_eq <- stats::quantile(YHY.boot,alpha)
  }



  # decision
  ifelse(cfi_eq > eq.bound,
         decision <-"REJECT Ho: We have evidence to reject the hypothesis that the specified model is not substantially better fitting than the baseline model.",
         decision <-"FAIL TO REJECT HO: We fail to reject the hypothesis that the specified model is not substantially better fitting than the baseline model.")


  ret <- data.frame(cfi_eq = cfi_eq,
                    eq.bound = eq.bound,
                    title1 = title1,
                    cfi_index = cfi_index,
                    decision = decision,
                    alpha = alpha,
                    ci.method = ci.method,
                    modif.eq.bound = modif.eq.bound,
                    round = round
                    )
  class(ret) <- "neg.cfi"
  return(ret)

}

#------- Run them all ---------
rmsea.res <- neg.rmsea(mod = mod, alpha = alpha, round = round, eq.bound = rmsea.eq.bound,
          modif.eq.bound = rmsea.modif.eq.bound, ci.method = rmsea.ci.method,
          nboot = rmsea.nboot)
cfi.res <- neg.cfi(mod = mod, alpha = alpha, round = round, eq.bound = cfi.eq.bound,
        modif.eq.bound = cfi.modif.eq.bound, ci.method = cfi.ci.method,
        nboot = cfi.nboot)
srmr.res <- neg.srmr(mod = mod, alpha = alpha, round = round, eq.bound = srmr.eq.bound,
         modif.eq.bound = srmr.modif.eq.bound, ci.method = srmr.ci.method,
         usrmr = usrmr, nboot = srmr.nboot)


# neg.rmsea organization
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

}
# neg.cfi organization
print.neg.cfi <- function(x, ...) {

  cat("----",x$title1, "----\n\n")
  cat("CFI index:", x$cfi_index,"\n")
  cat("*************************************\n")
  cat("Confidence Interval Method Selected:", x$ci.method, "\n")
  cat("Lower end of ", ifelse(x$ci.method == "yuan",(1-x$alpha)*100,(1-2*x$alpha)*100),sep ="", "% CI for CFI: ", round(x$cfi_eq, digits = x$round), "\n")
  cat("*************************************\n")
  cat("Modified Equivalence Bound:", ifelse(x$modif.eq.bound == TRUE, "yes","no"),"\n")
  cat("Equivalence Bound:", round(x$eq.bound, digits = x$round), "\n")
  cat("*************************************\n")
  cat("Test Decision (comparing confidence interval to equivalence bound):",x$decision, "\n")

}


# neg.srmr organization

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
cat("** Equivalence/Negligible Effect Tests for Evaluating Model Fit **", "\n\n\n")

cat("* RMSEA-Based Test: *","\n\n")
print.neg.rmsea(rmsea.res)

cat("\n\n\n")
cat("* CFI-Based Test: *","\n")
print.neg.cfi(cfi.res)

cat("\n\n\n")
cat("* SRMR-Based Test: *","\n")
print.neg.srmr(srmr.res)
}

# Try it out

library(lavaan)
d <- lavaan::HolzingerSwineford1939
hs.mod <- 'visual =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed =~ x7 + x8 + x9'
fit1 <- lavaan::cfa(hs.mod, data = d)
neg.semfit(mod = fit1)
