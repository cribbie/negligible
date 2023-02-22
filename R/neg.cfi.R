
#' Equivalence Tests for CFI
#'
#' Function performs one of six equivalence tests for CFI fit index.
#'
#' @aliases neg.cfi
#' @param alpha desired alpha level (default = .05)
#' @param mod lavaan model object
#' @param eq.bound lower end of equivalence interval for comparison; must be .99, .95, .92 or .90 if modif.eq.bound = TRUE
#' @param modif.eq.bound should the lower end of the equivalence interval be modified (default = FALSE)
#' @param ci.method method used to calculate confidence interval; options are "yuan", "equiv" or "yhy.boot"; "yuan" corresponds to (1-alpha) percent CI, "equiv" corresponds to (1-2alpha) percent CI, "yhy.boot" corresponds to (1-2alpha) percent boot CI (default = "equiv")
#' @param nboot number of bootstrap samples if "yhy.boot" is selected as ci.method (default = 250L)
#' @param round number of digits to round equivalence bound and confidence interval bounds (default = 3)
#' @param nbootpd number of bootstrap samples by "yhy.boot" for pd function
#' @param plot logical, plotting the results (default = TRUE)
#' @param saveplot saving plots (default = FALSE)
#' @param ... extra arguments
#'
#' @return returns a \code{list} containing analysis and respective statistics
#'   and decision.
#' \itemize{
#'    \item \code{title1} The title of the CFI equivalence test. The appropriate title of the test will be displayed depending on the ci.method chosen and whether modif.eq.bound is TRUE or FALSE.
#'    \item \code{cfi_index} The CFI index.
#'    \item \code{ci.method} The method for confidence interval calculation.
#'    \item \code{cfi_eq} The lower end of the confidence interval for the CFI index.
#'    \item \code{eq.bound} The equivalence bound.
#'    \item \code{PD} Proportional distance (PD).
#'    \item \code{cilpd} Lower bound of the 1-alpha CI for the PD.
#'    \item \code{ciupd} Upper bound of the 1-alpha CI for the PD.
#'}
#' @export
#' @details
#'#'
#' The user specifies the lavaan fitted model object, the desired equivalence bound, and method of confidence interval computation. By default, the function does not modify the equivalence bounds according to Yuan et al. (2016). The user can also choose to instead run an equivalence test using a modified equivalence bound if the equivalence bound to be modified is .01, .05, .08, or .10. Alpha level can also be modified.
#'
#' For information on modified equivalence bounds see Yuan, K. H., Chan, W., Marcoulides, G. A., & Bentler, P. M. (2016). Assessing structural equation models by equivalence testing with adjusted fit indexes. Structural Equation Modeling: A Multidisciplinary Journal, 23(3), 319-330. doi: https://doi.org/10.1080/10705511.2015.1065414.
#'
#' The proportional distance quantifies the proportional distance from 0 to the nearest negligible effect (equivalence) interval (here, eiU). As values get farther from 0 the relationship becomes more substantial, with values greater than 1 indicating that the effect falls outside of the negligible effect (equivalence) interval.
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

neg.cfi <- function(mod,alpha = .05, eq.bound, modif.eq.bound = FALSE,
                    ci.method = "equiv", nbootpd = 50,
                    nboot = 250L, round = 3, plot = TRUE,
                    saveplot = TRUE) {

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

    n <- lavaan::lavInspect(fit1, what = "nobs")-1

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

  # plots
  EIc <- eq.bound
  PD <- (cfi_index-1)/abs(EIc-1)

  YHY.boot.for.pd <- lavaan::bootstrapLavaan(mod, R = nbootpd, type = "yuan", FUN = lavaan::fitMeasures,
                                             fit.measures = "cfi")
  ci_lpd <- stats::quantile((YHY.boot.for.pd-1)/abs(EIc-1),alpha/2)
  ci_upd <- stats::quantile((YHY.boot.for.pd-1)/abs(EIc-1),1-alpha/2)
  ci_pd <- c(ci_lpd,ci_upd)

  ret <- data.frame(cfi_eq = cfi_eq,
                    eq.bound = eq.bound,
                    title1 = title1,
                    cfi_index = cfi_index,
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
                    oe="CFI",
                    il = 1)
  class(ret) <- "neg.cfi"
  return(ret)

}

#' @rdname neg.cfi
#' @param x object of class \code{neg.cfi}
#' @export
#'

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


  cat("**********************\n\n")
  cat("***", "Proportional Distance: The proportional distance from the effect of interest to the equivalence interval of the same sign","***","\n\n")
  cat("Proportional Distance (PD):", x$PD, "\n")
  cat(100*(1-x$alpha), "% CI for PD: ", "(",x$cilpd,", ",x$ciupd,")", "\n", sep="")
  cat("**********************\n")
  if (x$pl == TRUE) {
    neg.pd(effect=x$cfi_index, PD = x$PD, eil=x$eq.bound, eiu=x$eq.bound, PDcil=x$cilpd, PDciu=x$ciupd, cil=x$cfi_eq, ciu=x$cfi_eq, Elevel=100*(1-2*x$alpha), Plevel=100*(1-x$alpha), save = x$saveplot, oe=x$oe)
  }
}





