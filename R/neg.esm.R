#' Test for Evaluating Substantial Mediation
#'
#' Function computes the equivalence testing method (total effect) for evaluating substantial mediation and Kenny method for full mediation.
#'
#' @aliases neg.esm
#' @param X predictor variable
#' @param Y outcome variable
#' @param M mediator variable
#' @param alpha alpha level (default = .05)
#' @param minc minimum correlation between x and Y (default is .15)
#' @param eil lower bound of equivalence interval in standardized units(default is -.15)
#' @param eiu upper bound of equivalence interval in standardized units (default is .15)
#' @param nboot number of bootstraps (default = 500L)
#' @param data optional data argument
#' @param seed optional argument to set seed
#' @param plot logical, plotting the results (default = TRUE)
#' @param saveplot saving plots (default = FALSE)
#' @param ... extra arguments
#'
#' @return A \code{list} including the following:
#' \itemize{
#'   \item \code{minc} Minimum correlation between X and Y for a valid negligible effect (equivalence) test
#'   \item \code{corxy} Sample correlation between the IV (X) and DV (Y)
#'   \item \code{dir_eff} Sample standardized direct effect between the IV (X) and DV (Y) after controlling for the mediator (M)
#'   \item \code{eiL} Lower bound of the negligible effect (equivalence) interval
#'   \item \code{eiU} Upper bound of the negligible effect (equivalence) interval
#'   \item \code{cil} Lower bound of the 1-2*alpha CI for the standardized direct effect of X on Y
#'   \item \code{ciu} Upper bound of the 1-2*alpha CI for the standardized direct effect of X on Y
#'   \item \code{PD} Proportional distance (PD)
#'   \item \code{cilpd} Lower bound of the 1-alpha CI for the PD
#'   \item \code{ciupd} Upper bound of the 1-alpha CI for the PD
#'   \item \code{ab_par} Standardized indirect effect
#'   \item \code{abdivc_k} Proportion mediated: Standardized indirect effect divided by the standardized total effect
#'   \item \code{alpha} Nominal Type I error rate
#' }
#' @export
#' @details This function evaluates whether a negligible direct effect of X on Y exists after controlling for the mediator. Another way to word this is that the indirect effect accounts for a substantial proportion of the variability in X-Y relationship. See Beribisky, Mara, and Cribbie (https://doi.org/10.20982/tqmp.16.4.p424)
#'
#' The user specifies the IV (X), DV (Y) and mediator (M). The user can also specify the alpha level, the lower/upper bound of the negligible effect interval (eiL, eiU), the number of bootstrap samples (nboot), as well as the minimum correlation between X and Y that is permitted for a valid test of substantial mediation.
#'
#' The variables X, Y and M can be specified as stand-alone, or a data argument can be used if the data reside in an R dataset.
#'
#' For the Kenny method see: https://davidakenny.net/cm/mediate.htm
#'
#' The proportional distance quantifies the proportional distance from 0 to the nearest negligible effect (equivalence) interval (eiL, eiU). As values get farther from 0 the relationship becomes more substantial, with values greater than 1 indicating that the effect falls outside of the negligible effect (equivalence) interval.
#'
#' Note that the number of bootstrap samples (nboot) are low for the example since the example has a time limit of 5 seconds to pass CRAN testing; we recommend running a much higher number of bootstrap samples for analyses.
#'
#' @author Rob Cribbie \email{cribbie@@yorku.ca} and
#'   Nataly Beribisky \email{natalyb1@@yorku.ca}
#' @export neg.esm
#' @examples
#' #equivalence test for substantial mediation
#' #with an equivalence interval of -.15 to .15
#' X<-rnorm(100,sd=2)
#' M<-.5*X + rnorm(100)
#' Y<-.5*M + rnorm(100)
#' neg.esm(X,Y,M, eil = -.15, eiu = .15, nboot = 5)
neg.esm<-function(X,Y,M,alpha=.05,minc=.15,
                  eil=-.15,eiu=.15,nboot=1000L,
                  data=NULL, plot=TRUE, saveplot=FALSE,
                  seed = NA) {
  if (is.null(data)) {
    if(!is.numeric(X)) stop('Variable X must be a numeric variable!')
    if(!is.numeric(M)) stop('Variable M must be a numeric variable!')
    if(!is.numeric(Y)) stop('Variable Y must be a numeric variable!')
    X <- scale(X)
    Y <- scale(Y)
    M <- scale(M)
    dat<-data.frame(X,Y,M) # returns values with incomplete cases removed
    dat <- stats::na.omit(dat)

  }

  if (!is.null(data)) {
    X<-deparse(substitute(X))
    Y<-deparse(substitute(Y))
    M<-deparse(substitute(M))

    X<-as.numeric(data[[X]])
    Y<-as.numeric(data[[Y]])
    M<-as.numeric(data[[M]])

    X <- scale(X)
    Y <- scale(Y)
    M <- scale(M)

    dat <- data.frame(X,Y,M)
    dat <- stats::na.omit(dat)

  }

  if(is.na(seed)){
    seed <- sample(.Random.seed[1], size = 1)
  } else {
    seed <- seed
  }
  set.seed(seed)

  m <- '
    Y ~ c*X + b*M
    # mediator
    M ~ a*X
    # indirect effect (a*b)
    ab := a*b
    # total effect
    total := c + (a*b)
    Y ~~ Y
    M ~~ M
    X ~~ X
    '
  fit <- lavaan::lavaan(m, bootstrap=nboot,
                        se = "bootstrap", data = dat)
  pel<-data.frame(lavaan::parameterEstimates(fit, level=1-2*alpha,
                                             ci=TRUE))
  dir_eff <- pel[1,5]
  med_p<-pel$pvalue[pel$label=='ab']
  ab_par<-pel$est[pel$label=='ab']
  ab_par<-round(ab_par, 3)
  c_par<-pel$est[pel$label=='total']
  cil <-pel$ci.lower[pel$label=='c']
  cil <- round(cil, 3)
  ciu <- pel$ci.upper[pel$label=='c']
  ciu <- round(ciu, 3)
  abdivc_k <- ab_par/c_par
  abdivc_k <- round(abdivc_k, 3)
  corxy <- stats::cor(dat$X,dat$Y)
  corxy <- round(corxy, 3)
  ifelse(pel$ci.lower[pel$label=='c']>eil &
           pel$ci.upper[pel$label=='c']<eiu &
           abs(stats::cor(dat$Y,dat$X))>=minc,esm_dec<-"The null hypothesis that the direct effect (difference between the total and indirect effect) is non-negligible can be rejected. Substantial Mediation CAN be concluded. Be sure to interpret the magnitude (and precision) of the effect size.",
         esm_dec<-"The null hypothesis that the direct effect (difference between the total and indirect effect) is non-negligible cannot be rejected. Substantial Mediation CANNOT be concluded. Be sure to interpret the magnitude (and precision) of the effect size.")

  #Kenny Method for Full Mediation
  #"One rule of thumb is that if one wants to claim complete
  #mediation ab/c should be at least .80."
  ifelse(abs(ab_par/c_par) > .8 & abs(c_par)>.2,
         kenny_dec<-"Full Mediation CAN be concluded",
         kenny_dec<-"Full Mediation CANNOT be concluded")

  #Effect Size
  prop_med<-MBESS::mediation(x=dat$X,mediator=dat$M,dv=dat$Y)$Effect.Sizes[7,1]
  prop_med<-round(prop_med, 3)
  csie<-MBESS::upsilon(x=dat$X,mediator=dat$M,dv=dat$Y,B=nboot)[1,1]
  csie<-round(csie,3)
  csie_lb<-MBESS::upsilon(x=dat$X,mediator=dat$M,dv=dat$Y,B=nboot)[1,2]
  csie_lb<-round(csie_lb, 3)
  csie_ub<-MBESS::upsilon(x=dat$X,mediator=dat$M,dv=dat$Y,B=nboot)[1,3]
  csie_ub<-round(csie_ub, 3)
  # Calculate Proportional Distance
  if (dir_eff > 0) {
    EIc <- eiu
  }
  else {
    EIc <- eil
  }

  PD <- dir_eff/abs(EIc)
  PD <- round(PD, 3)

  #Calculate the PD CI
  propdis<-numeric(nboot)
  propmed<-numeric(nboot)

  for (i in 1:nboot) {
    #Creating a resampled dataset
    sample_d = dat[sample(1:nrow(dat), nrow(dat), replace = TRUE), ]

    #Running the regression on these data
    m <- '
    Y ~ c*X + b*M
    # mediator
    M ~ a*X
    # indirect effect (a*b)
    ab := a*b
    # total effect
    total := c + (a*b)
    Y ~~ Y
    M ~~ M
    X ~~ X
    '
    fit <- lavaan::lavaan(m, data = sample_d)
    pel<-data.frame(lavaan::parameterEstimates(fit))

    if (pel[1,5] > 0) {
      EIc <- eiu
    }
    else {
      EIc <- eil
    }

    propdis[i] <- pel[1,5]/abs(EIc)
    propmed[i] <- MBESS::mediation(x=sample_d$X,mediator=sample_d$M,dv=sample_d$Y)$Effect.Sizes[7,1]
  }

  ci_pd<-stats::quantile(propdis,probs=c(alpha/2,1-alpha/2))
  ci_pd <- round(ci_pd, 3)
  ci_pm<-stats::quantile(propmed,probs=c(alpha/2,1-alpha/2))
  ci_pm<-round(ci_pm, 3)
  #### Summary #####
  title0 <- "Effect Sizes for the Indirect Effect"
  title1 <- "Equivalence Testing Method for Substantial Mediation (ESM)"
  title2 <- "Kenny Method for Full Mediation"

  stats_esm <- c(minc, corxy, eil, eiu, nboot, cil, ciu) # resample stats

  stats_kenny <- c(abdivc_k, kenny_dec)
  pd_stats <- c( EIc, PD) # proportional distance stats

  ret <- data.frame(EIc = EIc,
                    minc = minc,
                    title1 = title1,
                    title2 = title2,
                    corxy = corxy,
                    dir_eff = dir_eff,
                    eil = eil,
                    eiu = eiu,
                    nboot = nboot,
                    cil = cil,
                    ciu = ciu,
                    cilpd = ci_pd[1],
                    ciupd = ci_pd[2],
                    esm_dec = esm_dec,
                    PD = PD,
                    pl = plot,
                    saveplot = saveplot,
                    ab_par = ab_par,
                    alpha = alpha,
                    abdivc_k = abdivc_k,
                    kenny_dec = kenny_dec,
                    prop_med = prop_med,
                    prop_medlb = ci_pm[1],
                    prop_medub = ci_pm[2],
                    csie = csie,
                    csielb = csie_lb,
                    csieub = csie_ub,
                    title1 = title1,
                    title2 = title2,
                    oe="Direct Effect",
                    seed = seed
  )
  class(ret) <- "neg.esm"
  return(ret)
}
#' @rdname neg.esm
#' @param x object of class \code{neg.esm}
#' @export
#'

print.neg.esm <- function(x, ...) {

  cat("***",x$title0, "***\n\n")
  cat("Proportion Mediated:", x$prop_med, "\n")
  cat(100*(1-x$alpha), "% CI for Proportion Mediated: ", "(",x$prop_medlb,", ",x$prop_medub,")", "\n", sep="")
  cat("Completely Standardized Indirect Effect (CSIE):", x$csie, "\n")
  cat(100*(1-x$alpha), "% CI for CSIE: ", "(",x$csielb,", ",x$csieub,")", "\n\n", sep="")
  cat("***",x$title1, "***\n\n")
  cat("Number of Bootstrap Iterations: ", x$nboot, "(random seed = ", x$seed, ")\n", sep="")
  cat("Indirect Effect:", x$ab_par,"\n")
  cat("Correlation between X and Y (must be greater in magnitude than ",x$minc,")",": ", x$corxy, sep="","\n")
  cat((1-2*x$alpha)*100, "% CI on Direct Effect: ", "(", x$cil,", ", x$ciu,")", sep="", "\n")
  cat("Equivalence Interval: ","Lower = ", x$eil, "; ", "Upper = ", x$eiu, "\n\n", sep="")
  cat("Decision from the ESM:", x$esm_dec,"\n\n")

  cat("**********************\n\n")
  cat("***",x$title2, "***\n\n")
  cat("ab/c (must be greater in magnitude than .80):", x$abdivc_k, "\n")
  cat("Correlation between X and Y (must be greater in magnitude than .2): ", x$corxy, sep="","\n\n")
  cat("Decision from Kenny Procedure:", x$kenny_dec, "\n\n")

  cat("**********************\n\n")
  cat("***", "Proportional Distance (the proportional distance from the effect of interest to the equivalence interval of the same sign)","***","\n\n")
  cat("Proportional Distance (PD):", x$PD, "\n")
  cat(100*(1-x$alpha), "% CI for PD: ", "(",x$cilpd,", ",x$ciupd,")", "\n", sep="")
  cat("**********************\n")

  if (x$pl == TRUE) {
    neg.pd(effect=x$dir_eff, PD = x$PD, eil=x$eil, eiu=x$eiu, PDcil=x$cilpd, PDciu=x$ciupd, cil=x$cil, ciu=x$ciu, Elevel=100*(1-2*x$alpha), Plevel=100*(1-x$alpha), save = x$saveplot, oe=x$oe)
  }


}
