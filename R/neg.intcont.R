#' @title Negligible Interaction Test for Continuous Predictors
#' @description Testing for the presence of a negligible interaction between two continuous predictor variables
#'
#' @param outcome continuous outcome variable
#' @param pred1 first continuous predictor variable
#' @param pred2 second continuous predictor variable
#' @param eiL lower limit of the negligible effect (equivalence) interval
#' @param eiU upper limit of the negligible effect (equivalence) interval
#' @param standardized logical; should the solution be based on standardized variables (and eiL/eiU)
#' @param nbootpd number of bootstrap samples for the calculation of the CI for the proportional distance
#' @param data optional data file containing the categorical variables
#' @param alpha nominal acceptable Type I error rate level
#' @param plot logical; should a plot be printed out with the effect and the proportional distance
#' @param save logical; should the plot be saved
#'
#' @return A \code{list} containing the following:
#' \itemize{
#'   \item \code{intcoef} Interaction coefficient
#'   \item \code{intcil} Lower bound of the 1-alpha CI for the interaction coefficient
#'   \item \code{intciu} Upper bound of the 1-alpha CI for the interaction coefficient
#'   \item \code{eiL} Lower bound of the negligible effect (equivalence) interval
#'   \item \code{eiU} Upper bound of the negligible effect (equivalence) interval
#'   \item \code{sprs} Semi-partial correlation squared for the interaction term
#'   \item \code{PD} Proportional distance
#'   \item \code{CI95L} Lower bound of the 1-alpha CI for the PD
#'   \item \code{CI95U} Upper bound of the 1-alpha CI for the PD
#'   \item \code{alpha} Nominal Type I error rate
#' }
#' @export
#' @details This function evaluates whether the interaction between two continuous predictor variables is negligible. This can be important for deciding whether to remove an interaction term from a model or to evaluate a hypothesis related to negligible interaction.
#'
#' eiL/eiU represent the bounds of the negligible effect (equivalence) interval (i.e., the minimally meaningful effect size, MMES) and should be set based on the context of the research. When standardized = TRUE, Acock (2014) suggests that the MMES for correlations can also be applied to standardized effects - Acock, A. C. (2014). A Gentle Introduction to Stata (4th ed.). Texas: Stata Press.
#'
#' User can input the outcome variable and two predictor variable names directly (i.e., without a data statement), or can use the data statement to indicate the dataset in which the variables can be found.
#'
#' The advantage of this approach when standardized = TRUE and there are only two predictors is that the Delta method is adopted. However, for general cases researchers can also use the neg.reg function.
#'
#' The proportional distance (interaction coefficient/negligible effect bound) estimates the proportional distance of the effect from 0 to negligible effect bound, and acts as an alternative effect size measure.
#'
#' The confidence interval for the proportional distance is computed via bootstrapping (percentile bootstrap).
#'
#' @examples
#' d<-perfectionism
#' neg.intcont(outcome = mpshfpre.sop, pred1 = cesdpre.total, pred2 = atqpre.total, data = d,
#' eiL = -.25, eiU = .25, standardized = TRUE, nbootpd = 100)
neg.intcont <- function (outcome = NULL, pred1 = NULL,
         pred2 = NULL, eiL, eiU, standardized = TRUE,
         nbootpd = 1000, data, alpha = .05,
         plot = TRUE, save = FALSE) {
  if (!is.null(data)) {
    outcome <- deparse(substitute(outcome))
    pred1 <- deparse(substitute(pred1))
    pred2 <- deparse(substitute(pred2))
    outcome <- as.numeric(data[[outcome]])
    pred1 <- as.numeric(data[[pred1]])
    pred2 <- as.numeric(data[[pred2]])
  }
  if (standardized == FALSE) {
    dat<-data.frame(outcome, pred1, pred2)
    dat<-dat[stats::complete.cases(dat),]
    mod <- stats::lm(outcome ~ pred1*pred2, data=dat)
    invisible(utils::capture.output(sprs<-rockchalk::getDeltaRsquare(mod)))
    cis <- stats::confint(mod, level = 1-2*alpha)
    intcil <- cis["pred1:pred2",1]
    intciu <- cis["pred1:pred2",2]
    cisr <- stats::confint(mod, level = 1-alpha)
    intcilr <- cisr["pred1:pred2",1]
    intciur <- cisr["pred1:pred2",2]
    intcoef <- stats::coef(mod)["pred1:pred2"]
    ifelse(intcil >= eiL & intciu <= eiU,
           decis <- "The null hypothesis that the interaction is not negligible can be rejected",
           decis <- "The null hypothesis that the interaction is not negligible CANNOT be rejected")
    test<-"Standard TOST method is used since standardized = FALSE"

    # Calculate Proportional Distance
    ifelse(intcoef<0, eiPD<-eiL, eiPD<-eiU)
    PD <- intcoef/eiPD
        # confidence interval for Proportional distance
    propd<-numeric(nbootpd)
    for (i in 1:nbootpd) {
      xx<-dplyr::sample_n(dat,size=nrow(dat),replace=TRUE)
      modpd <- stats::lm(outcome ~ pred1*pred2, data=xx)
      ic_pd <- stats::coef(modpd)["pred1:pred2"]
      ifelse(ic_pd<0, eipd<-eiL, eipd<-eiU)
      propd[i]<-ic_pd/eipd
    }
    CI95L<-stats::quantile(propd,.025,na.rm=TRUE)
    CI95U<-stats::quantile(propd,.975,na.rm=TRUE)
  }
  if (standardized == TRUE) {
    outcome <- scale(outcome)
    pred1 <- scale(pred1)
    pred2 <- scale(pred2)
    int <- pred1*pred2
    dat <- data.frame(outcome,pred1,pred2,int)
    dat <- dat[stats::complete.cases(dat),]
    m <- stats::lm(outcome ~ pred1*pred2, data=dat)
    invisible(utils::capture.output(sprs <- rockchalk::getDeltaRsquare(m)))
    invisible(utils::capture.output(SEs<-fungible::seBeta(X = dat[,2:4], y = dat[,1],
           cov.x = stats::cov(dat[,2:4]),
           cov.xy = stats::cov(dat[,1:4])[1,2:4],
           var.y = stats::var(dat[,1]),
           Nobs = nrow(dat),
           alpha = alpha*2)))
    intcil<-SEs$CIs[3,1]
    intciu<-SEs$CIs[3,3]
    intcoef<-SEs$CIs[3,2]
    ifelse(intcil >= eiL & intciu <= eiU,
           decis <- "The null hypothesis that the interaction is not negligible can be rejected",
           decis <- "The null hypothesis that the interaction is not negligible CANNOT be rejected")
    invisible(utils::capture.output(SEsr<-fungible::seBeta(X = dat[,2:4], y = dat[,1],
        cov.x = stats::cov(dat[,2:4]),
        cov.xy = stats::cov(dat[,1:4])[1,2:4],
        var.y = stats::var(dat[,1]),
        Nobs = nrow(dat),
        alpha = alpha)))
    intcilr<-SEsr$CIs[3,1]
    intciur<-SEsr$CIs[3,3]
        test<-"-- The Delta Method is used for Calculating the SEs with Standardized Variables --"

    # Calculate Proportional Distance
    ifelse(intcoef<0, eiPD<-eiL, eiPD<-eiU)
    PD <- intcoef/eiPD
    # confidence interval for Proportional distance
    propd<-numeric(nbootpd)
    for (i in 1:nbootpd) {
      xx<-dplyr::sample_n(dat,size=nrow(dat),replace=TRUE)
      modpd <- stats::lm(outcome ~ pred1*pred2, data=xx)
      ic_pd <- stats::coef(modpd)["pred1:pred2"]
      ifelse(ic_pd<0, eipd<-eiL, eipd<-eiU)
      propd[i]<-ic_pd/eipd
    }
    CI95L<-stats::quantile(propd,alpha/2,na.rm=TRUE)
    CI95U<-stats::quantile(propd,1-alpha/2,na.rm=TRUE)
  }
  ret <- data.frame(intcoef = intcoef,
                    intcil = intcil,
                    intciu = intciu,
                    intcilr = intcilr,
                    intciur = intciur,
                    eiL = eiL,
                    eiU = eiU,
                    sprs = sprs,
                    decis = decis,
                    test = test,
                    PD = PD,
                    eiPD = eiPD,
                    CI95L = CI95L,
                    CI95U = CI95U,
                    pl = plot,
                    alpha = alpha,
                    save = save)
  class(ret) <- "neg.intcont"
  return(ret)
}


print.neg.intcont <- function(x, ...) {
  cat("-- Evaluating Negligible Interaction with --\n")
  cat("-- Continuous Predictor and Outcome Variables --\n\n")
  cat(x$test, "\n\n")
  cat("Interaction Coefficient: ", "\n", x$intcoef, "\n\n")
  cat(100*(1-x$alpha), "% CI on the Interaction Coefficient: ", "\n", "(",x$intcilr,",",x$intciur,")","\n\n",sep="")
  cat("Interaction Semi-Partial Correlation Squared", "\n")
  cat(x$sprs, "\n\n")
  cat("A one unit increase in the first predictor results in a",
      x$intcoef, "unit change in the slope of the second predictor (or vice versa)","\n\n")
  cat("Note that if standardized = TRUE, one unit = one sd", "\n\n")
  cat("**********************\n\n")
  cat("-- Negligible Effect Testing via the 1-2*alpha CI --\n\n")
  cat(100*(1-2*x$alpha), "% CI on the Interaction Coefficient: ", "\n", "(",x$intcil,",",x$intciu,")","\n\n",sep="")
  cat("Lower Bound of the Negligible Effect (Equivalence) Interval: ", "\n", x$eiL, "\n\n")
  cat("Upper Bound of the Negligible Effect (Equivalence) Interval: ", "\n", x$eiU, "\n\n")
  cat("NHST Decision: ", "\n", x$decis, "\n\n")
  cat("**********************\n\n")

  if (x$pl == TRUE) {
    neg.pd(effect=x$intcoef, PD = x$PD, EIsign=x$eiPD, PDcil=x$CI95L, PDciu=x$CI95U, cil=x$intcil, ciu=x$intciu, Elevel=100*(1-2*x$alpha), Plevel=100*(1-x$alpha), save = x$save)
  }
}


