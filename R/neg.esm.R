#' Test for Evaluating Substantial Mediation
#'
#' Function computes the equivalence testing method (total effect) for evaluating substantial mediation and Kenny method for full mediation.
#'
#' @aliases neg.esm
#' @param x predictor variable
#' @param y outcome variable
#' @param m mediator variable
#' @param alpha alpha level (default = .05)
#' @param minc minimum correlation between x and y (default is .15)
#' @param eil lower bound of equivalence interval in standardized units(default is -.15)
#' @param eiu upper bound of equivalence interval in standardized units (default is .15)
#' @param nboot number of bootstraps (default = 500L)
#' @param data optional data argument
#' @param seed optional argument to set seed
#' @param plot logical, plotting the results (default = TRUE)
#' @param saveplot saving plots (default = FALSE)
#' @param ... extra arguments
#'
#' @return returns a \code{list} containing each analysis and their respective statistics
#'   and decision
#'
#' @author Rob Cribbie \email{cribbie@@yorku.ca} and
#'   Nataly Beribisky \email{natalyb1@@yorku.ca}
#' @export neg.esm
#' @examples
#' \dontrun{
#' #equivalence test for substantial mediation
#' with an equivalence interval of -.15 to .15
#' x<-rnorm(200,sd=2)
#' m<-.5*x + rnorm(100)
#' y<-.5*m + rnorm(100)
#' neg.esm(x,y,m, eil = -.15, eiu = .15)
#'}


neg.esm<-function(x,y,m,alpha=.05,minc=.15,
                 eil=-.15,eiu=.15,nboot=500L,
                 data=NULL, plot=TRUE, saveplot=FALSE,
                 seed = NA) {
  if (is.null(data)) {
    if(!is.numeric(x)) stop('Variable x must be a numeric variable!')
    if(!is.numeric(m)) stop('Variably m must be a numeric variable!')
    if(!is.numeric(y)) stop('Variably y must be a numeric variable!')
    x <- scale(x)
    y <- scale(y)
    m <- scale(m)
    dat<-data.frame(x,y,m) # returns values with incomplete cases removed
    dat <- stats::na.omit(dat)

  }

  if (!is.null(data)) {
    x<-deparse(substitute(x))
    y<-deparse(substitute(y))
    m<-deparse(substitute(m))

    x<-as.numeric(data[[x]])
    y<-as.numeric(data[[y]])
    m<-as.numeric(data[[m]])

    x <- scale(x)
    y <- scale(y)
    m <- scale(m)

    dat <- data.frame(x,y,m)
    dat <- stats::na.omit(dat)

  }

  if(is.na(seed)){
    seed <- sample(.Random.seed[1], size = 1)
  } else {
    seed <- seed
  }
  set.seed(seed)

  m <- '
    y ~ c*x + b*m
    # mediator
    m ~ a*x
    # indirect effect (a*b)
    ab := a*b
    # total effect
    total := c + (a*b)
    y ~~ y
    m ~~ m
    x ~~ x
    '
  fit <- lavaan::lavaan(m, bootstrap=nboot,
                        se = "bootstrap", data = dat)
  pel<-data.frame(lavaan::parameterEstimates(fit, level=1-2*alpha,
                                             ci=TRUE))
  dir_eff <- pel[1,5]
  med_p<-pel$pvalue[pel$label=='ab']
  ab_par<-pel$est[pel$label=='ab']
  c_par<-pel$est[pel$label=='total']
  cil <-pel$ci.lower[pel$label=='c']
  ciu <- pel$ci.upper[pel$label=='c']
  abdivc_k <- ab_par/c_par
  corxy <- stats::cor(dat$x,dat$y)
  ifelse(pel$ci.lower[pel$label=='c']>eil &
           pel$ci.upper[pel$label=='c']<eiu &
           abs(stats::cor(dat$y,dat$x))>=minc,esm_dec<-"Substantial Mediation CAN be concluded",
         esm_dec<-"Substantial Mediation CANNOT be concluded")

  #Kenny Method for Full Mediation
  #"One rule of thumb is that if one wants to claim complete
  #mediation ab/c should be at least .80."
  ifelse(abs(ab_par/c_par) > .8 & abs(c_par)>.2,
         kenny_dec<-"Full Mediation CAN be concluded",
         kenny_dec<-"Full Mediation CANNOT be concluded")


  # Calculate Proportional Distance
  if (dir_eff > 0) {
    EIc <- eiu
  }
  else {
    EIc <- eil
  }

  PD <- dir_eff/abs(EIc)

  #Calculate the PD CI
  propdis<-numeric(nboot)

  for (i in 1:nboot) {
    #Creating a resampled dataset
    sample_d = dat[sample(1:nrow(dat), nrow(dat), replace = TRUE), ]

    #Running the regression on these data
    m <- '
    y ~ c*x + b*m
    # mediator
    m ~ a*x
    # indirect effect (a*b)
    ab := a*b
    # total effect
    total := c + (a*b)
    y ~~ y
    m ~~ m
    x ~~ x
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
  }

  ci_pd<-stats::quantile(propdis,probs=c(alpha/2,1-alpha/2))

  #### Summary #####

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
                    title1 = title1,
                    title2 = title2,
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

  cat("***",x$title1, "***\n\n")
  cat("Number of bootstrap iterations:", x$nboot, "(random seed =", x$seed, ")\n")
  cat("Indirect effect:", x$ab_par,"\n")
  cat("Correlation between x and y (must be greater in magnitude than ",x$minc,")",": ", x$corxy, sep="","\n")
  cat((1-2*x$alpha)*100, "% CI on Direct Effect: ", "(", x$cil,", ", x$ciu,")", sep="", "\n")
  cat("Equivalence Interval:","Lower =", x$eil, ",", "Upper =", x$eiu, "\n")
  cat("Decision from the ESM:", x$esm_dec,"\n\n")

  cat("**********************\n\n")
  cat("***",x$title2, "***\n\n")
  cat("ab/c (must be greater in magnitude than .80):", x$abdivc_k, "\n")
  cat("Correlation between x and y (must be greater in magnitude than .2): ", x$corxy, sep="","\n")
  cat("Decision from Kenny procedure:", x$kenny_dec, "\n\n")

  cat("**********************\n\n")
  cat("***", "Proportional Distance: The proportional distance from the effect of interest to the equivalence interval of the same sign","***","\n\n")
  cat("Proportional Distance (PD):", x$PD, "\n")
  cat(100*(1-x$alpha), "% CI for PD: ", "(",x$cilpd,", ",x$ciupd,")", "\n", sep="")
  cat("**********************\n")

  if (x$pl == TRUE) {
    neg.pd(effect=x$dir_eff, PD = x$PD, EIsign=x$EIc, PDcil=x$cilpd, PDciu=x$ciupd, cil=x$cil, ciu=x$ciu, Elevel=100*(1-2*x$alpha), Plevel=100*(1-x$alpha), save = x$saveplot)
  }


}


