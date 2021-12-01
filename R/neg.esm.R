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
#' @param eil lower bound of equivalence interval(default is -.15)
#' @param eiu upper bound of equivalence interval (default is .15)
#' @param nboot number of boostraps (default = 50L)
#' @param data optional data argument
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
#' with an equivalenceinterval of -.15 to .15
#' x<-stats::rnorm(200,sd=2)
#' m<-.5*x + stats::rnorm(100)
#' y<-.5*m + stats::rnorm(100)
#' eq.esm(x,y,m, eil = -.15, eiu = .15)
#'}


neg.esm<-function(x,y,m,alpha=.05,minc=.15,
              eil=-.15,eiu=.15,nboot=50L,
              data=NULL, ...) {
  if (is.null(data)) {
    if(!is.numeric(x)) stop('Variable x must be a numeric variable!')
    if(!is.numeric(m)) stop('Variably m must be a numeric variable!')
    if(!is.numeric(y)) stop('Variably y must be a numeric variable!')
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

    dat<-data.frame(x,y,m)
    dat <- stats::na.omit(dat)

  }

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

  #### Summary #####

  title1 <- "Equivalence Testing Method for Substantial Mediation (ESM)"
  title2 <- "Kenny Method for Full Mediation"

  stats_esm <- c(minc, corxy, eil, eiu, nboot, cil, ciu) # resample stats
  names(stats_esm) <- c("Minimum Correlation Between x and y", "Equivalence Interval Lower Bound",
                       "Equivalence Interval Upper Bound",
                       "# of Bootstraps", "5th Percentile", "95th Percentile")
  stats_kenny <- c(abdivc_k, kenny_dec)
  pd_stats <- c( EIc, PD) # proportional distance stats

  ret <- data.frame(EIc = EIc,
                    minc = minc,
                    title1 = title1,
                    title2 = title2,
                    corxy = corxy,
                    eil = eil,
                    eiu = eiu,
                    nboot = nboot,
                    cil = cil,
                    ciu, ciu,
                    esm_dec = esm_dec,
                    PD = PD,
                    ab_par = ab_par,
                    alpha = alpha,
                    abdivc_k = abdivc_k,
                    kenny_dec = kenny_dec,
                    title1 = title1,
                    title2 = title2
  )
  class(ret) <- "neg.esm"
  return(ret)
}
#' @rdname neg.esm
#' @param x object of class \code{neg.esm}
#' @return
#' @export
#'

print.neg.esm <- function(x, ...) {

  cat("----",x$title1, "----\n\n")
  cat("Indirect effect:", x$ab_par,"\n")
  cat("Correlation between x and y (must be greater in magnitude than ",x$minc,")",": ", x$corxy, sep="","\n")
  cat((1-2*x$alpha)*100, "% CI on Direct Effect: ", "(", x$cil,", ", x$ciu,")", sep="", "\n")
  cat("Equivalence Interval:","Lower =", x$eil, ",", "Upper =", x$eiu, "\n")
  cat("Decision from the ESM:", x$esm_dec,"\n")
  cat("**********************\n\n")

  cat("----",x$title2, "----\n\n")
  cat("ab/c (must be greater in magnitude than .80):", x$abdivc_k, "\n")
  cat("Decision from Kenny procedure:", x$kenny_dec, "\n")

  cat("**********************\n")


}

