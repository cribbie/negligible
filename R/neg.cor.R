
#' Test for Lack of Association between Two Continuous Normally Distributed Variables: Equivalence-based correlation tests
#'
#' Function performs an equivalence based test of lack of association with resampling.
#'
#' From Goertzen, J. R., & Cribbie, R. A. (2010). Detecting a lack of association. British Journal of Mathematical and Statistical Psychology, 63(3), 527–537
#'
#' @aliases neg.cor
#' @param v1 the first variable of interest
#' @param v2 the second variable of interest
#' @param eiL the lower bound of the equivalence interval, in terms of the magnitude of a correlation
#' @param eiU the upper bound of the equivalence interval, in terms of the magnitude of a correlation
#' @param data data frame where two variables (v1 and y) are contained - optional
#' @param alpha desired alpha level
#' @param na.rm logical; remove missing values?
#' @param plot whether or not to print graphics of the results (default = TRUE)
#' @param seed optional argument to set seed
#' @param saveplot saving plots (default = FALSE)
#' @param ... additional arguments to be passed
#'
#' @return A \code{list} including the following:
#' \itemize{
#'   \item \code{corxy} Sample correlation value
#'   \item \code{eiL} Lower bound of the negligible effect (equivalence) interval
#'   \item \code{eiU} Upper bound of the negligible effect (equivalence) interval
#'   \item \code{nresamples} Number of resamples for the bootstrapping procedure
#'   \item \code{q1} Lower bound of the confidence interval for the correlation
#'   \item \code{q2} Upper bound of the confidence interval for the correlation
#'   \item \code{PD} Proportional distance
#'   \item \code{CIPDL} Lower bound of the 1-alpha CI for the PD
#'   \item \code{CIPDU} Upper bound of the 1-alpha CI for the PD
#'   \item \code{alpha} Nominal Type I error rate
#' }
#' @export
#' @details This function evaluates whether a negligible relationship exists among two continuous variables.
#'
#' The statistical test is based on a bootstrap-generated 1-2*alpha CI for the correlation; in other words, does the 1-2*alpha CI for the falls completely within the negligible effect (equivalence) interval.
#' 
#' The user needs to specify the lower and upper bounds of the negligible effect (equivalence) interval (eiL,eiU). Since we working in a correlation magnitude, setting these bounds requires estimating the minimally meaningful effect size (MMES); in this case, the minimally meaningful correlation (e.g., eiL = - .3, eiU = .3).
#' 
#' The 'plot' argument, if TRUE, will generate a plot of the observed effect (correlation) with the associated 1-2*alpha CI, along with a plot of the PD and the associated 1-alpha CI.
#' 
#' @author Rob Cribbie \email{cribbie@@yorku.ca}
#'   Phil Chalmers \email{rphilip.chalmers@@gmail.com} and
#'   Nataly Beribisky \email{natalyb1@@yorku.ca}
#' @export neg.cor
#' @examples
#' #Negligible correlation test between v1 and v2
#' #with an interval of ei=(-.2.2)
#' v1 <- rnorm(50)
#' v2 <- rnorm(50)
#' plot(v1, v2)
#' cor(v1, v2)
#' neg.cor(v1 = v1, v2 = v2, eiU = .2, eiL = -.2)
neg.cor <- function(v1, v2, eiU, eiL, alpha = 0.05, na.rm = TRUE,
                    plot = TRUE, data=NULL, saveplot=FALSE, seed = NA,...) {
  if (is.null(data)) {
    if(!is.numeric(v1)) stop('Variable v1 must be a numeric variable!')
    if(!is.numeric(v2)) stop('Variable v2 must be a numeric variable!')
    dat <- data.frame(v1,v2) # returns values with incomplete cases removed
    dat <- stats::na.omit(dat)
    var1 <- v1 <- dat[, 1] # takes first column (corresponding to first var)
    var2 <- v2 <- dat[, 2] # takes second column (corresponding to second var)
  }

  if (!is.null(data)) {
    v1<-deparse(substitute(v1))
    v2<-deparse(substitute(v2))
    v1<-as.numeric(data[[v1]])
    v2<-as.numeric(data[[v2]])
    dat<- data.frame(v1,v2)

    dat <- stats::na.omit(dat)
    var1 <- v1 <- dat[, 1] # takes first column (corresponding to first var)
    var2 <- v2 <- dat[, 2] # takes second column (corresponding to second var)
  }

  if(is.na(seed)){
    seed <- sample(.Random.seed[1], size = 1)
  } else {
    seed <- seed
  }
  set.seed(seed)

  corxy <- stats::cor(var1, var2) # get correlation between two variables
  n <- length(var1) # number of observations
  nresamples <- 10000
  #### Run the resampling version of the two t-test
  #### procedure for equivalence #####
  resamp <- function(x, m = 10000, theta, conf.level = 1-alpha,  # m is default to 10000 and conf level default to .95
                     ...) {
    n <- length(x) # n is number of observations, see below.
    Data <- matrix(sample(x, size = n * m, replace = T), # sample from x, with replacement, choose n*1000 items
                   nrow = m)
    thetastar <- apply(Data, 1, theta, ...) # apply theta function to rows of the matrix of resampled data
    M <- mean(thetastar) # mean
    S <- stats::sd(thetastar) # sd
    alpha <- 1 - conf.level #alpha
    CI <- stats::quantile(thetastar, c(alpha/2, 1 -
                                  alpha/2))
    return(list(ThetaStar = thetastar, Mean.ThetaStar = M,
                S.E.ThetaStar = S, Percentile.CI = CI)) # return list with thetastar, mean, sd, and CI
  }
  matr <- cbind(var1, var2)
  mat <- as.matrix(matr)
  theta <- function(x, mat) { # this is the function that is applied to the matrix, it takes correlation between rows of x which are later specified from 1:n
    stats::cor(mat[x, 1], mat[x, 2])
  }
  results <- resamp(x = 1:n, m = nresamples, theta = theta,
                    mat = mat)
  q1 <- stats::quantile(results$ThetaStar, alpha)
  q2 <- stats::quantile(results$ThetaStar, 1-alpha)
  q1negei <- q1 + eiL # first quantile plus lower ei
  q2negei <- q2 + eiL # second quantile plus lower ei
  q1posei <- q1 + eiU # first quantile plus upper ei
  q2posei <- q2 + eiU # second quantile plus upper ei
  ifelse(q2negei < 0 & q1posei > 0, decis_rs <- "The null hypothesis that the correlation between v1 and v2 falls outside of the equivalence interval can be rejected. A negligible association CAN be concluded. Be sure to interpret the magnitude (and precision) of the effect size.",
         decis_rs <- "The null hypothesis that the correlation between v1 and v2 falls outside of the equivalence interval cannot be rejected. A negligible association CANNOT be concluded. Be sure to interpret the magnitude (and precision) of the effect size.")
  #### Plots ####
  # Calculate Proportional Distance
  if (corxy > 0) {
    EIc <- eiU
  }
  else {
    EIc <- eiL
  }

  PD <- corxy/abs(EIc)

  # confidence interval for Proportional distance

  statfun <- function(x, data) {
    sample_cor <- stats::cor(data[x,1], data[x,2])
    propdis <- sample_cor/abs(EIc)
    return(propdis)
  }
  npbs <- nptest::np.boot(x = 1:length(v1), statistic = statfun, data = data.frame(v1,v2), level = c(1-alpha), method = c("perc"))

  CIPD <- npbs$perc
  CIPDL<-npbs$perc[1]
  CIPDU<-npbs$perc[2]




  #### Summary #####

  title <- "Equivalence Based Test of Lack of Association with Resampling"

  stats_rs <- c(corxy, eiL, eiU, nresamples, q1, q2) # resample stats
  names(stats_rs) <- c("Pearson r", "Equivalence Interval Lower Bound",
                       "Equivalence Interval Upper Bound",
                       "# of Resamples", "5th Percentile", "95th Percentile")
  pd_stats <- c( EIc, CIPDL, CIPDU) # proportional distance stats

  ret <- data.frame(EIc = EIc,
                    seed = seed,
                    alpha = alpha,
                    saveplot = saveplot,
                    pl = plot,
                    title = title,
                    corxy = corxy,
                    eiL = eiL,
                    eiU = eiU,
                    nresamples = nresamples,
                    q1 = q1,
                    q2 = q2,
                    decis_rs = decis_rs,
                    PD = PD,
                    oe="Correlation (r)",
                    CIPDL = CIPDL,
                    CIPDU = CIPDU)

  class(ret) <- "neg.cor"
  return(ret)

}

#' @rdname neg.cor
#' @param x object of class \code{neg.cor}
#' @export
#'
print.neg.cor <- function(x, ...) {


  cat("----",x$title, "----\n\n")
  cat("Random Seed =", x$seed, "\n")
  cat("Pearson's r:", x$corxy, "\n\n")

  cat("Equivalence Interval:","Lower =", x$eiL, ",", "Upper =", x$eiU, "\n\n")

  cat("# of Resamples:", x$nresamples,"\n")
  cat("Bootstrapped ", (1-2*x$alpha)*100, "% ", "Confidence Interval (",x$alpha*100,"th Percentile, ", (1-x$alpha)*100,"th Percentile", "):", "(", x$q1, ",", x$q2, ")", "\n\n", sep="")
  cat(x$decis_rs, "\n\n")
  cat("**********************\n")
  cat("Proportional Distance (PD):", x$PD, "\n")
  cat(1-x$alpha,"% CI for PD: (", x$CIPDL,",", x$CIPDU, ") \n")
  cat("**********************\n")


  if (x$pl == TRUE) {
    neg.pd(effect=x$corxy, PD = x$PD, eil=x$eiL, eiu=x$eiU, PDcil=x$CIPDL, PDciu=x$CIPDU, cil=x$q1, ciu=x$q2, Elevel=100*(1-2*x$alpha), Plevel=100*(1-x$alpha), save = x$save, oe=x$oe)
  }

}
