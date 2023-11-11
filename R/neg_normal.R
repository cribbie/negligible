#' @title Negligible Effect Test for Normality of a Univariate Distribution
#'
#' @description This function allows researchers to test whether the distribution
#' of scores in a distribution has a Shapiro-Wilk W statistic that is negligibly different from 1.
#'
#'
#' @param x Distribution of Interest
#' @param eiL Lower Bound of the Negligible Effect Interval for W
#' @param nboot Number of Bootstrap Samples for computing the CIs
#' @param alpha Nominal Type I Error Rate
#' @param plot If the user prefers plots to be generated
#' @param data Dataset containing x
#' @param ... Extra arguments
#'
#' @return A \code{list} including the following:
#' \itemize{
#'   \item \code{sw} Sample Shapiro-Wilk W statistic
#'   \item \code{sskew} Sample skewness
#'   \item \code{skurt} Sample kurtosis
#'   \item \code{sddiff_mn_mdn} Standardized difference between the sample mean and median
#'   \item \code{sddiff_mn_trmn} Standardized difference between the sample mean and trimmed mean
#'   \item \code{lb} Lower bound of 1-alpha CI for W
#'   \item \code{eiL} Maximum W for which the degree of nonnormality is considered extreme
#' }
#' @export
#'
#' @details #' This function allows researchers to test whether the distribution
#' of scores in a distribution has a Shapiro-Wilk W statistic that is negligibly different from 1.
#' I.e., we are testing the null hypothesis that W is less than or equal to some
#' prespecified lower bound for W (i.e., the least extreme value of W that is
#' non-negligibly different from 1). We recommend .95 and .975 as liberal and conservative bounds,
#' respectively
#'
#' @author Rob Cribbie \email{cribbie@@yorku.ca} and
#'   Linda Farmus \email{lifarm@@yorku.ca}
#' @export neg.normal
#'
#' @examples
#' #Normal Distribution
#' xx<-stats::rnorm(200)
#' neg.normal(xx)
#' #Positive Skewed Distribution
#' xx<-stats::rchisq(200, df=3)
#' neg.normal(xx)

neg.normal <- function(x, eiL = .95, nboot = 1000,
                       plot = TRUE, alpha = .05,
                       data = NULL) {
  if (!is.null(data)) {
    x <- deparse(substitute(x))
    if(x=="NULL") {x<-NULL}
  }
  if (is.null(data)) {
    x <- x[!is.na(x)]
  }
  smean <- mean(x)
  smed <- stats::median(x)
  strmean <- mean(x, tr=.2)
  sddiff_mn_mdn <- (smed - smean)/stats::sd(x)
  sddiff_mn_trmn <- (strmean - smean)/stats::sd(x)
  sskew <- e1071::skewness(x, type = 2)
  skurt <- e1071::kurtosis(x, type = 2)
  ssd<- stats::sd(x)
  smin<- min(x)
  smax<- max(x)
  sw <- stats::shapiro.test(x)$statistic
  swp <- stats::shapiro.test(x)$p.value
  ifelse(swp <= alpha, decist <- paste0("The (traditional NHST) null hypothesis that the distribution is normal in form can be rejected at alpha = ", alpha),
         decist <- paste0("The (traditional NHST) null hypothesis that the distribution is normal in form cannot be rejected at alpha = ", alpha))
  wilk<-numeric(nboot)
  for (j in 1:nboot) { # for each bootstrap sample
    xx<-sample(x,replace=TRUE)
    wilk[j]<-stats::shapiro.test(xx)$statistic
  }
  lb <- stats::quantile(wilk, alpha)
  ifelse(lb >= eiL,
         decis <- "The null hypothesis that the degree of nonnormality is extreme (i.e., W <= eiL) can be rejected",
         decis <- "The null hypothesis that the degree of nonnormality is extreme (i.e., W <= eiL) cannot be rejected")



  pdw <- (sw-1)/abs(eiL-1)
  propd <- numeric(nboot)
  for (i in 1:nboot) {
    xx2<-sample(x, replace=TRUE)
    w<- stats::shapiro.test(xx2)$statistic
    propd[i]<-(w-1)/abs(eiL-1)
  }
  CI95L<-stats::quantile(propd, .025)
  CI95U<-stats::quantile(propd, .975)


  pdCI<- c(CI95L, 0)
  ret <- data.frame(sw = sw,
                    swp = swp,
                    decist = decist,
                    sskew = sskew,
                    skurt = skurt,
                    smean = smean,
                    smed = smed,
                    ssd = ssd,
                    smin = smin,
                    smax = smax,
                    sddiff_mn_mdn = sddiff_mn_mdn,
                    sddiff_mn_trmn = sddiff_mn_trmn,
                    strmean = strmean,
                    lb = lb,
                    eiL = eiL,
                    decis = decis,
                    pdw = pdw,
                    CI95L = CI95L,
                    CI95U = CI95U,
                    alpha = alpha,
                    pl = plot,
                    oe = "Shapiro-Wilk W")
  class(ret) <- "neg.normal"
  return(ret)
}

#' @rdname neg.normal
#' @param x object of class \code{neg.normal}
#' @export
#'
print.neg.normal <- function(x, ...) {
  cat("*******************************", "\n\n")
  cat("Descriptive Measures Regarding the Distribution\n\n")
  cat("Sample Skewness (Type 2): ", x$sskew, "\n", sep="")
  cat("Sample Kurtosis (Type 2): ", x$skurt, "\n\n", sep="")
  cat("Sample Measures of Central Tendency:", "\n\n")
  cat("Sample Mean:", x$smean, "\n")
  cat("Sample Median:", x$smed, "\n")
  cat("Standardized Difference between the Sample Mean and Median ((Median - Mean)/SD):", x$sddiff_mn_mdn, "\n")
  cat("Sample Trimmed Mean:", x$strmean, "\n")
  cat("Standardized Difference between the Sample Mean and Trimmed Mean ((Trimmed Mean - Mean)/SD):", x$sddiff_mn_trmn, "\n\n")
  cat("Sample Measures of Variability:", "\n\n")
  cat("Sample Standard Deviation:", x$ssd, "\n")
  cat("Sample Min:", x$smin, "\n")
  cat("Sample Max:", x$smax, "\n")
  cat("\n\n")
  cat("*******************************", "\n\n")
  cat("Traditional Shapiro Wilk Test Statistic:", "\n", x$sw, "\n\n")
  cat("Traditional Shapiro Wilk Test p-Value:", "\n", x$swp, "\n\n")
  cat(x$decist, "\n\n")
  cat("*******************************", "\n\n")
  cat("Statistical Test of a Negligible Difference Between the Target/Population Distribution and a Theoretical Normal Distribution\n\n")
  cat("Shapiro-Wilk Statistic:", x$sw, "\n")
  cat("Lower Bound of the 100(1-2*", x$alpha, ")% CI: ", x$lb, "\n", sep="")
  cat("Lower Bound of the Negligible Effect (Equivalence) Interval:", x$eiL, "\n")
  cat("Negligible Effect (Equivalence) Testing Decision:", x$decis, "\n\n")
  cat("*******************************", "\n\n")
}
