#' @title Negligible Effect Test for Variances of Independent Populations
#'
#' @description This function allows researchers to test whether the difference
#' in the variances of independent populations is negligible, where
#' negligible represents the smallest meaningful effect size (MMES, where
#' in this case the effect is the difference in population variances)
#'
#'
#' @param dv Outcome Variable
#' @param iv Independent Variable
#' @param eps Used to Establish the Equivalence Bound (Conservative: .25; Liberal: .50, according to Wellek, 2010)
#' @param alpha Nominal Type I Error Rate
#' @param na.rm Missing Data Treatment
#' @param data Dataset containing dv and iv
#' @param ... Extra arguments
#'
#' @return A \code{list} including the following:
#' \itemize{
#'   \item \code{vars} Sample variances
#'   \item \code{sds} Sample standard deviations
#'   \item \code{mads} Sample median absolute deviations
#'   \item \code{ratio} Ratio of the largest to smallest variance
#'   \item \code{eps} Epsilon (e) can be described as the maximum difference in the variances that one would consider to be unimportant (see Details).
#'   \item \code{LWW_md} Levene-Wellek-Welch statistic based on the median.
#'   \item \code{crit_LWW_md} Critical value for the Levene-Wellek-Welch statistic based on the median.
#'   \item \code{alpha} Nominal Type I error rate
#' }
#' @export
#' @details This function evaluates whether the difference in the population variances of J independent groups can be considered negligible (i.e., the population variances can be considered equivalent).
#'
#' The user provides the name of the outcome/dependent variable (should be continuous) and the name of Independent Variable (predictor, should be a factor), as well as the epsilon value (eps) which determines the smallest difference in variances that can be considered non-negligible.
#'
#' Wellek (2010) suggests liberal and conservative values of eps = .50 and eps = .25, respectively. See Wellek, 2010, pp. 16, 17, 22, for details.
#'
#' See Mara & Cribbie (2018): https://doi.org/10.1080/00220973.2017.1301356
#'
#' @author Rob Cribbie \email{cribbie@@yorku.ca} and
#'   Constance Mara \email{Constance.Mara@@cchmc.org}
#' @export neg.indvars
#'
#' @examples
#' #Two Group Example
#' indvar<-rep(c("a","b"),c(10,12))
#' depvar<-rnorm(22)
#' d<-data.frame(indvar,depvar)
#' neg.indvars(depvar,indvar)
#' neg.indvars(dv=depvar,iv=indvar,eps=.25,data=d)
#' neg.indvars(dv=depvar,iv=indvar,eps=.5)
#'
#' #Four Group Example
#' indvar<-rep(c("a","b","c","d"),c(10,12,15,13))
#' depvar<-rnorm(50)
#' d<-data.frame(indvar,depvar)
#' neg.indvars(dv=depvar,iv=indvar,eps=.25,data=d)
#' neg.indvars(dv=depvar,iv=indvar)

neg.indvars<-function(dv, iv, eps = .5, alpha = .05,
            na.rm=TRUE, data = NULL, ...) {
  if (!is.null(data)) {
    dv<-deparse(substitute(dv))
    iv<-deparse(substitute(iv))
    dv<-as.numeric(data[[dv]])
    iv<-as.factor(data[[iv]])
  }
  medians <- tapply(dv, iv, stats::median)
  n <- tapply(dv, iv, length)
  resp.median <- abs(dv - medians[iv])
  ngroup<-length(unique(iv))
  vars<-(tapply(dv, iv, stats::var))
  ratio_lsvar<-max(vars)/min(vars)

  ## Equivalence test for Equivalence of Variances ##
  LWW_md<-stats::oneway.test(resp.median~iv)$statistic*((ngroup-1)/(mean(n)))
  crit_LWW_md<-((ngroup-1)/(mean(n)))*stats::qf(p=alpha, df1=ngroup-1,
                              df2=stats::oneway.test(resp.median~iv)$parameter[2],
                              ncp=(mean(n))*eps^2)
  ifelse (LWW_md <= crit_LWW_md,
          decis_equiv<-"The null hypothesis that the differences between the population variances falls outside the equivalence interval can be rejected. A negligible difference among the population variances can be concluded. Be sure to interpret the magnitude (and precision) of the effect size.",
          decis_equiv<-"The null hypothesis that the differences between the population variances falls outside the equivalence interval cannot be rejected. A negligible difference among the population variances cannot be concluded. Be sure to interpret the magnitude (and precision) of the effect size.")
  ret <- data.frame(vars = vars,
                    sds = tapply(dv,iv,stats::sd),
                    mads = tapply(dv,iv,stats::mad),
                    ratio = ratio_lsvar,
                    eps = eps,
                    LWW_md = LWW_md,
                    crit_LWW_md = crit_LWW_md,
                    decis = decis_equiv,
                    alpha = alpha)
  class(ret) <- "neg.indvars"
  return(ret)
}

#' @rdname neg.indvars
#' @param x object of class \code{neg.indvars}
#' @return
#' @export
#'
print.neg.indvars <- function(x, ...) {
  cat("-- Equivalence of Population Variances --\n")
  cat("-- Independent Groups --\n\n")
  cat("Group Variances: ", "\n", x$vars, "\n\n")
  cat("Group Standard Deviations: ", "\n", x$sds, "\n\n")
  cat("Group Median Absolute Deviations: ", "\n", x$mads,"\n\n")
  cat("**********************\n\n")
  cat("Ratio of Largest to Smallest Variances: ", "\n", x$ratio[1], "\n\n")
  cat("**********************\n\n")
  cat("Epsilon Value (establish the Equivalence Interval): ", "\n", x$eps[1], "\n\n")
  cat("Levene-Wellek-Welch (LWW) Statistic: ", "\n", x$LWW_md[1], "\n\n")
  cat("Critical Value for LWW: ", "\n", x$crit_LWW_md[1], "\n\n")
  cat("NHST Decision: ", "\n", x$decis[1], "\n\n")
  cat("**********************\n\n")
}

