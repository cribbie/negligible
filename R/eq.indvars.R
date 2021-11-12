#' Title Function for Evaluating the Equivalence of the Variances of Independent Populations
#'
#' @param dv Outcome Variable
#' @param iv Independent Variable
#' @param eps Used to Establish the Equivalence Interval
#' @param alpha Nominal Type I Error Rate
#' @param na.rm Missing Data Treatment
#' @param data Dataset containing dv and iv
#' @param ... Extra arguments
#'
#' @return returns a \code{list} containing each analysis and their respective statistics
#'   and decision
#' @author Rob Cribbie \email{cribbie@@yorku.ca} and
#'   Constance Mara \email{Constance.Mara@@cchmc.org}
#' @export eq.indvars
#'
#' @examples
#' \dontrun{
#' #Two Group Example
#' indvar<-rep(c("a","b"),c(10,12))
#' depvar<-rnorm(22)
#' d<-data.frame(indvar,depvar)
#' var_equiv(depvar,indvar)
#' var_equiv(dv=depvar,iv=indvar,eps=.25,data=d)
#' var_equiv(dv=depvar,iv=indvar,eps=.5)
#'
#' #Four Group Example
#' indvar<-rep(c("a","b","c","d"),c(10,12,15,13))
#' depvar<-rnorm(50)
#' d<-data.frame(indvar,depvar)
#' var_equiv(dv=depvar,iv=indvar,eps=.25,data=d)
#' var_equiv(dv=depvar,iv=indvar)
#' }

eq.indvars<-function(dv, iv, eps = .5, alpha = .05,
            na.rm=TRUE, data = NULL, ...) {
  if (!is.null(data)) {
    dv<-deparse(substitute(dv))
    iv<-deparse(substitute(iv))
    dv<-as.numeric(data[[dv]])
    iv<-as.factor(data[[iv]])
  }
  medians <- tapply(dv, iv, median)
  n <- tapply(dv, iv, length)
  resp.median <- abs(dv - medians[iv])
  ngroup<-length(unique(iv))
  vars<-(tapply(dv, iv, var))
  ratio_lsvar<-max(vars)/min(vars)

  ## Equivalence test for Equivalence of Variances ##
  LWW_md<-oneway.test(resp.median~iv)$statistic*((ngroup-1)/(mean(n)))
  crit_LWW_md<-((ngroup-1)/(mean(n)))*qf(p=alpha, df1=ngroup-1,
                              df2=oneway.test(resp.median~iv)$parameter[2],
                              ncp=(mean(n))*eps^2)
  ifelse (LWW_md <= crit_LWW_md,
          decis_equiv<-"The null hypothesis that the differences between the population variances falls outside the equivalence interval can be rejected.",
          decis_equiv<-"The null hypothesis that the differences between the population variances falls outside the equivalence interval cannot be rejected")
  ret <- data.frame(vars = vars,
                    sds = tapply(dv,iv,sd),
                    mads = tapply(dv,iv,mad),
                    ratio = ratio_lsvar,
                    eps = eps,
                    LWW_md = LWW_md,
                    crit_LWW_md = crit_LWW_md,
                    decis = decis_equiv,
                    alpha = alpha)
  class(ret) <- "eq.indvars"
  return(ret)
}

#' @rdname eq.indvars
#' @param x object of class \code{eq.indvars}
#' @return
#' @export
#'
print.eq.indvars <- function(x, ...) {
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

