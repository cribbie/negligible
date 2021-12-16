#' @title Negligible Effect Test on the Difference between the Means of Independent Populations
#'
#' @description This function allows researchers to test whether the difference
#' between the means of two independent populations is negligible, where
#' negligible represents the smallest meaningful effect size (MMES, which
#' in this case the effect is the mean difference)
#'
#' @param x Data for Group 1 (if dv and iv are omitted)
#' @param y Data for Group 2 (if dv and iv are omitted)
#' @param dv Dependent Variable (if x and y are omitted)
#' @param iv Dichotomous Predictor/Independent Variable (if x and y are omitted)
#' @param eil Lower Bound of the Equivalence Interval
#' @param eiu Upper Bound of the Equivalence Interval
#' @param varequiv Are the population variances assumed to be equal? Population variances are assumed to be unequal if normality=FALSE.
#' @param normality Are the population variances (and hence the residuals) assumed to be normally distributed?
#' @param tr Proportion of trimming from each tail (relevant if normality = FALSE)
#' @param nboot Number of bootstrap samples for calculating CIs
#' @param alpha Nominal Type I Error rate
#' @param plot Should a plot of the results be produced?
#' @param saveplot Should the plot be saved?
#' @param data Dataset containing x/y or iv/dv
#'
#' @return returns a \code{list} containing each analysis and their respective statistics
#'   and decision
#'
#' @author Rob Cribbie \email{cribbie@@yorku.ca}
#'   R. Philip Chalmers \email{chalmrp@@yorku.ca}
#'   Naomi Martinez Gutierrez \email{naomimg@@yorku.ca}
#' @export neg.twoindmeans
#'
#' @examples
#' \dontrun{
#' indvar<-rep(c("a","b"),c(10,12))
#' depvar<-rnorm(22)
#' d<-data.frame(indvar,depvar)
#' neg.twoindmeans(dv=depvar,iv=indvar,eil=-1,eiu=1,plot=TRUE,data=d)
#' neg.twoindmeans(dv=depvar,iv=indvar,eil=-1,eiu=1)
#' neg.twoindmeans(x=depvar[iv=="a"],y=depvar[iv=="b"],eil=-1,eiu=1)
#' xx<-neg.twoindmeans(dv=depvar,iv=indvar,eil=-1,eiu=1)
#' xx$decis
#' }
#'

neg.twoindmeans <- function(x = NULL, y = NULL, dv = NULL, iv = NULL,
                             eil, eiu, varequiv = FALSE, normality = FALSE,
                             tr = 0.2, nboot = 500, alpha = 0.05,
                             plot = TRUE, saveplot = FALSE, data=NULL) {
  if (is.null(x) & is.null(y) & !is.null(data)) {
    dv<-deparse(substitute(dv))
    iv<-deparse(substitute(iv))
    dv<-as.numeric(data[[dv]])
    iv<-as.factor(data[[iv]])
    dat<- data.frame(dv,iv)
    dat<- dat[stats::complete.cases(dat),]
    x<-dat$dv[dat$iv==levels(dat$iv)[1]]
    y<-dat$dv[dat$iv==levels(dat$iv)[2]]
  }
  if (is.null(dv) & is.null(iv) & !is.null(data)) {
    x<-deparse(substitute(x))
    y<-deparse(substitute(y))
    x<-as.numeric(data[[x]])
    y<-as.numeric(data[[y]])
  }
  if (is.null(dv) & is.null(iv) & is.null(data)) {
      x <- x[!is.na(x)]
      y <- y[!is.na(y)]
  }
  if (is.null(x) & is.null(y) & is.null(data)) {
    d<-data.frame(iv,dv)
    d$iv<-factor(d$iv)
    d <- d[stats::complete.cases(d),]
    x<-d$dv[d$iv==levels(d$iv)[1]]
    y<-d$dv[d$iv==levels(d$iv)[2]]
  }
  if (normality==TRUE) {
    if (varequiv == FALSE) {
      denom <- sqrt((stats::var(x)/length(x)) + (stats::var(y)/length(y)))
      t1 <- (mean(x) - mean(y) - eiu)/denom
      t2 <- (mean(x) - mean(y) - eil)/denom
      dft <- (((stats::var(x)/length(x)) + (stats::var(y)/length(y)))^2)/
        ((stats::var(x)^2/(length(x)^2 * (length(x) - 1))) +
           (stats::var(y)^2/(length(y)^2 *(length(y) - 1))))
      probt1 <- stats::pt(t1, dft, lower.tail = T)
      probt2 <- stats::pt(t2, dft, lower.tail = F)
      ifelse(probt1 <= alpha & probt2 <= alpha,
             decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval can be rejected. Be sure to interpret the magnitude of the effect size.",
             decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval cannot be rejected. Be sure to interpret the magnitude of the effect size.")
      effsize_raw<-mean(x)-mean(y)
      effsize_d<-(mean(x)-mean(y))/sqrt((stats::var(x)+stats::var(y))/2)
      ifelse(sign(mean(x)-mean(y))==sign(eil),
             ein<-eil,ein<-eiu)
      effsize_pd<-(mean(x)-mean(y))/abs(ein) #change here
      bsraw<-numeric(nboot)
      bsd<-numeric(nboot)
      bspd<-numeric(nboot)
      for (i in 1:nboot) {
        bsx<-sample(x,length(x),replace=TRUE)
        bsy<-sample(y,length(y),replace=TRUE)
        bsraw[i]<-mean(bsx)-mean(bsy)
        bsd[i]<-(mean(bsx)-mean(bsy))/sqrt((stats::var(bsx)+stats::var(bsy))/2)
        ifelse(sign(mean(bsx)-mean(bsy))==sign(eil),
               ein2<-eil,ein2<-eiu)
        bspd[i]<-(mean(bsx)-mean(bsy))/abs(ein2) #change here
      }
      ci_raw2<-stats::quantile(bsraw,probs=c(alpha/2,1-alpha/2))
      ci_raw<-stats::quantile(bsraw,probs=c(alpha,1-alpha))
      ci_d<-stats::quantile(bsd,probs=c(alpha,1-alpha))
      ci_pd<-stats::quantile(bspd,probs=c(alpha/2,1-alpha/2))

      title <- "Schuirmann-Welch Test of the Equivalence of Two Independent Groups"
    }
    if (varequiv == TRUE) {
      denom <- sqrt(((((length(x) - 1) * stats::sd(x)^2) +
                ((length(y) - 1) * stats::sd(y)^2))/(length(x) +
                length(y) - 2)) * (1/length(x) + 1/length(y)))
      t1 <- (mean(x) - mean(y) - eiu)/denom
      t2 <- (mean(x) - mean(y) - eil)/denom
      dft <- length(x) + length(y) - 2
      probt1 <- stats::pt(t1, dft, lower.tail = T)
      probt2 <- stats::pt(t2, dft, lower.tail = F)
      ifelse(probt1 <= alpha & probt2 <= alpha,
             decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval can be rejected. Be sure to interpret the magnitude of the effect size.",
             decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval cannot be rejected. Be sure to interpret the magnitude of the effect size.")
      effsize_raw<-mean(x)-mean(y)
      effsize_d<-(mean(x)-mean(y))/
        sqrt((((length(x) - 1) * stats::sd(x)^2) + ((length(y) - 1) * stats::sd(y)^2))/
                (length(x) + length(y) - 2))
      ifelse(sign(mean(x)-mean(y))==sign(eil),
             ein<-eil,ein<-eiu)
      effsize_pd<-(mean(x)-mean(y))/abs(ein)
      bsraw<-numeric(nboot)
      bsd<-numeric(nboot)
      bspd<-numeric(nboot)
      for (i in 1:nboot) {
        bsx<-sample(x, length(x),replace=TRUE)
        bsy<-sample(y, length(y),replace=TRUE)
        bsraw[i]<-mean(bsx)-mean(bsy)
        bsd[i]<-(mean(bsx)-mean(bsy))/
          sqrt((((length(bsx) - 1) * stats::sd(bsx)^2) + ((length(bsy) - 1) * stats::sd(bsy)^2))/
                 (length(bsx) + length(bsy) - 2))
        ifelse(sign(mean(bsx)-mean(bsy))==sign(eil),
               ein2<-eil,ein2<-eiu)
        bspd[i]<-(mean(bsx)-mean(bsy))/abs(ein2)
      }
      ci_raw2<-stats::quantile(bsraw,probs=c(alpha/2,1-alpha/2))
      ci_raw<-stats::quantile(bsraw,probs=c(alpha,1-alpha))
      ci_d<-stats::quantile(bsd,probs=c(alpha,1-alpha))
      ci_pd<-stats::quantile(bspd,probs=c(alpha/2,1-alpha/2))
      title <- "Schuirmann's Test of the Equivalence of Two Independent Groups"
    }

  }
  if (normality == FALSE) {
    h1 <- length(x) - 2 * floor(tr * length(x))
    h2 <- length(y) - 2 * floor(tr * length(y))
    q1 <- (length(x) - 1) * WRS2::winvar(x, tr)/(h1 *
             (h1 - 1))
    q2 <- (length(y) - 1) * WRS2::winvar(y, tr)/(h2 *
             (h2 - 1))
    dft <- (q1 + q2)^2/((q1^2/(h1 - 1)) + (q2^2/(h2 -
                                                   1)))
    crit <- stats::qt(1 - alpha/2, dft)
    dif1 <- mean(x, tr) - mean(y, tr) - eiu
    dif2 <- mean(x, tr) - mean(y, tr) - eil
    t1 <- dif1/sqrt(q1 + q2)
    t2 <- dif2/sqrt(q1 + q2)
    probt1 <- stats::pt(t1, dft)
    probt2 <- 1 - stats::pt(t2, dft)
    ifelse(probt1 <= alpha & probt2 <= alpha,
           decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval can be rejected. Be sure to interpret the magnitude of the effect size.",
           decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval cannot be rejected. Be sure to interpret the magnitude of the effect size.")
    effsize_raw<-mean(x,tr=tr)-mean(y,tr=tr)
    effsize_d <- .642*(mean(x,tr=tr) - mean(y,tr=tr))/sqrt((WRS2::winvar(x,tr=.2)+WRS2::winvar(y,tr=.2))/2)
    ifelse(sign(mean(x,tr=tr)-mean(y,tr=tr))==sign(eil),
           ein<-eil,ein<-eiu)
    effsize_pd<-(mean(x,tr=tr)-mean(y,tr=tr))/abs(ein)
    bsraw<-numeric(nboot)
    bsd<-numeric(nboot)
    bspd<-numeric(nboot)
    for (i in 1:nboot) {
      bsx<-sample(x, length(x),replace=TRUE)
      bsy<-sample(y,length(y),replace=TRUE)
      bsraw[i]<-mean(bsx,tr=tr)-mean(bsy,tr=tr)
      bsd[i]<-.642*(mean(bsx,tr=tr)-mean(bsy,tr=tr))/
        sqrt((WRS2::winvar(bsx,tr=.2)+WRS2::winvar(bsy,tr=.2))/2)
      ifelse(sign(mean(bsx,tr=tr)-mean(bsy,tr=tr))==sign(eil),
             ein2<-eil,ein2<-eiu)
      bspd[i]<-(mean(bsx,tr=tr)-mean(bsy,tr=tr))/abs(ein2)
    }
    ci_raw2<-stats::quantile(bsraw,probs=c(alpha/2,1-alpha/2))
    ci_raw<-stats::quantile(bsraw,probs=c(alpha,1-alpha))
    ci_d<-stats::quantile(bsd,probs=c(alpha,1-alpha))
    ci_pd<-stats::quantile(bspd,probs=c(alpha/2,1-alpha/2))
    title <- "Schuirmann-Yuen Test of the Equivalence of Two Independent Groups"
  }
  ret <- data.frame(meanx = mean(x),
                    meany = mean(y),
                    trmeanx = mean(x, tr),
                    trmeany = mean(y, tr),
                    sdx = stats::sd(x),
                    sdy = stats::sd(y),
                    madx = stats::mad(x),
                    mady = stats::mad(y),
                    eil = eil,
                    eiu = eiu,
                    eisign = ein,
                    effsizeraw = effsize_raw,
                    cilraw2 = ci_raw2[1],
                    ciuraw2 = ci_raw2[2],
                    cilraw = ci_raw[1],
                    ciuraw = ci_raw[2],
                    effsized = effsize_d,
                    cild = ci_d[1],
                    ciud = ci_d[2],
                    effsizepd = effsize_pd,
                    cilpd = ci_pd[1],
                    ciupd = ci_pd[2],
                    t1 = t1,
                    t2 = t2,
                    df1 = dft,
                    df2 = dft,
                    pval1 = round(probt1,5),
                    pval2 = round(probt2,5),
                    decis = decis,
                    alpha = alpha,
                    title1 = title,
                    pl = plot,
                    save = saveplot)
  class(ret) <- "neg.twoindmeans"
  return(ret)
}

#' @rdname neg.twoindmeans
#' @param x object of class \code{neg.twoindmeans}
#' @param ... extra arguments
#' @return
#' @export
#'

print.neg.twoindmeans <- function(x, ...) {
  cat("---- Equivalence of Two Independent Groups----\n\n")
  cat("Test Statistic:", "\n", x$title1, "\n\n")
  cat("**********************\n\n")
  cat("Group Means: ", x$meanx, ", ", x$meany, "\n", sep="")
  cat("Group Trimmed Means: ", x$trmeanx, ", ", x$trmeany, "\n", sep="")
  cat("Group SDs: ", x$sdx, ", ", x$sdy, "\n", sep="")
  cat("Group MADs: ", x$madx, ", ", x$mady, "\n\n", sep="")
  cat("**********************\n\n")
  cat("Standardized Mean Difference (SMD):", x$effsized, "\n")
  cat(100*(1-2*x$alpha), "% CI for SMD:", "(", x$cild, ", ", x$ciud, ")", "\n\n", sep="")
  cat("**********************\n\n")
  cat("Equivalence Interval:","Lower =", x$eil, ",", "Upper =", x$eiu, "\n\n")
  cat("**********************\n\n")
  cat("Raw Mean Difference (MD):", x$effsizeraw, "\n")
  cat(100*(1-2*x$alpha), "% CI for MD:", "(", x$cilraw, ", ", x$ciuraw, ")", "\n", sep="")
  cat(100*(1-x$alpha), "% CI for MD:", "(", x$cilraw2, ", ", x$ciuraw2, ")", "\n\n", sep="")
  cat("**********************\n\n")
  cat("Proportional Distance (PD):", x$effsizepd, "\n")
  cat(100*(1-x$alpha), "% CI for PD:", "(", x$cilpd, ", ", x$ciupd, ")", "\n\n", sep="")
  cat("**********************\n\n")
  cat("TOST Test Statistics:", "\n\n", "Ho: mu1-mu2>=eiu:", "\n", "t =", x$t1,
      " (df = ", x$df1, ")", ", p = ", x$pval1, "\n\n", "Ho: mu1-mu2<=eil:", "\n", "t =", x$t2,
      " (df =", x$df1, ")", ", p = ", x$pval2, "\n\n", sep="")
  cat("**********************\n\n")
  cat("NHST Decision:", "\n", x$decis, "\n\n", sep="")
  cat("**********************\n\n")

  if (x$pl == TRUE) {
    neg.pd(effect = x$effsizeraw,
                         PD = x$effsizepd,
                         EIsign=x$eisign,
                         PDcil=x$cilpd,
                         PDciu=x$ciupd,
                         cil=x$cilraw,
                         ciu=x$ciuraw,
                         Elevel=100*(1-2*x$alpha),
                         Plevel=100*(1-x$alpha),
                         save = x$save)
  }
}

