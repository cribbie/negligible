#' @title Negligible Effect Test on the Difference between the Means of Independent Populations
#'
#' @description This function allows researchers to test whether the difference
#' between the means of two independent populations is negligible, where
#' negligible represents the smallest meaningful effect size (MMES, which
#' in this case the effect is the mean difference)
#'
#' @param v1 Data for Group 1 (if dv and iv are omitted)
#' @param v2 Data for Group 2 (if dv and iv are omitted)
#' @param dv Dependent Variable (if v1 and v2 are omitted)
#' @param iv Dichotomous Predictor/Independent Variable (if v1 and v2 are omitted)
#' @param eiL Lower Bound of the Equivalence Interval
#' @param eiU Upper Bound of the Equivalence Interval
#' @param varequiv Are the population variances assumed to be equal? Population variances are assumed to be unequal if normality=FALSE.
#' @param normality Are the population variances (and hence the residuals) assumed to be normally distributed?
#' @param tr Proportion of trimming from each tail (relevant if normality = FALSE)
#' @param nboot Number of bootstrap samples for calculating CIs
#' @param alpha Nominal Type I Error rate
#' @param plot Should a plot of the results be produced?
#' @param saveplot Should the plot be saved?
#' @param data Dataset containing v1/v2 or iv/dv
#'
#' @return A \code{list} including the following:
#' \itemize{
#'   \item \code{meanx} Sample mean of the first population/group.
#'   \item \code{meany} Sample mean of the second population/group.
#'   \item \code{trmeanx} Sample trimmed mean of the first population/group.
#'   \item \code{trmeany} Sample trimmed mean of the second population/group.
#'   \item \code{sdx} Sample standard deviation of the first population/group.
#'   \item \code{sdy} Sample standard deviation of the second population/group.
#'   \item \code{madx} Sample median absolute deviation of the first population/group.
#'   \item \code{mady} Sample median absolute deviation of the second population/group.
#'   \item \code{eiL} Lower bound of the negligible effect (equivalence) interval.
#'   \item \code{eiU} Upper bound of the negligible effect (equivalence) interval.
#'   \item \code{effsizeraw} Simple difference in the means (or trimmed means if normality = FALSE)
#'   \item \code{cilraw2} Lower bound of the 1-alpha CI for the raw mean difference.
#'   \item \code{ciuraw2} Upper bound of the 1-alpha CI for the raw mean difference.
#'   \item \code{cilraw} Lower bound of the 1-2*alpha CI for the raw mean difference.
#'   \item \code{ciuraw} Upper bound of the 1-2*alpha CI for the raw mean difference.
#'   \item \code{effsized} Standardized mean (or trimmed mean if normality = FALSE) difference.
#'   \item \code{cild} Lower bound of the 1-alpha CI for the standardized mean (or trimmed mean if normality = FALSE) difference.
#'   \item \code{ciud} Upper bound of the 1-alpha CI for the standardized mean (or trimmed mean if normality = FALSE) difference.
#'   \item \code{effsizepd} Proportional distance statistic.
#'   \item \code{cilpd} Lower bound of the 1-alpha CI for the proportional distance statistic.
#'   \item \code{ciupd} Upper bound of the 1-alpha CI for the proportional distance statistic.
#'   \item \code{t1} First t-statistic from the TOST procedure.
#'   \item \code{t1} Second t-statistic from the TOST procedure.
#'   \item \code{df1} Degrees of freedom for the first t-statistic from the TOST procedure.
#'   \item \code{df2} Degrees of freedom for the second t-statistic from the TOST procedure.
#'   \item \code{p1} p value associated with the first t-statistic from the TOST procedure.
#'   \item \code{p2} p value associated with the second t-statistic from the TOST procedure.
#'   \item \code{alpha} Nominal Type I error rate
#' }
#' @export
#' @details This function evaluates whether the difference in the means of 2 independent populations can be considered negligible (i.e., the population means can be considered equivalent).
#'
#' The user specifies either the data associated with the first and second groups/populations (iv1, iv2, both should be continuous) or specifies the Indepedent Variable/Predictor (iv, should be a factor) and the Dependent Variable (outcome, should be continuous). A 'data' statement can be used if the variables are stored in an R dataset.
#'
#' The user must also specify the lower and upper bounds of the negligible effect (equivalence) interval. These are specified in the original units of the outcome variable.
#'
#' The arguments 'varequiv' and 'normality' control what test statistic is adopted. If varequiv = TRUE and normality = TRUE the ordinary Student t statistic is adopted. If varequiv = FALSE and normality = TRUE the Welch t statistic is adopted. If normality = FALSE  the ordinary Student t statistic is adopted.  d
#' @author Rob Cribbie \email{cribbie@@yorku.ca}
#'   R. Philip Chalmers \email{chalmrp@@yorku.ca}
#'   Naomi Martinez Gutierrez \email{naomimg@@yorku.ca}
#' @export neg.twoindmeans
#'
#' @examples
#' indvar<-rep(c("a","b"),c(10,12))
#' depvar<-rnorm(22)
#' d<-data.frame(indvar,depvar)
#' neg.twoindmeans(dv=depvar,iv=indvar,eiL=-1,eiU=1,plot=TRUE,data=d)
#' neg.twoindmeans(dv=depvar,iv=indvar,eiL=-1,eiU=1)
#' neg.twoindmeans(v1=depvar[indvar=="a"],v2=depvar[indvar=="b"],eiL=-1,eiU=1)
#' xx<-neg.twoindmeans(dv=depvar,iv=indvar,eiL=-1,eiU=1)
#' xx$decis
neg.twoindmeans <- function(v1 = NULL, v2 = NULL, dv = NULL, iv = NULL,
                            eiL, eiU, varequiv = FALSE, normality = FALSE,
                            tr = 0.2, nboot = 500, alpha = 0.05,
                            plot = TRUE, saveplot = FALSE, data=NULL) {
  if (!is.null(data)) {
    v1 <- deparse(substitute(v1))
    v2 <- deparse(substitute(v2))
    iv <- deparse(substitute(iv))
    dv <- deparse(substitute(dv))
    if(v1=="NULL") {v1<-NULL}
    if(v2=="NULL") {v2<-NULL}
    if(iv=="NULL") {iv<-NULL}
    if(dv=="NULL") {dv<-NULL}
  }
  if (is.null(v1) & is.null(v2) & !is.null(data)) {
    dv<-as.numeric(data[[dv]])
    iv<-as.factor(data[[iv]])
    dat<- data.frame(dv,iv)
    dat<- dat[stats::complete.cases(dat),]
    v1<-dat$dv[dat$iv==levels(dat$iv)[1]]
    v2<-dat$dv[dat$iv==levels(dat$iv)[2]]
    lev1<-levels(dat$iv)[1]
    lev2<-levels(dat$iv)[2]
  }
  if (is.null(dv) & is.null(iv) & !is.null(data)) {
    v1<-as.numeric(data[[v1]])
    v2<-as.numeric(data[[v2]])
    lev1<-"v1"
    lev2<-"v2"
  }
  if (is.null(dv) & is.null(iv) & is.null(data)) {
    v1 <- v1[!is.na(v1)]
    v2 <- v2[!is.na(v2)]
    lev1<-"v1"
    lev2<-"v2"
  }
  if (is.null(v1) & is.null(v2) & is.null(data)) {
    d<-data.frame(iv,dv)
    d$iv<-factor(d$iv)
    d <- d[stats::complete.cases(d),]
    v1<-d$dv[d$iv==levels(d$iv)[1]]
    v2<-d$dv[d$iv==levels(d$iv)[2]]
    lev1<-levels(d$iv)[1]
    lev2<-levels(d$iv)[2]
  }
  if (normality==TRUE) {
    if (varequiv == FALSE) {
      denom <- sqrt((stats::var(v1)/length(v1)) + (stats::var(v2)/length(v2)))
      t1 <- (mean(v1) - mean(v2) - eiU)/denom
      t2 <- (mean(v1) - mean(v2) - eiL)/denom
      dft <- (((stats::var(v1)/length(v1)) + (stats::var(v2)/length(v2)))^2)/
        ((stats::var(v1)^2/(length(v1)^2 * (length(v1) - 1))) +
           (stats::var(v2)^2/(length(v2)^2 *(length(v2) - 1))))
      probt1 <- stats::pt(t1, dft, lower.tail = T)
      probt2 <- stats::pt(t2, dft, lower.tail = F)
      ifelse(probt1 <= alpha & probt2 <= alpha,
             decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval can be rejected. A negligible difference in the population means can be concluded. Be sure to interpret the magnitude (and precision) of the effect size.",
             decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval cannot be rejected. There is insufficient evidence to conclude a negligible difference in the population means. Be sure to interpret the magnitude (and precision) of the effect size.")
      effsize_raw<-mean(v1)-mean(v2)
      effsize_d<-(mean(v1)-mean(v2))/sqrt((stats::var(v1)+stats::var(v2))/2)
      ifelse(sign(mean(v1)-mean(v2))==sign(eiL),
             ein<-eiL,ein<-eiU)
      effsize_pd<-(mean(v1)-mean(v2))/abs(ein) #change here
      bsraw<-numeric(nboot)
      bsd<-numeric(nboot)
      bspd<-numeric(nboot)
      for (i in 1:nboot) {
        bsx<-sample(v1,length(v1),replace=TRUE)
        bsy<-sample(v2,length(v2),replace=TRUE)
        bsraw[i]<-mean(bsx)-mean(bsy)
        bsd[i]<-(mean(bsx)-mean(bsy))/sqrt((stats::var(bsx)+stats::var(bsy))/2)
        ifelse(sign(mean(bsx)-mean(bsy))==sign(eiL),
               ein2<-eiL,ein2<-eiU)
        bspd[i]<-(mean(bsx)-mean(bsy))/abs(ein2) #change here
      }
      ci_raw2<-stats::quantile(bsraw,probs=c(alpha/2,1-alpha/2))
      ci_raw<-stats::quantile(bsraw,probs=c(alpha,1-alpha))
      ci_d<-stats::quantile(bsd,probs=c(alpha,1-alpha))
      ci_pd<-stats::quantile(bspd,probs=c(alpha/2,1-alpha/2))
      
      title <- "Schuirmann-Welch Test of the Equivalence of Two Independent Groups"
    }
    if (varequiv == TRUE) {
      denom <- sqrt(((((length(v1) - 1) * stats::sd(v1)^2) +
                        ((length(v2) - 1) * stats::sd(v2)^2))/(length(v1) +
                                                                 length(v2) - 2)) * (1/length(v1) + 1/length(v2)))
      t1 <- (mean(v1) - mean(v2) - eiU)/denom
      t2 <- (mean(v1) - mean(v2) - eiL)/denom
      dft <- length(v1) + length(v2) - 2
      probt1 <- stats::pt(t1, dft, lower.tail = T)
      probt2 <- stats::pt(t2, dft, lower.tail = F)
      ifelse(probt1 <= alpha & probt2 <= alpha,
             decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval can be rejected. A negligible difference in the population means can be concluded. Be sure to interpret the magnitude of the effect size.",
             decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval cannot be rejected. A negligible difference in the population means cannot be concluded. Be sure to interpret the magnitude of the effect size.")
      effsize_raw<-mean(v1)-mean(v2)
      effsize_d<-(mean(v1)-mean(v2))/
        sqrt((((length(v1) - 1) * stats::sd(v1)^2) + ((length(v2) - 1) * stats::sd(v2)^2))/
               (length(v1) + length(v2) - 2))
      ifelse(sign(mean(v1)-mean(v2))==sign(eiL),
             ein<-eiL,ein<-eiU)
      effsize_pd<-(mean(v1)-mean(v2))/abs(ein)
      bsraw<-numeric(nboot)
      bsd<-numeric(nboot)
      bspd<-numeric(nboot)
      for (i in 1:nboot) {
        bsx<-sample(v1, length(v1),replace=TRUE)
        bsy<-sample(v2, length(v2),replace=TRUE)
        bsraw[i]<-mean(bsx)-mean(bsy)
        bsd[i]<-(mean(bsx)-mean(bsy))/
          sqrt((((length(bsx) - 1) * stats::sd(bsx)^2) + ((length(bsy) - 1) * stats::sd(bsy)^2))/
                 (length(bsx) + length(bsy) - 2))
        ifelse(sign(mean(bsx)-mean(bsy))==sign(eiL),
               ein2<-eiL,ein2<-eiU)
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
    h1 <- length(v1) - 2 * floor(tr * length(v1))
    h2 <- length(v2) - 2 * floor(tr * length(v2))
    q1 <- (length(v1) - 1) * WRS2::winvar(v1, tr)/(h1 *
                                                     (h1 - 1))
    q2 <- (length(v2) - 1) * WRS2::winvar(v2, tr)/(h2 *
                                                     (h2 - 1))
    dft <- (q1 + q2)^2/((q1^2/(h1 - 1)) + (q2^2/(h2 -
                                                   1)))
    crit <- stats::qt(1 - alpha/2, dft)
    dif1 <- mean(v1, tr) - mean(v2, tr) - eiU
    dif2 <- mean(v1, tr) - mean(v2, tr) - eiL
    t1 <- dif1/sqrt(q1 + q2)
    t2 <- dif2/sqrt(q1 + q2)
    probt1 <- stats::pt(t1, dft)
    probt2 <- 1 - stats::pt(t2, dft)
    ifelse(probt1 <= alpha & probt2 <= alpha,
           decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval can be rejected. A negligible difference in population means can be concluded. Be sure to interpret the magnitude of the effect size.",
           decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval cannot be rejected. A negligible difference in means cannot be concluded. Be sure to interpret the magnitude of the effect size.")
    effsize_raw<-mean(v1,tr=tr)-mean(v2,tr=tr)
    effsize_d <- .642*(mean(v1,tr=tr) - mean(v2,tr=tr))/sqrt((WRS2::winvar(v1,tr=.2)+WRS2::winvar(v2,tr=.2))/2)
    ifelse(sign(mean(v1,tr=tr)-mean(v2,tr=tr))==sign(eiL),
           ein<-eiL,ein<-eiU)
    effsize_pd<-(mean(v1,tr=tr)-mean(v2,tr=tr))/abs(ein)
    bsraw<-numeric(nboot)
    bsd<-numeric(nboot)
    bspd<-numeric(nboot)
    for (i in 1:nboot) {
      bsx<-sample(v1, length(v1),replace=TRUE)
      bsy<-sample(v2,length(v2),replace=TRUE)
      bsraw[i]<-mean(bsx,tr=tr)-mean(bsy,tr=tr)
      bsd[i]<-.642*(mean(bsx,tr=tr)-mean(bsy,tr=tr))/
        sqrt((WRS2::winvar(bsx,tr=.2)+WRS2::winvar(bsy,tr=.2))/2)
      ifelse(sign(mean(bsx,tr=tr)-mean(bsy,tr=tr))==sign(eiL),
             ein2<-eiL,ein2<-eiU)
      bspd[i]<-(mean(bsx,tr=tr)-mean(bsy,tr=tr))/abs(ein2)
    }
    ci_raw2<-stats::quantile(bsraw,probs=c(alpha/2,1-alpha/2))
    ci_raw<-stats::quantile(bsraw,probs=c(alpha,1-alpha))
    ci_d<-stats::quantile(bsd,probs=c(alpha,1-alpha))
    ci_pd<-stats::quantile(bspd,probs=c(alpha/2,1-alpha/2))
    title <- "Schuirmann-Yuen Test of the Equivalence of Two Independent Groups"
  }
  ret <- data.frame(lev1 = lev1,
                    lev2 = lev2,
                    meanx = mean(v1),
                    meany = mean(v2),
                    trmeanx = mean(v1, tr),
                    trmeany = mean(v2, tr),
                    sdx = stats::sd(v1),
                    sdy = stats::sd(v2),
                    madx = stats::mad(v1),
                    mady = stats::mad(v2),
                    eiL = eiL,
                    eiU = eiU,
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
                    oe="Mean Difference",
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
  cat("Levels of the IV: ", x$lev1, ", ", x$lev2, "\n", sep="")
  cat("Group Means: ", x$meanx, ", ", x$meany, "\n", sep="")
  cat("Group Trimmed Means: ", x$trmeanx, ", ", x$trmeany, "\n", sep="")
  cat("Group SDs: ", x$sdx, ", ", x$sdy, "\n", sep="")
  cat("Group MADs: ", x$madx, ", ", x$mady, "\n\n", sep="")
  cat("**********************\n\n")
  cat("Standardized Mean Difference (SMD):", x$effsized, "\n")
  cat(100*(1-2*x$alpha), "% CI for SMD: ", "(", x$cild, ", ", x$ciud, ")", "\n\n", sep="")
  cat("**********************\n\n")
  cat("Equivalence Interval:","Lower =", x$eiL, ",", "Upper =", x$eiU, "\n\n")
  cat("**********************\n\n")
  cat("Mean Difference (MD):", x$effsizeraw, "\n")
  cat(100*(1-2*x$alpha), "% CI for MD: ", "(", x$cilraw, ", ", x$ciuraw, ")", "\n", sep="")
  cat(100*(1-x$alpha), "% CI for MD: ", "(", x$cilraw2, ", ", x$ciuraw2, ")", "\n\n", sep="")
  cat("**********************\n\n")
  cat("Proportional Distance (PD):", x$effsizepd, "\n")
  cat(100*(1-x$alpha), "% CI for PD: ", "(", x$cilpd, ", ", x$ciupd, ")", "\n\n", sep="")
  cat("**********************\n\n")
  cat("TOST Test Statistics:", "\n\n", "Ho: mu1-mu2>=eiU:", "\n", "t = ", x$t1,
      " (df = ", x$df1, ")", ", p = ", x$pval1, "\n\n", "Ho: mu1-mu2<=eiL:", "\n", "t = ", x$t2,
      " (df = ", x$df1, ")", ", p = ", x$pval2, "\n\n", sep="")
  cat("**********************\n\n")
  cat("NHST Decision:", "\n", x$decis, "\n\n", sep="")
  cat("**********************\n\n")
  
  if (x$pl == TRUE) {
    neg.pd(effect = x$effsizeraw,
           PD = x$effsizepd,
           eil=x$eiL,
           eiu=x$eiU,
           PDcil=x$cilpd,
           PDciu=x$ciupd,
           cil=x$cilraw,
           ciu=x$ciuraw,
           Elevel=100*(1-2*x$alpha),
           Plevel=100*(1-x$alpha),
           save = x$save,
           oe=x$oe)
  }
}

