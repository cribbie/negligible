#' @title Negligible Effect Test on the Difference between the Means of Dependent Populations
#'
#' @description This function allows researchers to test whether the difference
#' between the means of two dependent populations is negligible, where
#' negligible represents the smallest meaningful effect size (MMES)
#' @param var1 Data for Group 1 (if outcome, group and ID are omitted)
#' @param var2 Data for Group 2 (if outcome, group and ID are omitted)
#' @param outcome Dependent Variable (if var1 and var2 are omitted)
#' @param group Dichotomous Predictor/Independent Variable (if var1 and var2 are omitted)
#' @param ID participant ID (if var1 and var2 are omitted)
#' @param eil Lower Bound of the Equivalence Interval
#' @param eiu Upper Bound of the Equivalence Interval
#' @param normality Are the population variances (and hence the residuals) assumed to be normally distributed?
#' @param nboot Number of bootstrap samples for calculating CIs
#' @param alpha Nominal Type I Error rate
#' @param plot Should a plot of the results be produced?
#' @param saveplot Should the plot be saved?
#' @param data Dataset containing var1/var2 or outcome/group/ID
#' @param seed Seed number
#' @param ... Extra arguments
#'
#' @return A \code{list} including the following:
#' \itemize{
#'   \item \code{meanx} Sample mean of the first population/group.
#'   \item \code{meany} Sample mean of the second population/group.
#'   \item \code{medx} Sample median of the first population/group.
#'   \item \code{medy} Sample median second population/group.
#'   \item \code{sdx} Sample standard deviation of the first population/group.
#'   \item \code{sdy} Sample standard deviation of the second population/group.
#'   \item \code{madx} Sample median absolute deviation of the first population/group.
#'   \item \code{mady} Sample median absolute deviation of the second population/group.
#'   \item \code{eil} Lower bound of the negligible effect (equivalence) interval.
#'   \item \code{eiu} Upper bound of the negligible effect (equivalence) interval.
#'   \item \code{effsizeraw} Simple difference in the means (or medians if normality = FALSE)
#'   \item \code{cilraw2} Lower bound of the 1-alpha CI for the raw mean difference.
#'   \item \code{ciuraw2} Upper bound of the 1-alpha CI for the raw mean difference.
#'   \item \code{cilraw} Lower bound of the 1-2*alpha CI for the raw mean difference.
#'   \item \code{ciuraw} Upper bound of the 1-2*alpha CI for the raw mean difference.
#'   \item \code{effsized} Standardized mean (or median if normality = FALSE) difference.
#'   \item \code{cild} Lower bound of the 1-alpha CI for the standardized mean (or median if normality = FALSE) difference.
#'   \item \code{ciud} Upper bound of the 1-alpha CI for the standardized mean (or median if normality = FALSE) difference.
#'   \item \code{effsizepd} Proportional distance statistic.
#'   \item \code{cilpd} Lower bound of the 1-alpha CI for the proportional distance statistic.
#'   \item \code{ciupd} Upper bound of the 1-alpha CI for the proportional distance statistic.
#'   \item \code{t1} First t-statistic from the TOST procedure.
#'   \item \code{t1} Second t-statistic from the TOST procedure.
#'   \item \code{df1} Degrees of freedom for the first t-statistic from the TOST procedure.
#'   \item \code{df2} Degrees of freedom for the second t-statistic from the TOST procedure.
#'   \item \code{pval1} p value associated with the first t-statistic from the TOST procedure.
#'   \item \code{pval2} p value associated with the second t-statistic from the TOST procedure.
#'   \item \code{alpha} Nominal Type I error rate
#'   \item \code{seed} Seed number
#' }
#' @export
#' @details This function evaluates whether the difference in the means of 2 dependent populations can be considered negligible (i.e., the population means can be considered equivalent).
#'
#' The user specifies either the data associated with the first and second groups/populations (var1, var2, both should be continuous) or specifies the Indepedent Variable/Predictor (group, should be a factor) and the Dependent Variable (outcome, should be continuous). A 'data' statement can be used if the variables are stored in an R dataset.
#'
#' The user must also specify the lower and upper bounds of the negligible effect (equivalence) interval. These are specified in the original units of the outcome variable.
#'
#' @author Rob Cribbie \email{cribbie@@yorku.ca}
#'   Naomi Martinez Gutierrez \email{naomimg@@yorku.ca}
#' @export neg.paired
#'
#' @examples
#' #wide format
#' ID<-rep(1:20)
#' control<-rnorm(20)
#' intervention<-rnorm(20)
#' d<-data.frame(ID, control, intervention)
#' head(d)
#' neg.paired(var1=control,var2=intervention,eil=-1,eiu=1,plot=TRUE,
#'            data=d)
#' neg.paired(var1=d$control,var2=d$intervention,eil=-1,eiu=1,plot=TRUE)
#' neg.paired(var1=d$control,var2=d$intervention,eil=-1,eiu=1,normality=FALSE,
#'            plot=TRUE)
#'
#' \dontrun{
#' #long format
#' sample1<-sample(1:20, 20, replace=FALSE)
#' sample2<-sample(1:20, 20, replace=FALSE)
#' ID<-c(sample1, sample2)
#' group<-rep(c("control","intervention"),c(20,20))
#' outcome<-c(control,intervention)
#' d<-data.frame(ID,group,outcome)
#' neg.paired(outcome=outcome,group=group,ID=ID,eil=-1,eiu=1,plot=TRUE,data=d)
#' neg.paired(outcome=d$outcome,group=d$group,ID=d$ID,eil=-1,eiu=1,plot=TRUE)
#' neg.paired(outcome=d$outcome,group=d$group,ID=d$ID,eil=-1,eiu=1,plot=TRUE, normality=FALSE)
#'
#' #long format with multiple variables
#' sample1<-sample(1:20, 20, replace=FALSE)
#' sample2<-sample(1:20, 20, replace=FALSE)
#' ID<-c(sample1, sample2)
#' attendance<-sample(1:3, 20, replace=TRUE)
#' group<-rep(c("control","intervention"),c(20,20))
#' outcome<-c(control,intervention)
#' d<-data.frame(ID,group,outcome,attendance)
#' neg.paired(outcome=outcome,group=group,ID=ID,eil=-1,eiu=1,plot=TRUE,data=d)
#' neg.paired(outcome=d$outcome,group=d$group,ID=d$ID,eil=-1,eiu=1,plot=TRUE)
#'
#' #open a dataset
#' library(negligible)
#' d<-perfectionism
#' names(d)
#' head(d)
#' neg.paired(var1=atqpre.total,var2=atqpost.total,
#'            eil=-10,eiu=10,data=d)
#'
#' #Dataset with missing data
#' x<-rnorm(10)
#' x[c(3,6)]<-NA
#' y<-rnorm(10)
#' y[c(7)]<-NA
#' neg.paired(x,y,eil=-1,eiu=1, normality=FALSE)
#' }

neg.paired <- function(var1 = NULL, var2 = NULL,
                       outcome = NULL, group = NULL, ID = NULL,
                       eil, eiu, normality = TRUE,
                       nboot = 10000, alpha = 0.05,
                       plot = TRUE, saveplot = FALSE,
                       data=NULL, seed = NA,...) {
  if (!is.null(data)) {
    var1 <- deparse(substitute(var1))
    var2 <- deparse(substitute(var2))
    group <- deparse(substitute(group))
    outcome <- deparse(substitute(outcome))
    ID <- deparse(substitute(ID))
    if(var1=="NULL") {var1<-NULL}
    if(var2=="NULL") {var2<-NULL}
    if(group=="NULL") {group<-NULL}
    if(outcome=="NULL") {outcome<-NULL}
    if(ID=="NULL") {ID<-NULL}
  }
  
  if (is.null(var1) & is.null(var2) & !is.null(data)) {
    d<-data.frame(data$ID,data$group,data$outcome)
    names(d)[1] <- "ID"
    names(d)[2] <- "group"
    names(d)[3] <- "outcome"
    d <- d[stats::complete.cases(d),]
    d_wide <- tidyr::spread(d, group, outcome)
    var1<-as.numeric(unlist(d_wide[2]))
    var2<-as.numeric(unlist(d_wide[3]))
  }
  if (is.null(outcome) & is.null(group) & is.null(ID) & is.null(data)) {
    d <- data.frame(var1,var2)
    d <- d[stats::complete.cases(d),]
    var1 <- d$var1
    var2 <- d$var2
  }
  
  if (is.null(outcome) & is.null(group) & is.null(ID) & !is.null(data)) {
    var1 <- as.numeric(data[[var1]])
    var2 <- as.numeric(data[[var2]])
    d <- data.frame(var1,var2)
    d <- d[stats::complete.cases(d),]
    var1 <- d$var1
    var2 <- d$var2
  }
  
  if (is.null(var1) & is.null(var2) & is.null(data)) {
    d<-data.frame(ID,group,outcome)
    d <- d[stats::complete.cases(d),]
    d_wide <- tidyr::spread(d, group, outcome)
    var1<-as.numeric(unlist(d_wide[2]))
    var2<-as.numeric(unlist(d_wide[3]))
  }
  
  if(is.na(seed)){
    seed <- sample(.Random.seed[1], size = 1)
  } else {
    seed <- seed
  }
  set.seed(seed)
  
  if (normality==TRUE) {
    #TOST-P:
    r12<-stats::cor(var1,var2)
    sdif<-sqrt(stats::sd(var1)^2+stats::sd(var2)^2-2*r12*stats::sd(var1)*stats::sd(var2))
    se<-sdif/sqrt(length(var1)-1)
    dft<-length(var1) - 1
    t1<-(mean(var1) - mean(var2) - eiu)/se
    t2 <- (mean(var1) - mean(var2) - eil)/se
    probt1<-stats::pt(t1, dft, lower.tail=TRUE)
    probt2<-stats::pt(t2, dft, lower.tail=FALSE)
    ifelse(probt1 <= alpha & probt2 <= alpha,
           decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval can be rejected. Be sure to interpret the magnitude of the effect size.",
           decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval cannot be rejected. Be sure to interpret the magnitude of the effect size.")
    effsize_raw<-mean(var1)-mean(var2)
    difference = var1-var2
    effsize_d<-(mean(var1)-mean(var2))/stats::sd(difference)
    ifelse(sign(mean(var1)-mean(var2))==sign(eil),
           ein<-eil,ein<-eiu)
    effsize_pd<-(mean(var1)-mean(var2))/abs(ein)
    bsraw<-numeric(nboot)
    bsd<-numeric(nboot)
    bspd<-numeric(nboot)
    for (i in 1:nboot) {
      bsx<-sample(var1, length(var1),replace=TRUE)
      bsy<-sample(var2, length(var2),replace=TRUE)
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
    title <- "Two One-Sided Test of Equivalence for Paired-Samples (TOST-P)"
    title2 <- "Effect Size Measures use the Ordinary Mean and Standard Deviation"
    ret <- data.frame(meanx = mean(var1),
                      meany = mean(var2),
                      medx = stats::median(var1),
                      medy = stats::median(var2),
                      sdx = stats::sd(var1),
                      sdy = stats::sd(var2),
                      madx = stats::mad(var1),
                      mady = stats::mad(var2),
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
                      t1 = t1, #the NPAR results are z-scores
                      t2 = t2,
                      df1 = dft,
                      df2 = dft,
                      pval1 = round(probt1,5),
                      pval2 = round(probt2,5),
                      decis = decis,
                      alpha = alpha,
                      title1 = title,
                      title2 = title2,
                      pl = plot,
                      norm = normality,
                      save = saveplot,
                      oe = "Mean Difference",
                      seed = seed)
    class(ret) <- "neg.paired"
    return(ret)
  }
  if (normality == FALSE) {
    ##NPAR
    s1 <- var1 - var2 - eiu
    s2 <- var1 - var2 - eil
    dft <- length(var1)
    r1 <- rank(abs(s1))
    r2 <- rank(abs(s2))
    sr1 <- sum(r1[s1 < 0]) #sr1 is the absolute value of the sum of the negative ranks associated with s1.
    sr2 <- sum(r2[s2 > 0])  #sr2 is the absolute value of the sum of the positive ranks associated with s2.
    z1 <- (sr1 - (dft * (dft + 1)/4)) / (sqrt(dft * (dft + 1) * (2 * dft + 1) / 24))
    z2 <- (sr2 - (dft * (dft + 1)/4)) / (sqrt(dft * (dft + 1) * (2 * dft + 1) / 24))
    probt1<-stats::pnorm(z1, lower.tail=FALSE) #change to z later
    probt2<-stats::pnorm(z2, lower.tail=FALSE) #change to z later
    ifelse(probt1 <= alpha & probt2 <= alpha,
           decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval can be rejected. Be sure to interpret the magnitude of the effect size.",
           decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval cannot be rejected. Be sure to interpret the magnitude of the effect size.")
    effsize_raw<-stats::median(var1)-stats::median(var2)
    effsize_d<-(stats::median(var1)-stats::median(var2))/mean(c(stats::mad(var1),stats::mad(var2)))
    ifelse(sign(mean(var1)-mean(var2))==sign(eil),
           ein<-eil,ein<-eiu)
    effsize_pd<-(mean(var1)-mean(var2))/abs(ein)
    bsraw<-numeric(nboot)
    bsd<-numeric(nboot)
    bspd<-numeric(nboot)
    for (i in 1:nboot) {
      bsx<-sample(var1, length(var1),replace=TRUE)
      bsy<-sample(var2, length(var2),replace=TRUE)
      bsraw[i]<-stats::median(bsx)-stats::median(bsy)
      bsd[i]<-(stats::median(bsx)-stats::median(bsy))/
        mean(c(stats::mad(bsx),stats::mad(bsy)))
      ifelse(sign(mean(bsx)-mean(bsy))==sign(eil),
             ein2<-eil,ein2<-eiu)
      bspd[i]<-(mean(bsx)-mean(bsy))/abs(ein2)
    }
    ci_raw2<-stats::quantile(bsraw,probs=c(alpha/2,1-alpha/2))
    ci_raw<-stats::quantile(bsraw,probs=c(alpha,1-alpha))
    ci_d<-stats::quantile(bsd,probs=c(alpha,1-alpha))
    ci_pd<-stats::quantile(bspd,probs=c(alpha/2,1-alpha/2))
    title <- "Nonparametric Two One-sided Test of Equivalence for Paired Samples (NPAR)"
    title2 <- "Effect Size Measures use the Median and Median Absolute Deviation"
    ret <- data.frame(meanx = mean(var1),
                      meany = mean(var2),
                      medx = stats::median(var1),
                      medy = stats::median(var2),
                      sdx = stats::sd(var1),
                      sdy = stats::sd(var2),
                      madx = stats::mad(var1),
                      mady = stats::mad(var2),
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
                      z1 = z1,
                      z2 = z2,
                      df1 = dft,
                      df2 = dft,
                      pval1 = round(probt1,5),
                      pval2 = round(probt2,5),
                      decis = decis,
                      alpha = alpha,
                      title1 = title,
                      title2 = title2,
                      pl = plot,
                      norm = normality,
                      save = saveplot,
                      oe = "Mean Difference",
                      seed = seed)
    class(ret) <- "neg.paired"
    return(ret)
  }
}


print.neg.paired <- function(x, ...) {
  cat("---- Evaluating the Equivalence of the Central Tendencies of Paired Samples----\n\n")
  cat("Random Seed =", x$seed, "\n")
  cat("Test Statistic:", "\n", x$title1, "\n\n")
  cat("**********************\n\n")
  cat("Means: ", x$meanx, ", ", x$meany, "\n", sep="")
  cat("Medians: ", x$medx, ", ", x$medy, "\n", sep="")
  cat("SDs: ", x$sdx, ", ", x$sdy, "\n", sep="")
  cat("MADs: ", x$madx, ", ", x$mady, "\n\n", sep="")
  cat("**********************\n\n")
  cat("Equivalence Interval:","Lower =", x$eil, ",", "Upper =", x$eiu, "\n\n")
  cat("**********************\n\n")
  cat("Effect Sizes:", x$title2, "\n")
  cat("Standardized Mean Difference (SMD):", x$effsized, "\n")
  cat(100*(1-2*x$alpha), "% CI for SMD:", "(", x$cild, ", ", x$ciud, ")", "\n\n", sep="")
  if(x$norm == TRUE) {
    cat("Raw Mean Difference (MD):", x$effsizeraw, "\n")
    cat(100*(1-2*x$alpha), "% CI for MD:", "(", x$cilraw, ", ", x$ciuraw, ")", "\n", sep="")
    cat(100*(1-x$alpha), "% CI for MD:", "(", x$cilraw2, ", ", x$ciuraw2, ")", "\n\n", sep="")
  }
  if(x$norm == FALSE) {
    cat("Raw Median Difference (MD):", x$effsizeraw, "\n")
    cat(100*(1-x$alpha), "% CI for MD:", "(", x$cilraw2, ", ", x$ciuraw2, ")", "\n\n", sep="")
  }
  cat("**********************\n\n")
  cat("Proportional Distance (PD):", x$effsizepd, "\n")
  cat(100*(1-x$alpha), "% CI for PD:", "(", x$cilpd, ", ", x$ciupd, ")", "\n\n", sep="")
  cat("**********************\n\n")
  if(x$norm == TRUE) {
    cat("TOST Test Statistics:", "\n\n", "Ho: \u03BC1 - \u03BC2 \u2265 eiu:", "\n", "t = ", x$t1,
        " (df = ", x$df1, ")", ", p = ", x$pval1, "\n\n", "Ho: \u03BC1 - \u03BC2 \u2264 eil:", "\n", "t = ", x$t2,
        " (df = ", x$df1, ")", ", p = ", x$pval2, "\n\n", sep="")
  }
  else {
    cat("TOST Test Statistics:", "\n\n", "Ho: \u03BC1 - \u03BC2 \u2265 eiu:", "\n", "z = ", x$z1,
        " (df = ", x$df1, ")", ", p = ", x$pval1, "\n\n", "Ho: \u03BC1 - \u03BC2 \u2264 eil:", "\n", "z = ", x$z2,
        " (df = ", x$df1, ")", ", p = ", x$pval2, "\n\n", sep="")
  }
  cat("**********************\n\n")
  cat("NHST Decision:", "\n", x$decis, "\n\n", sep="")
  cat("**********************\n\n")
  
  if (x$pl == TRUE) {
    neg.pd(effect = x$effsizeraw,
           PD = x$effsizepd,
           eil=x$eil,
           eiu=x$eiu,
           PDcil=x$cilpd,
           PDciu=x$ciupd,
           cil=x$cilraw,
           ciu=x$ciuraw,
           Elevel=100*(1-2*x$alpha),
           Plevel=100*(1-x$alpha),
           save = x$save,
           oe = x$oe)
  }
}

#wide format
ID<-rep(1:20)
control<-rnorm(20)
intervention<-rnorm(20)
d<-data.frame(ID, control, intervention)
head(d)
neg.paired(var1=control,var2=intervention,eil=-1,eiu=1,plot=TRUE,data=d)

