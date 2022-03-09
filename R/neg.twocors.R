#' Test for Evaluating Negligible Effects of Two Independent or Dependent Correlation Coefficients: Based on Counsell & Cribbie (2015)
#'
#' This function tests whether two correlation coefficients can be considered equivalent according to a predefined interval (i.e., SESOI/MMES/delta) based on the Anderson-Hauck (1983) test of equivalence
#'
#' @aliases neg.twocors
#' @param r1v1 the name of the 1st variable included in the 1st correlation coefficient (r1, variable 1)
#' @param r1v2 the name of the 2nd variable included in the 1st correlation coefficient (r1, variable 2)
#' @param r2v1 the name of the 1st variable included in the 2nd correlation coefficient (r2, variable 1)
#' @param r2v2 the name of the 2nd variable included in the 2st correlation coefficient (r2, variable 2)
#' @param data a data.frame or matrix which includes the variables in r1 and r2
#' @param eiu upper bound of the equivalence interval measured as the largest difference between the two correlations for which the two coefficients would still be considered equivalent
#' @param eil lower bound of the equivalence interval measured as the largest difference between the two correlations for which the two coefficients would still be considered equivalent
#' @param alpha desired alpha level, defualt is .05
#' @param r1 entered 1st correlation coefficient manually, without a dataset
#' @param r2 entered 2nd correlation coefficient manually, without a dataset
#' @param n1 entered sample size associated with r1 manually, without a dataset
#' @param n2 entered sample size associated with r2 manually, without a dataset
#' @param dep are the correlation coefficients dependent (overlapping)?
#' @param r3 if the correlation coefficients are dependent and no datasets were entered, specify the correlation between the two, non-intersecting variables (e.g. if r1 = r12 and r2 = r13, then r3 = r23)
#' @param test 'AH' is the default based on recommendation in Counsell & Cribbie (2015), 'TOST' is an additional (albeit, more conservative) option.
#' @param bootstrap logical, default is TRUE, incorporating bootstrapping when calculating regression coefficients, SE, and CIs
#' @param nboot 1000 is the default. indicate if other number of bootstrapping iterations is desired
#' @param plots logical, plotting the results. TRUE is set as default
#' @param saveplots FALSE for no, "png" and "jpeg" for different formats
#' @param seed to reproduce previous analyses using bootstrapping, the user can set their seed of choice
#' @param ... additional arguments to be passed
#' @return returns a \code{list} containing each analysis and their respective statistics
#'   and decision
#'
#' @author Rob Cribbie \email{cribbie@@yorku.ca} and
#'   Alyssa Counsell \email{a.counsell@@ryerson.ca}
#' @export neg.twocors
#' @examples
#' \dontrun{
#' # Negligible difference between two correlation coefficients
#' # Equivalence interval: -.15 to .15
#' v1a<-stats::rnorm(10)
#' v2a<-stats::rnorm(10)
#' v1b <- stats::rnorm(10)
#' v2b <- stats::rnorm(10)
#' dat<-data.frame(v1a, v2a, v1b, v2b)
#' (r1 <- stats::cor(v1a, v2a))
#' (r2 <- stats::cor(v1b, v2b))
#' (r1-r2)
#' # dataset available (independent correlation coefficients):
#' neg.twocors(r1v1=v1a,r1v2=v2a,r2v1=v1b,r2v2=v2b,data=dat,eiu=.15,eil=-0.15)
#' # dataset available (dependent correlation coefficients):
#' neg.twocors(r1v1=v1a,r1v2=v1b,r2v1=v1a,r2v2=v2b,data=dat,eiu=.15,eil=-0.15,dep=TRUE)
#' # using the TOST procedure
#'  neg.twocors(r1v1=v1a,r1v2=v2a,r2v1=v1b,r2v2=v2b,data=dat,eiu=.15,eil=-0.15,test="TOST")
#' # without full data:
#'  neg.twocors(r1=.01,n1=300,r2=.015,n2=300,eil=-.1,eiu=.1)
#' # end
#'}


neg.twocors <- function(data=NULL, r1v1=NULL, r1v2=NULL, r2v1=NULL, r2v2=NULL, # specific arguments when data is available
                   r1=NULL, n1=NULL, r2=NULL, n2=NULL, # specific arguments when data is NOT available
                   dep=FALSE, r3=NA, # arguments related to dependent correlations
                   test = "AH", eiu, eil, alpha=.05, bootstrap=TRUE, nboot=1000, seed=NA, # testing-related arguments
                   plots=TRUE, saveplots=FALSE,...) # graphics-related arguments
  {

  ################################  delta / SESOI / Equivalence Interval ##########################
  if (is.null(eiu)) {
    stop("please enter the upper bound of the equivalence interval using the eiu= argument")}

  if (is.null(eil)) {
    stop("please enter the lower bound of the equivalence interval using the eil= argument")}

  if (sign(eiu) == sign(eil)) {
    stop("The equivalence interval must include zero")
  }
  if (eiu < 0 & eil > 0) {
    temp.eiu <- eil
    temp.eil <- eiu

    eiu <- temp.eiu
    eil <- temp.eil
  }


# FULL DATA SECTION, WHEN DATA INSERTED

          if(!is.null(data)) {
            withdata <- TRUE

            r1v1 <- deparse(substitute(r1v1))
            r1v2 <- deparse(substitute(r1v2))
            r2v1 <- deparse(substitute(r2v1))
            r2v2 <- deparse(substitute(r2v2))
            if(r1v1=="NULL") {r1v1<-NULL}
            if(r1v2=="NULL") {r1v2<-NULL}
            if(r2v1=="NULL") {r2v1<-NULL}
            if(r2v2=="NULL") {r2v2<-NULL}

              if(is.null(r1v1) | is.null(r1v2)){
              stop("Please choose two variables from your data to include in r1 caluclation")}
              if(is.null(r2v1) | is.null(r2v2)){
               stop("Please choose two variables from your data to include in r2 caluclation")}

            if(!(r1v1 %in% colnames(data))){
              stop(r1v1, " was not found in your data")
            } else if(!(r1v2 %in% colnames(data))){
              stop(r1v2, " was not found in your data")
            } else if(!(r2v1 %in% colnames(data))){
              stop(r2v1, " was not found in your data")
            } else if(!(r2v2 %in% colnames(data))){
              stop(r2v2, " was not found in your data")
            }

           # if(!is.numeric(r1v1)& !is.numeric(r1v2)& !is.numeric(r2v1) & !is.numeric(r2v2)){
              dr1v1 <- eval(substitute(data$r1v1), data)
              dr1v2 <- eval(substitute(data$r1v2), data)
              dr2v1 <- eval(substitute(data$r2v1), data)
              dr2v2 <- eval(substitute(data$r2v2), data)
            #}
              dat1 <- data.frame(dr1v1, dr1v2)
              names(dat1)[1] <- r1v1
              names(dat1)[2] <- r1v2
              dat1<- stats::na.omit(dat1)

              dat2 <- data.frame(dr2v1, dr2v2)
              names(dat2)[1] <- r2v1
              names(dat2)[2] <- r2v2
              dat2<- stats::na.omit(dat2)

              # r1 calculations
              if(!is.numeric(dr1v1) | !is.numeric(dr1v2)) {
               stop('Both variables for r1 must be numeric')}
              r1 <- stats::cor(dat1[,1], dat1[,2])
              n1 <- nrow(dat1)

              # r2 calculations
              if(!is.numeric(dr2v1) | !is.numeric(dr2v2)) {
                stop('Both variables for r2 must be numeric')}
              r2 <- stats::cor(dat2[,1], dat2[,2])
              n2 <- nrow(dat2)

              # r3 calculation
              shared.vars <- FALSE
              if (identical(dr1v1, dr2v1) | identical(r1v1,r2v1)) {
                r3 <- stats::cor(dr1v2, dr2v2)
                shared.vars <- TRUE}
              if (identical(dr1v2, dr2v2) | identical(r1v2,r2v2)){
                r3 <- stats::cor(dr1v1, dr2v1)
                shared.vars <- TRUE}
              if (identical(dr1v1, dr2v2)| identical(r1v1,r2v2)){
                r3 <- stats::cor(dr1v2, dr2v1)
                shared.vars <- TRUE}
              if (identical(dr1v2, dr2v1)| identical(r1v2,r2v1)){
                r3 <- stats::cor(dr1v1, dr2v2)
                shared.vars <- TRUE}

              if (dep==TRUE){
                if (shared.vars ==FALSE & is.na(r3)) {
                stop("The two correlation coefficients do not share a common variable and/or you did not specify r3")}
                }
                if(dep==FALSE & shared.vars==TRUE){
                  stop("you specified dep=FALSE (independent correlation coefficients), but looks like you are using a shared variable to calculate both r1 and r2, please adjust the dep argument to TRUE")
                }

            # RAW EFFECT SIZE
            rawdiff <-  r1 - r2
            se <- sqrt((((1 - r1^2)^2)/(n1 - 2)) + (((1 - r2^2)^2)/(n2 - 2))) # standard error for correlation coefficients difference, according to Kraatz (2007)
            temp.u <- rawdiff + stats::qnorm(alpha)*se # upper bound of confidence interval
            temp.l <- rawdiff - stats::qnorm(alpha)*se # lower bound of confidence interval
            temp.u.2a <- rawdiff + stats::qnorm(2*alpha)*se # upper bound of confidence interval
            temp.l.2a <- rawdiff - stats::qnorm(2*alpha)*se # lower bound of confidence interval
            if (temp.l < temp.u) {
                                u.ci.a <- temp.u
                                l.ci.a <- temp.l
                                u.ci.2a <- temp.u.2a
                                l.ci.2a <- temp.l.2a
            } else {
              u.ci.a <- temp.l
              l.ci.a <- temp.u
              u.ci.2a <- temp.l.2a
              l.ci.2a <- temp.u.2a
              } # if lower CI bound value is actually larger than upper value, switch definitions of lower and upper

# BOOTSTRAPPING SECTION

            # creating an empty matrix to soon hold our results
            rs.diffs<-numeric(nboot)
            propdis.rdif <-numeric(nboot)

            if(is.na(seed)){
              seed <- sample(.Random.seed[1], size = 1)
            } else {
              seed <- seed
            }
            set.seed(seed)
            for (i in 1:nboot) {
              temp.dat1 <- dplyr::sample_n(dat1,n1, replace = TRUE)
              temp.dat2 <- dplyr::sample_n(dat2,n2, replace = TRUE)
              temp.rl <- stats::cor(temp.dat1[,1], temp.dat1[,2])
              temp.r2 <- stats::cor(temp.dat2[,1], temp.dat2[,2])
              # raw diff
              rs.diffs[i]<- temp.rl - temp.r2
              # Proportional Distance
              ifelse(sign(rs.diffs[i])==sign(eiu), temp.EIc<-eiu, temp.EIc<-eil)
              propdis.rdif[i] <- rs.diffs[i]/abs(temp.EIc)
            } # end of bootstrapping For loop

# Proportional Distance and confidence intervals for PD
            ifelse(sign(rawdiff)==sign(eiu), EIc<-eiu, EIc<-eil)
            pd <- rawdiff/abs(EIc) # simple PD, without bootstrapping
            #pdbs <- mean(propdis.rdif) # PD from bootstrapping
            pd.se <- stats::sd(propdis.rdif)
            pd.u.ci <- stats::quantile(propdis.rdif,1-alpha/2)
            pd.l.ci <- stats::quantile(propdis.rdif,alpha/2)
            pd.u.ci.2a <- stats::quantile(propdis.rdif,1-alpha)
            pd.l.ci.2a <- stats::quantile(propdis.rdif,alpha)

            if (bootstrap==TRUE){
              #rawdiff.boot <- mean(rs.diffs) # the bootstrapping-generated rs difference
              ### se <- stats::sd(rs.diffs) ####### NEED TO CONFIRM
              u.ci.a <- stats::quantile(rs.diffs,1-alpha/2)
              l.ci.a <- stats::quantile(rs.diffs,alpha/2)
              u.ci.2a <- stats::quantile(rs.diffs,1-alpha)
              l.ci.2a <- stats::quantile(rs.diffs,alpha)
              #ifelse(sign(rawdiff.boot)==sign(eiu), EIc<-eiu, EIc<-eil)
              #pd <- rawdiff.boot/abs(EIc)
              #rawdiff <- rawdiff.boot
                }
              #}
            } # END OF FULL DATA SECTION

# NO DATA SECTION
  if(is.null(data) & !is.null(r1v1) & !is.null(r1v2) & !is.null(r2v1) & !is.null(r2v2)) {
    stop("if the data argument is not specified, please use the r1, n1, r2, and n2 arguments instead")
  }

    if(is.null(data) & is.null(r1v1) & is.null(r1v2) & is.null(r2v1) & is.null(r2v2)) {
     withdata <- FALSE
     r1v1 <- NA
     r1v2 <- NA
     r2v1 <- NA
     r2v2 <- NA

           if (!is.null(r1) & !is.null(n1)) {
                if (r1 >= 1 | r1 <= -1){ # making sure the correlation coefficient is valid, i.e., anywhere between -1 to 1
                  stop("Invalid correlation input for r1")
                  } else {
                    r1 <- as.numeric(r1)
                    n1 <- as.numeric(n1)}
          } else {
            stop("Please specify both n1 and r1")
          }

           if (!is.null(r2) & !is.null(n2)) {
      if (r2 >= 1 | r2 <= -1){ # making sure the correlation coefficient is valid, i.e., anywhere between -1 to 1
        stop("Invalid correlation input for r2")
        } else {
          r2 <- as.numeric(r2)
          n2 <- as.numeric(n2)}
    } else {
      stop("Please specify both n2 and r2")
    }

  # RAW EFFECT SIZE, NO DATA
  rawdiff <-  r1 - r2
  se <- sqrt((((1 - r1^2)^2)/(n1 - 2)) + (((1 - r2^2)^2)/(n2 - 2))) # standard error for correlation coefficients difference, according to Kraatz (2007)
  temp.u <- rawdiff + stats::qnorm(alpha)*se # upper bound of confidence interval
  temp.l <- rawdiff - stats::qnorm(alpha)*se # lower bound of confidence interval
  temp.u.2a <- rawdiff + stats::qnorm(2*alpha)*se # upper bound of confidence interval
  temp.l.2a <- rawdiff - stats::qnorm(2*alpha)*se # lower bound of confidence interval
  if (temp.l < temp.u) {
    u.ci.a <- temp.u
    l.ci.a <- temp.l
    u.ci.2a <- temp.u.2a
    l.ci.2a <- temp.l.2a
  } else {
    u.ci.a <- temp.l
    l.ci.a <- temp.u
    u.ci.2a <- temp.l.2a
    l.ci.2a <- temp.u.2a
  } # if lower CI bound value is actually larger than upper value, switch definitions of lower and upper

  # PD
  ifelse(sign(rawdiff)==sign(eiu), EIc<-eiu, EIc<-eil)
  pd <- rawdiff/abs(EIc)
  #pd.u.ci <- u.ci.a/abs(eiu) # because u.ci will be positive and l.ci will be negative
  #pd.l.ci <- l.ci.a/abs(eil)
  #pd.u.ci.2a <-  u.ci.2a/abs(eiu)
  #pd.l.ci.2a <- l.ci.2a/abs(eil)

  pd.u.ci <- u.ci.a/abs(EIc)
  pd.l.ci <- l.ci.a/abs(EIc)
  pd.u.ci.2a <-  u.ci.2a/abs(EIc)
  pd.l.ci.2a <- l.ci.2a/abs(EIc)

}

  # BOTH DATA AND NO DATA
  title <- "Test for Evaluating Negligible Difference Between Two Correlation Coefficients"

  if (dep == FALSE) { # Independent correlations tests
  subtitle <- "Comparison of Independent Correlation Coefficients"
  # AH-rho (Counsell & Cribbie, 2015)
  testID <- sprintf("AH-%s: Counsell-Cribbie Test for Comparing Two Independent Correlation Coefficients", "\u03C1")
  p.value <- stats::pnorm((abs(rawdiff) - eiu)/se) -
    stats::pnorm((-abs(rawdiff) - eiu)/se)
  r3 <- NA
  # KTOST-rho calculation as specified by to Counsell & Cribbie (2015)
  z1 <- (r1-r2-eil)/se
  z2 <- (r1-r2-eiu)/se
  p.value.1 <-stats::pnorm(z1, lower.tail=FALSE)
  p.value.2 <-stats::pnorm(z2, lower.tail=TRUE)
  testID.tost <- sprintf("KTOST-%s: Counsell-Cribbie Test for Comparing Two Independent Correlation Coefficients", "\u03C1")
  t1 <- NA
  t2 <- NA
  degfree <- NA
  }

  if (dep == TRUE) { # Dependent correlations tests
    if (withdata==FALSE & is.na(r3)){
        stop("you specified dep=TRUE (dependent correlation coefficients), please specify a value for r3, or include a dataset")
    }
  subtitle <- "Comparison of Dependent Correlation Coefficients"
  # AH-rho-D calculation as specified by to Counsell & Cribbie (2015)
  testID <- sprintf("AH-%s-D: Counsell-Cribbie Test for Comparing Two Dependent Correlation Coefficients", "\u03C1")
  R <- (1-r1^2-r2^2-r3^2)+(2*r1*r2*r3)
  N <- n1+n2
  inside.sqrt <- ((N-1)*(1+r3))/(2*R*(N-1)/(N-3) + ((1-r3)^3 * (r1+r2)^2)/4)
  pval.cmpnnt1 <- (abs(r1-r2)-eiu)*(sqrt(inside.sqrt))
  pval.cmpnnt2 <-  (-abs(r1-r2)-eiu)*(sqrt(inside.sqrt))
  p.value <- stats::pnorm(pval.cmpnnt1) - stats::pnorm(pval.cmpnnt2)
  # KTOST-rho-D calculation as specified by to Counsell & Cribbie (2015)
  t1 <- (r1-r2-eil)*(sqrt(inside.sqrt))
  t2 <- (r1-r2-eiu)*(sqrt(inside.sqrt))
  degfree <- N-3
  p.value.1 <-stats::pt(t1, degfree, lower.tail=FALSE)
  p.value.2 <-stats::pt(t2, degfree, lower.tail=TRUE)
  testID.tost <- sprintf("TOST-%s-D: Counsell-Cribbie Test for Comparing Two Dependent Correlation Coefficients", "\u03C1")
  z1 <- NA
  z2 <- NA
  }

# DECISION AH
  ifelse(p.value < alpha, decision <- "The null hypothesis that the difference between the two correlation coefficients is non-negligible (i.e., beyond the specified equivalence interval), can be rejected",
         decision <- "The null hypothesis that the difference between the two correlation coefficients is non-negligible (i.e., beyond the equivalence interval), was NOT rejected")
  # DECISION KTOST
  ifelse(p.value.1 < alpha & p.value.2 < alpha, decision.tost <- "The null hypothesis that the difference between the two correlation coefficients is non-negligible (i.e., beyond the equivalence interval), can be rejected",
         decision.tost <-"The null hypothesis that the difference between the two correlation coefficients is non-negligible (i.e., beyond the equivalence interval), was NOT rejected")

  # is the p value smaller than the accepted Type I error rate? If yes, then the null

# SUMMARY OF RESULTS
  res <- data.frame(title = title,
                    subtitle = subtitle,
                    test = test,
                    testID = testID,
                    testID.tost = testID.tost,
                    alpha = alpha,
                    r1v1 = r1v1,
                    r1v2 = r1v2,
                    r2v1 = r2v1,
                    r2v2 = r2v2,
                    r1 = r1,
                    n1 = n1,
                    r2 = r2,
                    n2 = n2,
                    r3 = r3,
                    rawdiff = rawdiff,
                    se = se,
                    dep = dep,
                    u.ci.a = u.ci.a,
                    l.ci.a = l.ci.a,
                    u.ci.2a = u.ci.2a,
                    l.ci.2a = l.ci.2a,
                    eiu = eiu,
                    eil = eil,
                    pd = pd,
                    pd.u.ci = pd.u.ci,
                    pd.l.ci = pd.l.ci,
                    pd.u.ci.2a = pd.u.ci.2a,
                    pd.l.ci.2a = pd.l.ci.2a,
                    EIc = EIc,
                    p.value = p.value,
                    p.value.1= p.value.1,
                    p.value.2 = p.value.2,
                    z1 = z1,
                    z2 = z2,
                    t1 = t1,
                    t2 = t2,
                    degfree = degfree,
                    decision = decision,
                    decision.tost = decision.tost,
                    perc.a = (1-alpha)*100,
                    perc.2a = (1-2*alpha)*100,
                    withdata = withdata,
                    seed = seed,
                    bootstrap = bootstrap,
                    nboot = nboot,
                    plots = plots,
                    saveplots = saveplots)
  class(res) <- "neg.twocors"
  return(res)
}
#' @rdname neg.twocors
#' @param x object of class \code{neg.twocors}
#' @export
#'

print.neg.twocors<- function(x,...) {
  cat("\n\n")
  cat(x$title, "\n\n")
  cat("***",x$subtitle, "***\n\n")

  if (x$withdata==TRUE){
    cat("Correlation coefficients:", "\n",
        " Variables: ", x$r1v1," & ",x$r1v2,", r1 = ",round(x$r1,3), "\n",
        " Variables: ", x$r2v1," & ",x$r2v2,", r2 = ",round(x$r2,3), "\n", sep = "")
    if(x$dep == TRUE){
      cat(" r3 = ",x$r3,"\n", sep = "")
    }

    cat("**********************\n\n")
    if (x$bootstrap == TRUE) {
      cat("Correlation coefficients' difference and confidence interval using ", x$nboot," bootstrap iterations (seed=",x$seed,"):","\n",
          " r1-r2 = ", round(x$rawdiff,3), ", ",x$perc.a, "% CI: [",round(x$l.ci.a,3),", ",round(x$u.ci.a,3),"]" ,"\n",
          " std. error = ", round(x$se,3), "\n", sep = "")
    } else {
      cat("Correlation coefficients' raw difference:","\n",
          " r1-r2 = ", round(x$rawdiff,3), ", ",x$perc.a, "% CI [",round(x$l.ci.a,3),", ",round(x$u.ci.a,3),"]" ,"\n",
          " std. error = ", round(x$se,3), "\n", sep = "")
    }
  } # end of withdata
  else {
    cat("Correlation coefficients:", "\n",
        " r1 = ",round(x$r1,3), "\n",
        " r2 = ",round(x$r2,3), "\n", sep = "")
    if(x$dep == TRUE){
       cat(" r3 = ",x$r3,"\n", sep = "")
      }
    cat("**********************\n\n")
    cat("Correlation coefficients' raw difference:","\n",
        " r1-r2 =", round(x$rawdiff,3), ", ",x$perc.a, "% CI [",round(x$l.ci.a,3),", ",round(x$u.ci.a,3),"]" ,"\n",
        " std. error = ", round(x$se,3), "\n", sep = "")
    }
  cat("**********************\n\n")

  # NHST results
  if (x$test == "AH"){
  ifelse(round(x$p.value,3) == 0, p.val <- "< 0.001", p.val <- paste("= ", round(x$p.value,3), sep = ""))
  cat(x$testID)
  cat("\nEquivalence Interval: ","Lower = ", x$eil, ", ", "Upper = ", x$eiu, "\n", sep = "")
  cat("p value ",p.val,"\n", sep = "")
  cat("NHST Decision: ", x$decision,"\n", sep = "")
  } else { # if KTOST
    cat(x$testID.tost)
    cat("\nEquivalence Interval: ","Lower = ", x$eil, ", ", "Upper = ", x$eiu, "\n", sep = "")
    if(x$test == "TOST" & x$bootstrap == TRUE){
      cat(x$perc.2a,"% bootstrap-generated CI for r1-r2: [",round(x$l.ci.2a,3),", ",round(x$u.ci.2a,3),"]\n", sep = "")
    }
    if(x$test == "TOST" & x$bootstrap == FALSE){
      cat(x$perc.2a,"% CI for r1-r2: [",round(x$l.ci.2a,3),", ",round(x$u.ci.2a,3),"]\n", sep = "")
    }
    if (x$dep == FALSE) {
      cat("Z1 value = ",round(x$z1,3),"\n", sep = "")
      cat("Z2 value = ",round(x$z2,3),"\n", sep = "")
      }
    if (x$dep == TRUE) {
      cat("t1 value = ",round(x$t1,3),"\n", sep = "")
      cat("t2 value = ",round(x$t2,3),"\n", sep = "")
      cat("with ",x$degfree," degrees of freedom\n", sep = "")
    }
    ifelse(round(x$p.value.1,3) == 0, p1.val <- "< 0.001", p1.val <- paste("= ", round(x$p.value.1,3), sep = ""))
    ifelse(round(x$p.value.2,3) == 0, p2.val <- "< 0.001", p2.val <- paste("= ", round(x$p.value.2,3), sep = ""))
    cat("p1 value ",p1.val,"\n", sep = "")
    cat("p2 value ",p2.val,"\n", sep = "", "\n")
    cat("NHST Decision: ", x$decision.tost,"\n", sep = "")
    cat("\n*Note that when using the KTOST-\u03C1 or TOST-\u03C1-D procedures, the null hypothesis of a non-negligible difference between the two population correlation coefficients is only rejected if BOTH p values are less than \u03B1.\n")
    }
    if(x$plots == TRUE) {
    neg.pd(effect=x$rawdiff, PD = x$pd, EIsign=x$EIc, PDcil=x$pd.l.ci, PDciu=x$pd.u.ci, cil=x$l.ci.2a,
           ciu=x$u.ci.2a, Elevel=100*(1-2*x$alpha), Plevel=100*(1-x$alpha), save = x$saveplots)
      if (x$test=="AH"){
      cat("\n*Note that NHST decisions using the AH-\u03C1 and AH-\u03C1-D procedures may not match KTOST-\u03C1 and TOST-\u03C1-D results or the Symmetric CI Approach at 100*(1-2\u03B1)% illustrated in the plot. \n")
        }
  }
  cat("**********************\n\n")
}
