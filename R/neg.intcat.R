#' @title Test for Negligible Interaction between Two Categorical Variables
#' with a Continuous Outcome
#'
#' @description This function allows researchers to test whether the
#' interaction effect among two categorical independent variables, with a
#' continuous outcome variable, is negligible.
#'
#'
#' @param iv1 Levels of the first independent variable
#' @param iv2 Levels of the second independent variable
#' @param dv Score on the continuous dependent/outcome variable
#' @param neiL Lower bound of the negligible effect interval
#' @param neiU Upper bound of the negligible effect interval
#' @param nboot Number of bootstrap samples for calculating CIs
#' @param alpha Nominal Type I Error rate
#' @param data Dataset containing iv1, iv2 and dv
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
#' @details This function allows researchers to test whether the
#' interaction effect among two categorical independent variables, with a
#' continuous outcome variable, is negligible. In this case, 'negligible'
#' represents the minimum meaningful interaction effect.
#'
#' This test uses an intersection union approach, where a decision regarding
#' the omnibus interaction effect is inferred from the decision regarding all
#' simple (2 x 2) interaction effects; in other words, if all simple interaction
#' effects are deemed negligible, then the omnibus interaction is also
#' deemed negligible.
#'
#' The test also uses the percentile bootstrap to determine confidence
#' intervals, an approach that has been found to be robust to violations
#' of normality and variance homogeneity.
#'
#' See Cribbie, R. A., Ragoonanan, C., & Counsell, A. (2016). Testing for
#' negligible interaction: A coherent and robust approach. British Journal of
#' Mathematical and Statistical Psychology, 69, 159-174.
#'
#' @author Rob Cribbie \email{cribbie@@yorku.ca}
#'
#' @export neg.intcat
#'
#' @examples
#' outcome<-rnorm(60,mean=50,sd=10)
#' iv_1<-rep(c("male","female"),each=30)
#' iv_2<-rep(c("young","middle","old"),each=10,times=2)
#' d<-data.frame(iv_1,iv_2,outcome)
#' neg.intcat(iv1=iv_1,iv2=iv_2,dv=outcome,neiL=-15,neiU=15,nboot=10)
#' neg.intcat(iv1=iv_1,iv2=iv_2,dv=outcome,neiL=-15,neiU=15,nboot=10,data=d)

neg.intcat <- function(iv1 = NULL, iv2 = NULL, dv = NULL,
                            neiL, neiU, nboot = 50, alpha = 0.05,
                            data=NULL) {
  if (!is.null(data)) {
    iv1 <- deparse(substitute(iv1))
    iv2 <- deparse(substitute(iv2))
    dv <- deparse(substitute(dv))
  }
  if (!is.null(iv1) & !is.null(iv2) &
      !is.null(dv) & !is.null(data)) {
    iv1<-as.factor(data[[iv1]])
    iv2<-as.factor(data[[iv2]])
    dv<-as.numeric(data[[dv]])
    dat<- data.frame(iv1,iv2,dv)
    dat<- dat[stats::complete.cases(dat),]
  }
  if (is.null(data)) {
    iv1<-as.factor(iv1)
    iv2<-as.factor(iv2)
    dv<-as.numeric(dv)
    dat<- data.frame(iv1,iv2,dv)
    dat<- dat[stats::complete.cases(dat),]
  }

  a <- length(levels(dat$iv1))
  b <- length(levels(dat$iv2))
  iv1_lev <- levels(dat$iv1)
  iv2_lev <- levels(dat$iv2)

  S<-(a*(a-1)/2)*(b*(b-1)/2)
  div<-2*S/3

  means<-tapply(dat$dv, dat$iv1:dat$iv2, mean)
  rownames(means)<-NULL

  ##################################################################
  ######################### means table ############################

  grouped_data <- dplyr::group_by(dat, iv1, iv2)
  summarized_data <- dplyr::summarise(grouped_data, mean_value = mean(dv, na.rm = TRUE), .groups = "keep")
  ungrouped_data <- dplyr::ungroup(summarized_data)
  means_data <- ungrouped_data

  means_table <- stats::xtabs(
  formula = mean_value ~ iv1 + iv2,
  data = means_data
)

  means_table <- as.data.frame.matrix(means_table)
  ##################################################################


  dat$id<-factor(1:nrow(dat))
  gen_eta2 <- suppressMessages(ez::ezANOVA(dat,between=c(iv1,iv2),dv=dv,wid=id,
                          white.adjust = FALSE)$ANOVA$ges[3])
  dat$iv1 <- as.numeric(dat$iv1)
  dat$iv2 <- as.numeric(dat$iv2)
  bootb<-matrix(NA,nrow=nboot,ncol=S)
  bootres<-matrix(NA,nrow=1,ncol=S)
  bootresmcp<-matrix(NA,nrow=1,ncol=S)
  omn_ci<-matrix(NA,nrow=S,ncol=2)
  ic_ci<-matrix(NA,nrow=S,ncol=2)
  omn_res<-matrix(NA,nrow=S,ncol=3)
  ic_res<-matrix(NA,nrow=S,ncol=3)
  ## Intersection-Union/Bootstrap Method
  q<-0
  for (j in 1:(a-1)) {
    for (k in (j+1):a) {
      for (l in 1:(b-1)) {
        for (m in (l+1):b) {
          q<-q+1
          for (e in 1:nboot) {
            dvbdat<-dplyr::slice_sample(dat,n=length(dat$dv),replace=TRUE)
            bootb[e,q]<-(mean(dvbdat$dv[dvbdat$iv1==j & dvbdat$iv2==l])-
                           mean(dvbdat$dv[dvbdat$iv1==j & dvbdat$iv2==m]))-
              (mean(dvbdat$dv[dvbdat$iv1==k & dvbdat$iv2==l])-
                 mean(dvbdat$dv[dvbdat$iv1==k & dvbdat$iv2==m]))
          }
  ifelse (q == 1, res<-c(iv1_lev[j],iv1_lev[k],iv2_lev[l],iv2_lev[m]),
                res<-rbind(res,c(iv1_lev[j],iv1_lev[k],iv2_lev[l],iv2_lev[m])))
  ifelse (stats::quantile(bootb[,q],alpha)>neiL &
            stats::quantile(bootb[,q],1-alpha)<neiU,
          bootres[1,q]<-1,bootres[1,q]<-0)
  ifelse(bootres[1,q]==1, omn_res[q,1]<-"Sig",omn_res[q,1]<-"NS")
  ifelse (stats::quantile(bootb[,q],alpha/div)>neiL &
            stats::quantile(bootb[,q],1-alpha/div)<neiU,
          bootresmcp[1,q]<-1,bootresmcp[1,q]<-0)
  ifelse(bootresmcp[1,q]==1, ic_res[q,1]<-"Sig",ic_res[q,1]<-"NS")

  omn_ci[q,1]<-stats::quantile(bootb[,q],alpha)
  omn_ci[q,2]<-stats::quantile(bootb[,q],1-alpha)
  omn_res[q,2]<-round(omn_ci[q,1],3)
  omn_res[q,3]<-round(omn_ci[q,2],3)
  ic_ci[q,1]<-stats::quantile(bootb[,q],alpha/div)
  ic_ci[q,2]<-stats::quantile(bootb[,q],1-(alpha/div))
  ic_res[q,2]<-round(ic_ci[q,1],3)
  ic_res[q,3]<-round(ic_ci[q,2],3)
        }
      }
    }
  }
  res_nomcp <- cbind(res,omn_res)
  colnames(res_nomcp)<-c("iv1a","iv1b","iv2a","iv2b","Neg_Int","CI_lower","CI_upper")
  res_mcp <- cbind(res,ic_res)
  colnames(res_mcp)<-c("iv1a","iv1b","iv2a","iv2b","Neg_Int","CI_lower","CI_upper")

  ifelse(all(bootres[1,]==1),omnibus<-1,omnibus<-0)
  ifelse(omnibus==1, outomn<-"The bootstrap-based intersection union test of negligible interaction among two categorical varaibles is statistically significant; negligible interaction can be concluded. In other words, all interaction contrast negligible effect tests were statistically significant.",
         outomn<-"The bootstrap-based intersection union test of negligible interaction among two categorical varaibles is NOT statistically significant; negligible interaction cannot be concluded. In other words, NOT all interaction contrast negligible effect tests were statistically significant.")

  iv1levs<-rep(iv1_lev, each=b)
  iv2levs<-rep(iv2_lev, times=a)
  cells<-a*b
  title<-"Bootstrap-based Intersection Union Test of Negligible Interaction Among Two Categorical Variables"

  ret <- data.frame(title=title,
                    iv1 = iv1,
                    iv2 = iv2,
                    dv = dv,
                    gen_eta2 = gen_eta2,
                    iv1levs = iv1levs,
                    iv2levs = iv2levs,
                    means = means,
                    cells = cells,
                    leviv1 = a,
                    leviv2 = b,
                    omnibus = omnibus,
                    neiL = neiL,
                    neiU = neiU,
                    nomcp_iv1a<-res_nomcp[,1],
                    nomcp_iv1b<-res_nomcp[,2],
                    nomcp_iv2a<-res_nomcp[,3],
                    nomcp_iv2b<-res_nomcp[,4],
                    nomcp_nhst<-res_nomcp[,5],
                    nomcp_cil<-res_nomcp[,6],
                    nomcp_ciu<-res_nomcp[,7],
                    mcp_iv1a<-res_mcp[,1],
                    mcp_iv1b<-res_mcp[,2],
                    mcp_iv2a<-res_mcp[,3],
                    mcp_iv2b<-res_mcp[,4],
                    mcp_nhst<-res_mcp[,5],
                    mcp_cil<-res_mcp[,6],
                    mcp_ciu<-res_mcp[,7],
                    outomn<-outomn,
                    q<-q
                    )
  ret <- list(ret=ret, means_table=as.data.frame(means_table))
  class(ret) <- "neg.intcat"
  return(ret)
}

#' @rdname neg.intcat
#' @param x object of class \code{neg.twointcat}
#' @param ... extra arguments
#' @export
#'

print.neg.intcat <- function(x, ...) {
  cat("\n\n")
  cat("---- Test of Negligible Interaction Among Two Categorical Variables----\n\n")
  cat("\n", x$ret$title[1], "\n\n", sep="")
  cat("**********************************************\n\n")
  cat("\n", "Cell Means for the Two-Way Design", "\n\n", sep="")
  ##################################################################
  ######################### means table ############################
  print(round(x$means_table,3))
  ##################################################################
  cat("\n\n", "Generalized Eta-Squared for the Interaction:", "\n\n", sep="")
  cat( round(x$ret$gen_eta2[1],3), sep="", "\n")
  cat("\n\n", "Negligible Effect Interval:", "\n\n", sep="")
  cat("{", x$ret$neiL[1], ", ", x$ret$neiU[1], "}", "\n", sep="")
  cat("\n\n", "Omnibus Test Result (Negligible Interaction):","\n\n", x$ret$outomn[1], sep="")
  cat("\n\n\n", "Negligible Interaction Contrast Results (No Multiplicity Control):", "\n", sep="")
  cat("\n", "Each row is one interaction contrast involving the 2 x 2 matrix made up of IV1 (iv1a, iv1b) and IV2 (iv2a, iv2b)", "\n\n", sep="")
  xx<-data.frame(x$ret$nomcp_iv1a[1:x$ret$q[1]],x$ret$nomcp_iv1b[1:x$ret$q[1]],
           x$ret$nomcp_iv2a[1:x$ret$q[1]],x$ret$nomcp_iv2b[1:x$ret$q[1]],
           x$ret$nomcp_nhst[1:x$ret$q[1]],x$ret$nomcp_cil[1:x$ret$q[1]],
           x$ret$nomcp_ciu[1:x$ret$q[1]])
  colnames(xx)<-c("iv1a","iv1b","iv2a","iv2b","Sig?","CI_l","CI_u")
  print(xx)
  cat("\n\n", "Negligible Interaction Contrast Results (Multiplicity Control):", "\n", sep="")
  cat("\n", "See Cribbie, Ragoonanan & Cousell (2016) for details regarding the multiplicity control", "\n", sep="")
  cat("\n", "Each row is one interaction contrast involving the 2 x 2 matrix made up of IV1 (iv1a, iv1b) and IV2 (iv2a, iv2b)", "\n\n", sep="")
  yy<-data.frame(x$ret$mcp_iv1a[1:x$ret$q[1]],x$ret$mcp_iv1b[1:x$ret$q[1]],
           x$ret$mcp_iv2a[1:x$ret$q[1]],x$ret$mcp_iv2b[1:x$ret$q[1]],
           x$ret$mcp_nhst[1:x$ret$q[1]],x$ret$mcp_cil[1:x$ret$q[1]],
           x$ret$mcp_ciu[1:x$ret$q[1]])
  colnames(yy)<-c("iv1a","iv1b","iv2a","iv2b","Sig?","CI_l","CI_u")
  print(yy)
}

