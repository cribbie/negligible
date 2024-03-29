% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neg.reg.R
\name{neg.reg}
\alias{neg.reg}
\alias{print.neg.reg}
\title{Test for Evaluating Negligible Effect Between a Predictor and Outcome in a Multiple Regression Model}
\usage{
neg.reg(
  data = NULL,
  formula = NULL,
  predictor = NULL,
  b = NULL,
  se = NULL,
  nop = NULL,
  n = NULL,
  eil,
  eiu,
  alpha = 0.05,
  test = "AH",
  std = FALSE,
  bootstrap = TRUE,
  nboot = 1000,
  plots = TRUE,
  saveplots = FALSE,
  seed = NA,
  ...
)

\method{print}{neg.reg}(x, ...)
}
\arguments{
\item{data}{a data.frame or matrix which includes the variables considered in the regression model}

\item{formula}{an argument of the form y~x1+x2...xn which defines the regression model}

\item{predictor}{name of the variable/predictor upon which the test will be applied}

\item{b}{effect size of the regression coefficient of interest, can be in standardized or unstandardized units}

\item{se}{standard error associated with the above regression coefficient effect size, pay close attention to standardized vs. unstandardized}

\item{nop}{number of predictors (excluding intercept) in the regression model}

\item{n}{the sample size used in the regression analysis}

\item{eil}{lower bound of the equivalence interval measured in the same units as the regression coefficients (can be either standardized or unstandardized)}

\item{eiu}{upper bound of the equivalence interval measured in the same units as the regression coefficients (can be either standardized or unstandardized)}

\item{alpha}{desired alpha level, default is .05}

\item{test}{AH is the default based on recommendation in Alter & Counsell (2020), TOST is an additional option}

\item{std}{indicate if eil and eiu along with b (when dataset is not entered) are in standardized units}

\item{bootstrap}{logical, default is TRUE, incorporating bootstrapping when calculating regression coefficients, SE, and CIs}

\item{nboot}{1000 is the default. indicate if other number of bootstrapping iterations is desired}

\item{plots}{logical, plotting the results. TRUE is set as default}

\item{saveplots}{FALSE for no, "png" and "jpeg" for different formats}

\item{seed}{to reproduce previous analyses using bootstrapping, the user can set their seed of choice}

\item{...}{extra arguments}

\item{x}{object of class \code{neg.reg}}
}
\value{
A \code{list} containing the following:
\itemize{
  \item \code{formula} The regression model
  \item \code{effect} Specifying if effect size is in standardized or unstandardized units
  \item \code{test} Test type, i.e., Anderson-Hauck (AH) or Two One-Sided Tests (TOST)
  \item \code{t.value} t test statistic. If TOST was specified, only the smaller of the t values will be presented
  \item \code{df} Degrees of freedom associated with the test statistic
  \item \code{n} Sample size
  \item \code{p.value} p value associated with the test statistic. If TOST was specified, only the larger of the p values will be presented
  \item \code{eil} Lower bound of the negligible effect (equivalence) interval
  \item \code{eiu} Upper bound of the negligible effect (equivalence) interval
  \item \code{predictor} Variable name of the predictor in question
  \item \code{b} Effect size of the specified predictor
  \item \code{se} Standard error associated with the effect size point estimate (in the same units)
  \item \code{l.ci} Lower bound of the 1-alpha CI for the effect size
  \item \code{u.ci} Upper bound of the 1-alpha CI for the effect size
  \item \code{pd} Proportional distance
  \item \code{pd.l.ci} Lower bound of the 1-alpha CI for the PD
  \item \code{pd.u.ci} Upper bound of the 1-alpha CI for the PD
  \item \code{seed} Seed identity if bootstrapping is used
  \item \code{decision} NHST decision
  \item \code{alpha} Nominal Type I error rate
}
}
\description{
This function evaluates whether a certain predictor variable in a multiple regression model can be considered statistically and practically negligible according to a predefined interval. i.e., minimally meaningful effect size (MMES)/smallest effect size of interest (SESOI). Where the effect tested is the relationship between the predictor of interest and the outcome variable, holding all other predictors constant.
}
\details{
This function evaluates whether a certain predictor variable in a multiple regression model can be considered statistically and practically negligible according to a predefined interval. i.e., minimally meaningful effect size (MMES)/smallest effect size of interest (SESOI). Where the effect tested is the relationship between the predictor of interest and the outcome variable, holding all other predictors constant.

Unlike the most common null hypothesis significance tests looking to detect a difference or the existence of an effect statistically different than zero, in negligible effect testing, the hypotheses are flipped: In essence, H0 states that the effect is non-negligible, whereas H1 states that the effect is in fact statistically and practically negligible.

The statistical tests are based on Anderson-Hauck (1983) and Schuirmann's (1987) Two One-Sided Test (TOST) equivalence testing procedures; namely addressing the question of whether the estimated effect size (and its associated uncertainty) of a predictor variable in a multiple regression model is smaller than the what the user defines as negligible effect size. Defining what is considered negligible effect is done by specifying the negligible (equivalence) interval: its upper (eiu) and lower (eil) bounds.

The negligible (equivalence) interval should be set based on the context of the research. Because the predictor's effect size can be in either standardized or unstandardized units, setting eil and eiu is a matter of determining what magnitude of the relationship between predictor and outcome in either standardized or unstandardized units is the minimally meaningful effect size (MMES) given the context of the research.

It is necessary to be consistent with the units of measurement. For example, unstandardized negligible interval bounds (i.e., eil and eiu) must only be used when std = FALSE (default). If the effect size (b), standard error (se), and sample size (n) are entered manually as arguments (i.e., without the dataset), these should also be in the same units of measurements. Whereas if the user prefers to specify eiu and eil in standardized unites, std = TRUE should be specified. In which case, any units entered into the function must also be in standardized form. Mixing unstandardized and standardized units would yield inaccurate results and likely lead to invalid conclusions. Thus, users must be cognizant of the measurement units of the negligible interval.

There are two main approaches to using neg.reg. The first (and more recommended) is by inserting a dataset ('data' argument) into the function. If the user/s have access to the dataset, they should use the following set of arguments: data, formula, predictor, bootstrap (optional), nboot (optional), and seed (optional). However, this function also accommodates cases where no dataset is available. In this case, users should use the following set of arguments instead: b, se, n, and nop. In either situation, users must specify the negligible interval bounds (eiu and eil). Other optional arguments and features include: alpha, std, test, plots, and saveplots.

The proportional distance (PD; effect size/eiu) estimates the proportional distance of the estimated effect to eiu, and acts as an alternative effect size measure.

The confidence interval for the PD is computed via bootstrapping (percentile bootstrap), unless the user does not insert a dataset. In which case, the PD confidence interval is calculated by dividing the upper and lower CI bounds for the effect size estimate by the absolute value of the negligible interval bounds.
}
\examples{
# Negligible Regression Coefficient (equivalence interval: -.1 to .1)
pr1 <- stats::rnorm(20, mean = 0, sd= 1)
pr2 <- stats::rnorm(20, mean = 0, sd= 1)
dp <- stats::rnorm(20, mean = 0, sd= 1)
dat <- data.frame(pr1,pr2,dp)

# dataset available (unstandardized coefficients, AH procedure, using bootstrap-generated CIs):
neg.reg(formula=dp~pr1+pr2,data=dat,predictor=pr1,eil=-.1,eiu=.1,nboot=5)
neg.reg(formula=dat$dp ~ dat$pr1 + dat$pr2, predictor= pr1, eil= -.25, eiu= .25, nboot=5)

# dataset unavailable (standardized coefficients, TOST procedure):
neg.reg(b= .017, se= .025, nop= 3, n= 500, eil=-.1,eiu=.1, std=TRUE, test="TOST")
# end.



}
\author{
Udi Alter \email{udialter@gmail.com} and
  Alyssa Counsell \email{a.counsell@torontomu.ca} and
  Rob Cribbie \email{cribbie@yorku.ca}
}
