\name{isoph}
\alias{isoph}
\title{Fit Isotonic Proportional Hazards Model}
\description{Nonparametric estimation of an isotonic covariate effect for Cox's partial likelihood.}
\usage{
  isoph(formula, data, shape, K, maxiter, eps, maxdec)
}
\arguments{
  \item{formula}{a formula object: a response ~ a univariate covariate. The response must be survival outcome unsing the Surv function in the survival package. For interval data, each subject's starting time must be set to 0.}
  \item{data}{data.frame or list that includes variables named in the formula argument.}
  \item{shape}{direction of the covariate effect on the hazard function, "increasing" or "decreasing"}
  \item{K}{anchor constraint (default is 0).}
  \item{maxiter}{maximum number of iteration (default is 10^5).}
  \item{eps}{stopping convergence criteria (default is 10^-3).}
  \item{maxdec}{maximum number of decisimal for output (default is 2).}
}
\details{The isoph function allows to analyze isotonic proportional hazards model, defined as
\deqn{ \lambda(t|z)=\lambda0(t)exp(\psi(Z)), }
where \eqn{ \lambda0 } is a baseline hazard function and \eqn{ \psi } is an isotonic function. Correspondingly, the partial likelihood is defined, so that, \eqn{ \psi } can be estimated by maximizing the partial likelihood over an isotonic constraint.

As similar to set a reference group in the standard Cox's proportional hazard model, one point has to be fixed with \eqn{ \psi(K)=0 }, where \eqn{ K } is an anchor point. A direction of \eqn{ \psi } should be defined as monotone increasing or monotone decreasing in \eqn{Z}, prior to data analysis. Time-dependent covariate is allowed with interval data, where each subject's starting time must be set to zero.

Pseudo iterative convex minorant algorithm is used.
}

\references{Yunro Chung, Anastasia Ivanova, Michael M. Hudgens, Jason P. Fine, Partial likelihood estimation of isotonic proportional hazards models (under revision).}
\author{Yunro Chung [cre], Anastasia Ivanova, Michael G. Hudgens and Jason P. Fine}

\examples{
#require(survival)
#require(Iso)

###
# 1. time-independent covariate with monotone increasing effect
###
# 1.1. create a test data set 1
test1=list(
  time=  c(2, 5, 1, 7, 9, 5, 3, 6, 8, 9, 7, 4, 5, 2, 8),
  status=c(0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1),
  z=     c(2, 1, 1, 3, 5, 6, 7, 9, 3, 0, 2, 7, 3, 9, 4)
)

# 1.2. Fit isotonic proportional hazards model
res1 = isoph(Surv(time, status)~z, data=test1, shape="increasing")

# 1.3. print result
res1

# 1.4 Figure
plot(res1, which=1) #which=1 for psi.hat (default)
plot(res1, which=2) #which=2 for hazard ratio, i.e. exp(psi.hat)

###
# 2. time-dependent covariate with monotone decreasing effect
###
# 2.1. create a test data set 2
test2=list(
  start= c(0, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 2),
  stop=  c(5, 2, 6, 3, 8, 1, 2, 5, 7, 2 ,2, 8),
  status=c(1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1),
  z=     c(1, 2, 1, 3, 5, 3, 1, 3, 2, 8, 5, 7)
)

# 2.2. fit isotonic proportional hazards model
res2 = isoph(Surv(start,stop,status)~z, data=test2, shape="decreasing")

# 2.3. print result
res2

# 2.4. Figure
plot(res2, which=1) #which=1 for psi.hat (default)
plot(res2, which=2) #which=2 for hazard ratio, i.e. exp(psi.hat)

###
#3. more arguments for plot.isoph
###
#3.1. renames labels
plot(res2, main="Iso PH", ylab="Iso", xlab="Cov", lglab="Cov wt obs", lgloc="center", lgcex=1.5)

#3.2. removes labels and changes line and point parameters
plot(res2, main=NA, ylab=NA, xlab=NA, lglab=NA, lty=2, lcol=2, lwd=2, pch=3, pcol=4, pcex=1.5)
}

\keyword{Isotonic regression, Survival analysis, Constrained estimation}

