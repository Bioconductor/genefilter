\name{coxfilter}
\alias{coxfilter}

\title{A filter function for univariate Cox regression. }
\description{
  A function that performs Cox regression with bindings for \code{surt},
  \code{cens}, and \code{p} is returned. This function filters genes
  according to the attained p-value from a Cox regression using
  \code{surt} as the survival times, and \code{cens} as the censoring
  indicator. It requires \code{survival}.
}
\usage{
coxfilter(surt, cens, p)
}
\arguments{
  \item{surt}{Survival times.}
  \item{cens}{Censoring indicator. }
  \item{p}{The p-value to use in filtering. }
}

\value{
  Calls to the \code{\link{coxph}} function in the \code{survival}
    library are used to fit a Cox model. The filter function returns
    \code{TRUE} if the p-value in the fit is less than \code{p}.
}

\author{R. Gentleman }
\seealso{\code{\link{Anova}}}

\examples{
   set.seed(-5)
   sfun<- coxfilter(rexp(10), ifelse(runif(10) < .7, 1, 0), .05)
   dat <- matrix(rnorm(1000), ncol=10)
   out <- genefilter(dat, list(sfun))
}
\keyword{manip}

