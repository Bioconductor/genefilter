\name{rowFtests}
\alias{rowFtests}
\alias{rowttests}
\title{F-test and t-test for rows of a matrix}
\description{F-test and t-test for rows of a matrix}

\usage{
rowFtests(x, fac,var.equal = TRUE)
rowttests(x, fac)} 

\arguments{
  \item{x}{Matrix}
  \item{fac}{Factor, with \code{length(fac)=ncol(x)}. For
    \code{rowttests}, \code{fac} must have exactly two levels.}
  \item{var.equal}{a logical variable indicating whether to treat the
    variances in the samples as equal.  If 'TRUE', a simple F test for
    the equality of means in a one-way analysis of variance is
    preformed.  If 'FALSE', an approximate method of Welch (1951) is
    used, which generalizes the commonly known 2-sample Welch test to the case of arbitrarily many samples.}
}
\details{If \code{fac} is specified, \code{rowttests} performs for each
  row of \code{x} a two-sided, two-class t-test with equal variances.
  \code{fac} must be a factor of length \code{ncol(x)} with two levels,
  corresponding to the two groups.
  
  If \code{fac} is missing, \code{rowttests} performs for each row of
  \code{x} a two-sided one-class t-test against the null hypothesis 'mean=0'.

  With the current implementation, \code{NA}s are not allowed. This
  could be changed if desirable.

  If \code{var.equal} is 'FALSE' \code{nrow(x)} +1 degree of freedoms
  are given, the first one is the first degree of freedom (it is the
  same for each row) and the other ones are the second degree of freedom (one for
  each row). 
}
\value{  
A list with the test statistics, p-values, degrees of freedom.
}
\references{B. L. Welch (1951), On the comparison of several mean values: an
     alternative approach. Biometrika, *38*, 330-336}
\author{Wolfgang Huber <w.huber@dkfz.de>}
\seealso{\code{\link[multtest:mt.teststat]{mt.teststat}}}
\examples{
   x  = matrix(runif(9700), ncol=97)
   k2 = floor(runif(ncol(x))*2)
   k7 = floor(runif(ncol(x))*7)

   t1 = rowttests(x)
   t2 = rowttests(x, factor(k2))
   ft = rowFtests(x, factor(k7))

   if(require(multtest)) {
     fs = mt.teststat(x, k7, test="f")
     stopifnot(all(abs(fs - ft$statistic) < 1e-6))

     ts2 = mt.teststat(x, k2, test="t.equalvar")
     stopifnot(all(abs(ts2 - t2$statistic) < 1e-6))

     ts1 = apply(x, 1, function(x) t.test(x)$statistic)
     stopifnot(all(abs(ts1 - t1$statistic) < 1e-6))
   }
}
\keyword{math}
