\name{rowFtests}
\alias{colFtests}
\alias{rowFtests}
\alias{rowttests}
\alias{colttests}
\alias{fastTNew}

\title{t-tests and F-tests for rows or columns of a matrix}
\description{t-tests and F-tests for rows or columns of a matrix}

\usage{
rowttests(x, fac, tstatOnly = FALSE) 
colttests(x, fac, tstatOnly = FALSE)
fastTNew(x, ig1, ig2)

rowFtests(x, fac, var.equal = TRUE)
colFtests(x, fac, var.equal = TRUE)
}
\arguments{
  \item{x}{Numeric matrix. Must not contain \code{NA} values.}
  \item{fac}{Factor. For the \code{row*} functions, its length must be
    the same as the number of columns of \code{x}.
    For the \code{col*} functions, as the number of rows of \code{x}.
    For the \code{*ttests} function, \code{fac} must have exactly
    two levels.}
  \item{tstatOnly}{a logical variable indicating whether to calculate
    parametric p-values.  If \code{FALSE}, just the t-statistics are
    returned. This can be considerably faster.}
  \item{ig1}{The indices of the columns of \code{x} that correspond
    to group 1.}
  \item{ig2}{The indices of the columns of \code{x} that correspond
    to group 2. }
  \item{var.equal}{a logical variable indicating whether to treat the
    variances in the samples as equal.  If 'TRUE', a simple F test for
    the equality of means in a one-way analysis of variance is
    preformed.  If 'FALSE', an approximate method of Welch (1951) is
    used, which generalizes the commonly known 2-sample Welch test to
    the case of arbitrarily many samples.}
}
\details{
  \code{rowttests} and \code{colttests} are implemented in C and
  are reasonably fast and memory-efficient.
  \code{fastTNew} is a wrapper for \code{rowttests}. It is intended
  to replace \code{fastT}, which we need around to deal with legacy
  code.

  If \code{fac} is specified, \code{rowttests} performs for each
  row of \code{x} a two-sided, two-class t-test with equal variances.
  \code{fac} must be a factor of length \code{ncol(x)} with two levels,
  corresponding to the two groups. The sign of the resulting t-statistic
  corresponds to "group 1 minus group 2".
  
  If \code{fac} is missing, \code{rowttests} performs for each row of
  \code{x} a two-sided one-class t-test against the null hypothesis 'mean=0'.

  \code{rowFtests} and \code{colFtests} are currently implemented using
  matrix algebra in R. Compared to the \code{*ttests} functions,
  they are slower and use more memory.

  If \code{var.equal} is 'FALSE', \code{nrow(x)+1} degree of freedoms
  are given, the first one is the first degree of freedom (it is the
  same for each row) and the other ones are the second degree of freedom
  (one for each row). 
}
\value{
  For the \code{row*} and \code{col*} functions, 
  a list with the test statistics, 
  p-values, and degrees of freedom.
  Additonally, the t-test functions return the differences
  of group means.
}
\references{B. L. Welch (1951), On the comparison of several mean values: an
     alternative approach. Biometrika, *38*, 330-336}
\author{Wolfgang Huber <huber@ebi.ac.uk>}
\seealso{\code{\link[multtest:mt.teststat]{mt.teststat}}}
\examples{
   x  = matrix(runif(970), ncol=97)
   f2 = factor(floor(runif(ncol(x))*2))
   f7 = factor(floor(runif(ncol(x))*7))

   r1 = rowttests(x)
   r2 = rowttests(x, f2)
   r7 = rowFtests(x, f7)

   ## compare with pedestrian tests
   about.equal = function(x,y,tol=1e-10)
     stopifnot(all(abs(x-y) < tol))

   s1 = t.test(x[1,])
   about.equal(s1$statistic, r1$statistic[1])
   about.equal(s1$p.value,   r1$p.value[1])

   s2 = t.test(x[1,] ~ f2, var.equal=TRUE)
   about.equal(s2$statistic, r2$statistic[1])
   about.equal(s2$p.value,   r2$p.value[1])

   dm = -diff(tapply(x[1,], f2, mean))
   about.equal(dm, r2$dm[1])

   s7 = summary(lm(x[1,]~f7))
   about.equal(s7$statistic$value, r7$statistic[1])

   ## colttests
   c2 = colttests(t(x), f2)
   stopifnot(identical(r2, c2))

   ## this can be removed once the old fastT is gone
   x   = matrix(runif(1e7), nrow=1e4)
   ig1 = sample(ncol(x), ncol(x)/2)
   ig2 = (1:ncol(x))[-ig1]

   print(system.time(t1<-fastTNew(x, ig1, ig2)))
   print(system.time(t2<-fastT(x, ig1, ig2, ratio=FALSE)))

   about.equal(t1$statistic, t2$z,  1e-6)
   about.equal(t1$dm,        t2$dm, 1e-6)

}


}
\keyword{math}
