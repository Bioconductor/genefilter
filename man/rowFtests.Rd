\name{rowFtests}
\alias{colFtests}
\alias{rowFtests}
\alias{rowFtests,matrix,factor-method}
\alias{rowFtests,exprSet,factor-method}
\alias{rowFtests,exprSet,character-method}
\alias{rowFtests,ExpressionSet,factor-method}
\alias{rowFtests,ExpressionSet,character-method}

\alias{rowttests}
\alias{rowttests,matrix,factor-method}
\alias{rowttests,matrix,missing-method}
\alias{rowttests,exprSet,factor-method}
\alias{rowttests,exprSet,character-method}
\alias{rowttests,exprSet,missing-method}
\alias{rowttests,ExpressionSet,factor-method}
\alias{rowttests,ExpressionSet,character-method}
\alias{rowttests,ExpressionSet,missing-method}

\alias{colttests}
\alias{fastT}

\title{t-tests and F-tests for rows or columns of a matrix}
\description{t-tests and F-tests for rows or columns of a matrix}

\usage{
rowttests(x, fac, tstatOnly = FALSE) 
colttests(x, fac, tstatOnly = FALSE)
fastT(x, ig1, ig2, var.equal = TRUE)

rowFtests(x, fac, var.equal = TRUE)
colFtests(x, fac, var.equal = TRUE)
}
\arguments{
  \item{x}{\code{exprSet} or numeric matrix. The matrix must not
    contain \code{NA} values.}
  \item{fac}{Factor; if \code{x} is an \code{exprSet}, this may also be
    a character vector of length 1 with the name of a covariate
    variable in \code{x}.

    For the \code{row*} functions, the length of the factor must be
    the same as the number of columns of \code{x}.
    For the \code{col*} functions, it must be the same as the number
    of rows of \code{x}.
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
  \code{fastT} is a wrapper for \code{rowttests}, useful
  for some legacy code.

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

   for (j in 1:nrow(x)) {
     s1 = t.test(x[j,])
     about.equal(s1$statistic, r1$statistic[j])
     about.equal(s1$p.value,   r1$p.value[j])

     s2 = t.test(x[j,] ~ f2, var.equal=TRUE)
     about.equal(s2$statistic, r2$statistic[j])
     about.equal(s2$p.value,   r2$p.value[j])

     dm = -diff(tapply(x[j,], f2, mean))
     about.equal(dm, r2$dm[j])

     s7 = summary(lm(x[j,]~f7))
     about.equal(s7$statistic$value, r7$statistic[j])
   }

   ## colttests
   c2 = colttests(t(x), f2)
   stopifnot(identical(r2, c2))
 
   ## missing values
   f2n = f2
   f2n[sample(length(f2n), 3)] = NA
   r2n = rowttests(x, f2n)
   for(j in 1:nrow(x)) {
     s2n = t.test(x[j,] ~ f2n, var.equal=TRUE)
     about.equal(s2n$statistic, r2n$statistic[j])
     about.equal(s2n$p.value,   r2n$p.value[j])
  }

}


}
\keyword{math}
