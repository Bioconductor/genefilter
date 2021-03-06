\name{rowFtests}
\alias{rowFtests}
\alias{rowFtests,matrix,factor-method}
\alias{rowFtests,ExpressionSet,factor-method}
\alias{rowFtests,ExpressionSet,character-method}

\alias{colFtests}
\alias{colFtests,matrix,factor-method}
\alias{colFtests,ExpressionSet,factor-method}
\alias{colFtests,ExpressionSet,character-method}

\alias{rowttests}
\alias{rowttests,matrix,factor-method}
\alias{rowttests,matrix,missing-method}
\alias{rowttests,ExpressionSet,factor-method}
\alias{rowttests,ExpressionSet,character-method}
\alias{rowttests,ExpressionSet,missing-method}

\alias{colttests}
\alias{colttests,matrix,factor-method}
\alias{colttests,matrix,missing-method}
\alias{colttests,ExpressionSet,factor-method}
\alias{colttests,ExpressionSet,character-method}
\alias{colttests,ExpressionSet,missing-method}

\alias{fastT}

\title{t-tests and F-tests for rows or columns of a matrix}
\description{t-tests and F-tests for rows or columns of a
  matrix, intended to be speed efficient.}
\usage{
rowttests(x, fac, tstatOnly = FALSE, na.rm = FALSE)
colttests(x, fac, tstatOnly = FALSE, na.rm = FALSE)
fastT(x, ig1, ig2, var.equal = TRUE)

rowFtests(x, fac, var.equal = TRUE)
colFtests(x, fac, var.equal = TRUE)
}
\arguments{
  \item{x}{Numeric matrix. The matrix must not contain \code{NA} values.
    For \code{rowttests} and \code{colttests}, \code{x} can also be an
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}.}

  \item{fac}{Factor which codes the grouping to be tested.
    There must be 1 or 2 groups for the t-tests (corresponding to one-
    and two-sample t-test), and 2 or more for the F-tests. If \code{fac}
    is missing, this is taken as a one-group test (i.e. is only allowed
    for the t-tests). The length
    of the factor needs to correspond to the sample size:
    for the \code{row*} functions, the length of the factor must be
    the same as the number of columns of \code{x},
    for the \code{col*} functions, it must be the same as the number
    of rows of \code{x}.

    If \code{x} is an
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}, then 
    \code{fac} may also be a character vector of length 1 with the
    name of a covariate in \code{x}.}

  \item{tstatOnly}{A logical variable indicating whether to calculate
    p-values from the t-distribution with appropriate degrees of
    freedom. If \code{TRUE}, just the t-statistics are
    returned. This can be considerably faster.}

  \item{na.rm}{A logical variable indicating whether to remove NA values
    prior to calculation test statistics.}

  \item{ig1}{The indices of the columns of \code{x} that correspond
    to group 1.}
  \item{ig2}{The indices of the columns of \code{x} that correspond
    to group 2.}

  \item{var.equal}{A logical variable indicating whether to treat the
    variances in the samples as equal.  If 'TRUE', a simple F test for
    the equality of means in a one-way analysis of variance is
    performed.  If 'FALSE', an approximate method of Welch (1951) is
    used, which generalizes the commonly known 2-sample Welch test to
    the case of arbitrarily many samples.}
}

\details{
  If \code{fac} is specified, \code{rowttests} performs for each
  row of \code{x} a two-sided, two-class t-test with equal variances.
  \code{fac} must be a factor of length \code{ncol(x)} with two levels,
  corresponding to the two groups. The sign of the resulting t-statistic
  corresponds to "group 1 minus group 2".
  If \code{fac} is missing, \code{rowttests} performs for each row of
  \code{x} a two-sided one-class t-test against the null hypothesis 'mean=0'.
  
  \code{rowttests} and \code{colttests} are implemented in C and
  should be reasonably fast and memory-efficient.
  \code{fastT} is an alternative implementation, in Fortran, possibly useful
  for certain legacy code.
  \code{rowFtests} and \code{colFtests} are currently implemented using
  matrix algebra in R. Compared to the \code{rowttests} and
  \code{colttests} functions,
  they are slower and use more memory.
}

\value{
  A \code{data.frame} with columns \code{statistic},
  \code{p.value} (optional in the case of the t-test functions) and
  \code{dm}, the difference of the group means (only in the
  case of the t-test functions).
  The \code{row.names} of the data.frame are taken from the
  corresponding dimension names of \code{x}.
  
  The degrees of freedom are provided in the attribute \code{df}.
  For the F-tests, if \code{var.equal} is 'FALSE', \code{nrow(x)+1}
  degree of freedoms 
  are given, the first one is the first degree of freedom (it is the
  same for each row) and the other ones are the second degree of freedom
  (one for each row). 
}
\references{B. L. Welch (1951), On the comparison of several mean values: an
     alternative approach. Biometrika, *38*, 330-336}
\author{Wolfgang Huber <whuber@embl.de>}
\seealso{\code{\link[multtest:mt.teststat]{mt.teststat}}}
\examples{
   ##
   ## example data
   ##
   x  = matrix(runif(40), nrow=4, ncol=10)
   f2 = factor(floor(runif(ncol(x))*2))
   f4 = factor(floor(runif(ncol(x))*4))

   ##
   ## one- and two group row t-test; 4-group F-test
   ##
   r1 = rowttests(x)
   r2 = rowttests(x, f2)
   r4 = rowFtests(x, f4)

   ## approximate equality
   about.equal = function(x,y,tol=1e-10)
     stopifnot(is.numeric(x), is.numeric(y), length(x)==length(y), all(abs(x-y) < tol))

   ##
   ## compare with the implementation in t.test
   ##
   for (j in 1:nrow(x)) {
     s1 = t.test(x[j,])
     about.equal(s1$statistic, r1$statistic[j])
     about.equal(s1$p.value,   r1$p.value[j])

     s2 = t.test(x[j,] ~ f2, var.equal=TRUE)
     about.equal(s2$statistic, r2$statistic[j])
     about.equal(s2$p.value,   r2$p.value[j])

     dm = -diff(tapply(x[j,], f2, mean))
     about.equal(dm, r2$dm[j])

     s4 = summary(lm(x[j,] ~ f4))
     about.equal(s4$fstatistic["value"], r4$statistic[j])
   }

   ##
   ## colttests
   ##
   c2 = colttests(t(x), f2)
   stopifnot(identical(r2, c2))

   ##
   ## missing values
   ##
   f2n = f2
   f2n[sample(length(f2n), 3)] = NA
   r2n = rowttests(x, f2n)
   for(j in 1:nrow(x)) {
     s2n = t.test(x[j,] ~ f2n, var.equal=TRUE)
     about.equal(s2n$statistic, r2n$statistic[j])
     about.equal(s2n$p.value,   r2n$p.value[j])
   }

   ##
   ## larger sample size
   ##
   x  = matrix(runif(1000000), nrow=4, ncol=250000)
   f2 = factor(floor(runif(ncol(x))*2))
   r2 = rowttests(x, f2) 
   for (j in 1:nrow(x)) {
     s2 = t.test(x[j,] ~ f2, var.equal=TRUE)
     about.equal(s2$statistic, r2$statistic[j])
     about.equal(s2$p.value,   r2$p.value[j])
   }

   ## single row matrix
   rowFtests(matrix(runif(10),1,10),as.factor(c(rep(1,5),rep(2,5))))
   rowttests(matrix(runif(10),1,10),as.factor(c(rep(1,5),rep(2,5))))
}
\keyword{math}
