\name{Anova}
\alias{Anova}
\title{A filter function for Analysis of Variance }
\description{
  \code{Anova} returns a function of one argument with bindings for
  \code{cov} and \code{p}.
  The function, when evaluated, performs an ANOVA using \code{cov} as
  the covariate. It returns \code{TRUE} if the p value for a difference
  in means is less than \code{p}. 
}
\usage{
Anova(cov, p=0.05, na.rm=TRUE)
}
\arguments{
  \item{cov}{The covariate. It must have length equal to the number of
    columns of the array that \code{Anova} will be applied to. }
  \item{p}{ The p-value for the test. }
  \item{na.rm}{If set to \code{TRUE} any \code{NA}'s will be removed. }
}

\details{
  The function returned by \code{Anova} uses \code{lm} to fit a linear
  model of the form
  \code{lm(x ~ cov)}, where \code{x} is the set of gene expressions.
  The F statistic for an overall effect is computed and if it has a
  \emph{p}-value less than \code{p} the function returns \code{TRUE},
  otherwise it returns \code{FALSE} for that gene.
  }
\value{
 \code{Anova} returns a function with bindings for \code{cov} and
 \code{p} that will perform a one-way ANOVA.

  The covariate can be continuous, in which case the test is for a linear 
  effect for the covariate.
}
\author{R. Gentleman }
\seealso{\code{\link{kOverA}}, \code{\link{lm}} }

\examples{
  set.seed(123)
  af <- Anova(c(rep(1,5),rep(2,5)), .01)
  af(rnorm(10))
}
\keyword{manip}
