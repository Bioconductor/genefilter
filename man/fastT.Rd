\name{fastT}
\alias{fastT}
\title{ Compute t-tests quickly. }
\description{
  A function to compute two-sample t-tests quickly for each row of a
  data matrix. 
}
\usage{
fastT(data, ig1, ig2, var.equal=TRUE)
}
\arguments{
  \item{data}{A matrix of numbers. }
  \item{ig1}{ The indices of the columns of \code{data} that correspond
    to group 1.}
  \item{ig2}{The indices of the columns of \code{data} that correspond
    to group 2. }
  \item{var.equal}{a logical variable indicating whether to treat the two
    variances as being equal. Currently, this parameter is
    \emph{ignored}, and the pooled variance is always used.}
}
\details{
  Implemented in Fortran, by R. Gray. This is basically the same code as
  is found in the permax package. The implementation is currently in
  single precision.
}
\value{
A list with three components
  \item{z }{A vector of row-wise t-statistics.}
  \item{dm }{A vector of row-wise difference in mean.}
  \item{var.equal}{The value used for var.equal.}
}
\author{ R. Gray }

\seealso{ \code{\link{permax}} }

\examples{
 x = matrix(rnorm(1000), nc=10)
 gp1 = 1:4
 gp2 = 5:10
 ans = fastT(x, gp1, gp2)

 for(i in 1:nrow(x))
   stopifnot(abs(ans$z[i] -  t.test(x[i, gp1], x[i, gp2],
       var.equal=TRUE)$statistic) < 1e-6)
}
\keyword{htest }
