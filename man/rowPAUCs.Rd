\name{rowPAUCs}
\alias{rowPAUCs}

\title{Fixme.}
\description{Fixme.}

\usage{
rowPAUCs(x, fac, cutpts, p=0.1)
}
\arguments{
  \item{x}{\code{exprSet} or numeric matrix. The matrix must not
    contain \code{NA} values.}
  \item{fac}{Factor; if \code{x} is an \code{exprSet}, this may also be
    a character vector of length 1 with the name of a covariate
    variable in \code{x}.}
  \item{cutpts}{Fixme..}
  \item{p}{Fixme..}
}
\details{
  Fixme.
}
\value{
  A list with ...
}
\references{TO DO: Pepe}
\author{Florian Hahne <f.hahne@dkfz.de>}
\seealso{\code{\link[ROC:Fixme.]{Fixme.}}}
\examples{
  
   x  = matrix(runif(970), ncol=97)
   f2 = floor(runif(ncol(x))*2)

   r1 = rowttests(x, f2)
   r2 = rowPAUCs(x, f2)

 if(interactive())
   plot(r1$statistic, r2$BLA)
}


}
\keyword{math}
