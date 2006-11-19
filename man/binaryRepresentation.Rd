\name{binaryRepresentation}
\alias{binaryRepresentation}
\title{Binary representation of integer vectors}
\description{Binary representation of integer vectors by a logical matrix}
\usage{binaryRepresentation(n, ndigits)}
\arguments{
  \item{n}{Integer vector, only positive values are expected.}
  \item{ndigits}{Integer scalar. How many digits should the output
    use. If missing, the mimimum possible value is chosen necessary to represent
    all numbers in \code{n}.}
}

\details{The implementation of this function is naive and slow (it
  really should be done in C), but I have found
  it useful in several different places and so decided to put it in this
  package. 
}

\value{Logical matrix whose columns correspond to the elements of
  \code{n} and whose rows to the digits.}
\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\examples{
  binaryRepresentation(0:15)
  binaryRepresentation(64:69, ndigits=8)
}
\keyword{arith}

