\name{genefilter}
\alias{genefilter}
\title{ A function to filter genes in an expressin array. }
\description{
 Genefilter takes a dataframe or matrix of expression data and a list of 
filtering functions and sequentially applies the filtering functions to 
the array. The reduced expression array (containing genes that passed all
filters) is returned.
}
\usage{
genefilter(expr, flist)
}

\arguments{
  \item{expr}{ A matrix or dataframe of expressions. }
  \item{flist}{ A list of filtering functions. }
}
\details{
As above.
}
\value{
 The subset of \code{expr} that passes all the filters.
}

\author{Robert Gentleman }

\seealso{  \code{\link{allNA}} }

\examples{
   #generate some data
   set.seed(3333)
   ngenes <- 1000
   m1 <- rnorm(ngenes,4,1)
   m2 <- rnorm(ngenes,4,1)
   Exp1 <- cbind(matrix(exp(rnorm(ngenes*5,m1,1)),nrow=ngenes),
                    matrix(exp(rnorm(ngenes*10,m2,1)),nrow=ngenes))
   Exp1[Exp1<20] <- 20
   sub <- Exp1>20 & Exp1<150
   Exp1[sub] <- ifelse(runif(length(sub[sub]))<.5,20,Exp1[sub])
   dimnames(Exp1) <- list(paste('x',format(1:ngenes,justify='l'),sep=''),
         paste('sample',format(1:ncol(Exp1),justify='l'),sep=''))
        dimnames(Exp1) <- list(paste('x',1:ngenes,sep=''),
                          paste('sample',1:ncol(Exp1),sep=''))
   Exp1 <- round(Exp1)
   s1.flist<-list(maxA, sdom, poverA)
   Expsub <- genefilter(Exp1, s1.flist)
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
