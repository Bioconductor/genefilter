\name{eSetFilter}
\alias{eSetFilter}
\alias{getFilterNames}
\alias{getFuncDesc}
\alias{getRdAsText}
\alias{parseDesc}
\alias{parseArgs}
\alias{setESetArgs}
\alias{isESet}
\alias{showESet}
\title{A function to filter an eSet object}
\description{
  Given a Bioconductor's ExpressionSet object, this function filters
  genes using a set of selected filters.
}
\usage{
eSetFilter(eSet)
getFilterNames()
getFuncDesc(lib = "genefilter", funcs = getFilterNames())
getRdAsText(lib)
parseDesc(text)
parseArgs(text)
showESet(eSet)
setESetArgs(filter)
isESet(eSet)
}
\arguments{
  \item{eSet}{\code{eSet} an ExpressionSet object}
  \item{lib}{\code{lib} a character string for the name of an R library
    where functions of interests reside}
  \item{funcs}{\code{funcs} a vector of character strings for names of
    functions of interest}
  \item{text}{\code{text} a character of string from a filed
    (e. g. description, argument, ..) filed of an Rd file for a
    fucntion}
  \item{filter}{\code{filter} a character string for the name of a
  filter function} 
}
\details{
  A set of filters may be selected to filter genes in through each of
  the filters in the order the filters have been selected
}
\value{
  A logical vector of length equal to the number of rows of 'expr'.
  The values in that vector indicate whether the corresponding row
  of 'expr' passed the set of filter functions.
}
\author{Jianhua Zhang}

\seealso{\code{\link{genefilter}}}
\examples{
 if( interactive() ) {
   data(sample.ExpressionSet)      
   res <- eSetFilter(sample.ExpressionSet)
 }
}
\keyword{manip}

