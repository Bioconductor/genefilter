\name{filtered_p}

\Rdversion{1.1}

\alias{filtered_p}
\alias{filtered_R}

\title{
  Compute and adjust p-values, with filtering
}

\description{
  Given filter and test statistics in the form of unadjusted p-values,
  or functions able to compute these statistics from the data, filter
  and then correct the p-values across a range of filtering
  stringencies. 
}

\usage{
filtered_p(filter, test, theta, data, method = "none")
filtered_R(alpha, filter, test, theta, data, method = "none")
}

\arguments{

  \item{alpha}{
    A cutoff to which p-values, possibly adjusted for multiple testing,
    will be compared.
  }
  
  \item{filter}{
    A vector of stage-one filter statistics, or a function which is able
    to compute this vector from \code{data}, if \code{data} is supplied.
  }
  
  \item{test}{
    A vector of unadjusted p-values, or a function which is able
    to compute this vector from the filtered portion of \code{data}, if
    \code{data} is supplied. The option to supply a function is useful
    when the value of the test statistic depends on which hypotheses are
    filtered out at stage one. (The \pkg{limma} t-statistic is an
    example.) 
  }
  
  \item{theta}{
    A vector with one or more filtering fractions to consider. Actual
    cutoffs are then computed internally by applying
    \code{\link{quantile}} to the filter statistics contained in (or
    produced by) the \code{filter} argument.
  }
  
  \item{data}{
    If \code{filter} and/or \code{test} are functions rather than
    vectors of statistics, they will be applied to \code{data}. The
    functions will be passed the whole \code{data} object, and must work
    over rows, etc. themselves as appropriate.
  }
  
  \item{method}{
    The unadjusted p-values contained in (or produced by) \code{test}
    will be adjusted for multiple testing after filtering, using the
    \code{\link{p.adjust}} function in the \pkg{stats} package. See the
    \code{method} argument there for options.
  }p
  
}

\value{
  For \code{filtered_p}, a matrix of p-values, possible adjusted for
  multiple testing, with one row per null hypothesis and one column per
  filtering fraction given in \code{theta}. For a given column, entries
  which have been filtered out are \code{NA}.

  For \code{filtered_R}, a count of the entries in the \code{filtered_p}
  result which are less than \code{alpha}.
}

\author{Richard Bourgon <bourgon@ebi.ac.uk>}

\examples{
# See the vignette: Diagnostic plots for independent filtering
}

\seealso{
  See \code{\link{rejection_plot}} for visualization of
  \code{filtered_p} results.
}
