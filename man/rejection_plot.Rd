\name{rejection_plot}

\Rdversion{1.1}

\alias{rejection_plot}

\title{
  Plot rejections vs. p-value cutoff
}

\description{
  Plot the number, or fraction, of null hypotheses rejected as a
  function of the p-value cutoff. Multiple sets of p-values are
  accepted, in a list or in the columns of a matrix, in order to permit
  comparisons. 
}

\usage{
rejection_plot(p,
               col, lty = 1, lwd = 1,   
               xlab = "p cutoff", ylab = "number of rejections",
               xlim = c(0, 1), ylim,
               legend = names(p),
               at = c("all", "sample"),
               n_at = 100,
               probability = FALSE,
               ...
               )
}

\arguments{
  
  \item{p}{
    The p-values to be used for plotting. These may be in the columns of
    a matrix, or in the elements of a list. One curve will be generated
    for each column/element, and all \code{NA} entries will be
    dropped. If column or element names are supplied, they are used by
    default for a plot legend.
  }

  \item{col}{
    Colors to be used for each curve plotted. Recycled if necessary. If
    \code{col} is omitted, \code{\link{rainbow}} is used to generate a
    set of colors.
  }
  
  \item{lty}{
    Line styles to be used for each curve plotted. Recycled if necessary.
  }
  
  \item{lwd}{
    Line widths to be used for each curve plotted. Recycled if necessary.
  }
  
  \item{xlab}{
    X-axis text label.
  }
  
  \item{ylab}{
    Y-axis text label.
  }
  
  \item{xlim}{
    X-axis limits.
  }
  
  \item{ylim}{
    Y-axis limits.
  }
  
  \item{legend}{
    Text for legend. Matrix column names or list element names (see
    \code{p} above) are used by default. If \code{NULL}, no legend is
    plotted.
  }
  
  \item{at}{
    Should step functions be plotted with a step at every value in
    \code{p}, or should linear interpolation be used at a sample of
    points spanning \code{xlim}? The latter looks when there are many
    p-values.
  }
  
  \item{n_at}{
    When \code{at = "sample"} is given, how many sample points should be
    used for interpolation and plotting?
  }
  
  \item{probability}{
    Should the fraction of null hypotheses rejected be reported instead
    of the count? See the \code{probability} argument to
    \code{\link{hist}}. 
  }
  
  \item{\dots}{
    Other arguments to pass to the \code{\link{plot}} call which sets up
    the axes. Note that the \code{...} argument will not be passed to
    the \code{\link{lines}} calls which actually generate the curves.
  }
  
}

\value{
  A list of the step functions used for plotting is returned invisibly.
}

\author{Richard Bourgon <bourgon@ebi.ac.uk>}

\examples{
# See the vignette: Diagnostic plots for independent filtering
}
