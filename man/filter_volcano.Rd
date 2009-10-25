\name{filter_volcano}

\Rdversion{1.1}

\alias{filter_volcano}

\title{Volcano plot for overall variance filtering}

\description{
  Generate a volcano plot contrasting p-value with fold change (on the
  log scale), in order to visualize the effect of filtering on overall
  variance and also assign significance via p-value.
}

\usage{
filter_volcano(
               d, p, S,
               n1, n2,
               alpha, S_cutoff,
               cex = 0.5, pch = 19,
               xlab = expression(paste(log[2], " fold change")),
               ylab = expression(paste("-", log[10], " p")),
               cols = c("grey80", "grey50", "black"),
               ltys = c(1, 3),
               use_legend = TRUE,
               ...
               )
}

\arguments{
  \item{d}{Fold changes, typically on the log scale, base 2.}
  
  \item{p}{The p-values}

  \item{S}{
    The overall standard deviation filter statistics, i.e., the
    square roots of the overall variance filter statistics.
  }

  \item{n1}{Sample size for group 1.}

  \item{n2}{Sample size for group 2.}
  
  \item{alpha}{Significance cutoff used for p-values.}
  
  \item{S_cutoff}{
    Filter cutoff used for the overall standard deviation in \code{S}.
  }
  
  \item{cex}{Point size for plotting.}

  \item{pch}{Point character for plotting.}
  
  \item{xlab}{Label for x-axis.}
  
  \item{ylab}{Label for y-axis.}
  
  \item{cols}{
    A vector of three colors used for plotting. These correspond to
    filtered data, data which pass the filter but are insignificant, and
    data pass the filter and are also statistically significant.
  }
  
  \item{ltys}{
    The induced bound on log-scale fold change is plotted, as is the
    significance cutoff for data passing the filter. The \code{ltys}
    argument gives line styles for these drawing these two thresholds on
    the plot.
  }

  \item{use_legend}{Should a legend for point color be produced?}

  \item{\dots}{Other arguments for \code{plot}.}
  
}

\author{Richard Bourgon <bourgon@ebi.ac.uk>}

\examples{
# See the vignette: Diagnostic plots for independent filtering
}
