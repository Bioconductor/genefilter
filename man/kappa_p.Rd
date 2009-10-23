\name{kappa_p}

\Rdversion{1.1}

\alias{kappa_p}
\alias{kappa_t}

\title{
  Compute proportionality constant for fold change bound.
}

\description{
  Filtering on overall variance induces a lower bound on fold
  change. This bound depends on the significance of the evidence against
  the null hypothesis, an is a multiple of the cutoff used for an
  overall variance filter. It also depends on sample size in both of the
  groups being compared. These functions compute the multiplier for the
  supplied p-values or t-statistics.
}

\usage{
kappa_p(p, n1, n2 = n1)
kappa_t(t, n1, n2 = n1)
}

\arguments{
  \item{p}{The p-values at which to compute the multiplier.}
  \item{t}{The t-statistics at which to compute the multiplier.}
  \item{n1}{Sample size for class 1.}
  \item{n2}{Sample size for class 2.}
}

\value{
  A vector of multipliers: one per p-value or t-static in
  \code{p} or \code{t}.
}

\author{Richard Bourgon <bourgon@ebi.ac.uk>}
