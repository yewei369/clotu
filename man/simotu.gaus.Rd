% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simotu.gaus.R
\name{simotu.gaus}
\alias{simotu.gaus}
\title{Simulation data with Gaussian distribution}
\usage{
simotu.gaus(
  n,
  p,
  t,
  nref = 3,
  unif.min = 0,
  unif.max = 0.35,
  full.mean = 20000,
  full.sd = 500,
  seed = F
)
}
\arguments{
\item{n, }{sample size for each target}

\item{p, }{number of OTUs}

\item{t, }{number of targets}

\item{nref, }{number of reference OTUs}

\item{unif.min, }{min limit of Uniform dist for non-zero probability}

\item{unif.max, }{max limit of Uniform dist for non-zero probability}

\item{full.mean, }{full range of mean value of total counts}

\item{full.sd, }{full range of standard deviation of total counts}

\item{seed, }{index of the seed, if F (default) no seed is set}
}
\value{
a list of OTU table and meta information
}
\description{
Generate simulation data with Gaussian distribution
}
\examples{
simotu.gaus(50,700,3,nref=5,full.mean=10000,unif.min=0,unif.max=0.4,seed=1234)
}
