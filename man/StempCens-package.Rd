\name{StempCens-package}
\alias{StempCens-package}
%\alias{StempCens}
\docType{package}
\title{
\packageTitle{StempCens}
}
\description{
\packageDescription{StempCens}
}
\details{
%The DESCRIPTION file:
%\packageDESCRIPTION{StempCens}
%\packageIndices{StempCens}

The functions provided are:

- \code{\link{CovarianceM}}: computes the spatio-temporal covariance matrix for balanced data.

- \code{\link{EffectiveRange}}: computes the effective range for an isotropic spatial correlation function.

- \code{\link{EstStempCens}}: returns the maximum likelihood estimates of the unknown parameters.

- \code{\link{PredStempCens}}: performs spatio-temporal prediction in a set of new \code{S} spatial locations for fixed time points.

- \code{\link{CrossStempCens}}: performs cross-validation, which measure the performance of the predictive model on new test dataset.

- \code{\link{DiagStempCens}}: returns measures and graphics for diagnostic analysis.
}
\author{

Larissa A. Matos (\href{https://orcid.org/0000-0002-2635-0901}{ORCID}), Katherine L. Valeriano (\href{https://orcid.org/0000-0001-6388-4753}{ORCID}) and Victor H. Lachos (\href{https://orcid.org/0000-0002-7239-2459}{ORCID})

\bold{Maintainer}: Larissa A. Matos (\email{larissa.amatos@gmail.com}).

}
\references{
  \insertRef{cook1977detection}{StempCens}

  \insertRef{delyon1999convergence}{StempCens}

  \insertRef{louis1982finding}{StempCens}

  \insertRef{zhu2001case}{StempCens}
}
