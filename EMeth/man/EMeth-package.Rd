\name{EMeth-package}
\alias{EMeth-package}
\alias{EMeth}
\docType{package}
\title{
\packageTitle{EMeth}
}
\description{
\packageDescription{EMeth}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{EMeth}
\packageIndices{EMeth}
~~ This package is an implementation of the EMeth algorithm. The function implemented are as follows:
- deconvEM: Implementation of normal version of EMeth algorithm.
- deconvEM_CV: Implementation of cross validation of regularized version of normal EMeth.
- deconvEM_laplace: Implementation of Laplace distribution version of EMeth algorithm.
- deconvEM_CV_laplace: Cross validation of regularized version of laplace EMeth.

Users can specify the W parameters in these functions to specify heterodrastic variance across probes. This is available for both normal distribution and laplace distribution version of EMeth.~~
}
\author{
\packageAuthor{EMeth}

Maintainer: \packageMaintainer{EMeth}
}
\references{
~~ To Be Written~~
}

\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
~~ See pipeline folders. ~~
}
