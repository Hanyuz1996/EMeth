\name{cv.emeth}
\alias{cv.emeth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
cv.emeth: cross validation for emeth.
}
\description{
An implementation of cross validation for the ridge penalty parameter in EMeth.
}
\usage{
cv.emeth(Y, eta, mu, aber, V, init = "default", nu = 0, family = "laplace", folds = 5, usesubset = TRUE, maxiter = 50, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  Most can be found in emeth. The different parameters are
  \item{folds}{
     Specify the number of folds for the cross validation.
}
  \item{usesubset}{
     Logic variable, if it is true, a random sampled subset of all probes are used to perform cross-validation.
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

  \item{result }{The result of EMeth algorithm after fitting the best selected parameter.}
  \item{choosenu }{The chosen parameter of nu.}
  \item{losslist}{A matrix saving the loss for each fold and each choice of nu.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Hanyu Zhang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## See examples folder.
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
