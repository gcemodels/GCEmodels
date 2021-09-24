\name{Objective-and-Gradent_Functions}
\alias{GCElin_objFunct}
\alias{GCElin_gradFunct}
\title{Set of functions in example RcppEigen package}
\description{
  These two functions are created for optimization.
}
\usage{
GCElin_objFunct()
GCElin_gradFunct()
}
\arguments{
}
\value{
  \code{GCElin_objFunct()} does not return a value, but displays a
  message to the console.

  \code{GCElin_gradFunct()} returns a numeric matrix computed as the
  outer (vector) product of \code{x}.
  
}
\details{
  These are example functions which should be largely
  self-explanatory. Their main benefit is to demonstrate how to write a
  function using the Eigen C++ classes, and to have to such a
  function accessible from R.
}
\references{
  See the documentation for Eigen, and RcppEigen, for more details.
}
\examples{
  x <- sqrt(1:4)
  rcppeigen_innerproduct(x)
  rcppeigen_outerproduct(x)
}
\author{Marco Sandri}

