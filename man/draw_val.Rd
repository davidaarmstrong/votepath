% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/votepath.r
\name{draw_val}
\alias{draw_val}
\title{Take a random draw from a variable variable}
\usage{
draw_val(obj, data, ret = c("draw", "expected"), incl_b_var = TRUE, ...)
}
\arguments{
\item{obj}{Model object of class \code{lm}, \code{glm}, \code{polr} or \code{multinom}}

\item{data}{Data frame giving variables used to estimate model}

\item{ret}{For return a draw from the variable or return the probabilities.}

\item{incl_b_var}{Logical indicating whether the sampling variable of
the coefficients should be incorporated into the simulation.}

\item{...}{Other arguments to be passed down, currently not implemented.}
}
\description{
Take a random draw from a variable variable
}
