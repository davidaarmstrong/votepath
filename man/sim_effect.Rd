% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/votepath.r
\name{sim_effect}
\alias{sim_effect}
\title{Simulate Variable Effect}
\usage{
sim_effect(
  obj,
  data,
  varname,
  diffchange = c("unit", "sd"),
  vals = NULL,
  b_var = TRUE,
  R = 100,
  lastMod = c("expected", "draw"),
  ...
)
}
\arguments{
\item{obj}{An object of class \code{votepath}}

\item{data}{A data frame that contains all of the variables form the analysis.}

\item{varname}{The name of a variable whose effect will be evaluated.}

\item{diffchange}{The amount to change \code{varname}.  Changes will be \code{x-.5*diffchange} and \code{x+.5*diffchange}.
For categorical variables, the first and last categories will be chosen.  The \code{diffchange} parameter
defines the change so long as \code{vals=NULL}.  If \code{vals} is not \code{NULL}, then those values will be used for everyone.}

\item{vals}{A vector of length 2 giving the values that will be used to evaluate the effect size.  This will override
\code{diffchange}.  The values must be of the same class as the variable being changed.  For example, if the variable being
changed is a factor, the \code{vals} vector also has to be a factor with the same levels as the variable in \code{varname}.}

\item{b_var}{Logical indicating whether sampling variability on the coefficients should be incorporated in the simulation.}

\item{R}{Number of simulations to be conducted.}

\item{lastMod}{Should the prediction from the last model be a draw or a the expected value?}

\item{...}{Other arguments to be passed down.}
}
\description{
Simulate a variable's effect in a vote path analysis.
}
