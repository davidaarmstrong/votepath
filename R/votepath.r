pGumbel <- function (q, mu = 0, sigma = 1){
  stopifnot(sigma > 0)
  exp(-exp(-((q - mu)/sigma)))
}


#' Take a random draw from a variable variable
#' @param obj Model object of class `lm`, `glm`, `polr` or `multinom`
#' @param data Data frame giving variables used to estimate model
#' @param ret For return a draw from the variable or return the probabilities.
#' @param incl_b_var Logical indicating whether the sampling variable of
#' the coefficients should be incorporated into the simulation.
#' @param ... Other arguments to be passed down, currently not implemented.
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats model.matrix coef vcov rmultinom rbinom rnorm plogis family formula pcauchy pnorm
#' @export
draw_val <- function(obj,
                     data,
                     ret=c("draw", "expected"),
                     incl_b_var = TRUE,
                     ...){
  ret <- match.arg(ret)
  X <- model.matrix(formula(obj), data)
  if(inherits(obj, "multinom")){
    b <- c(t(coef(obj)))
    if(incl_b_var){
      B <- mvrnorm(1, b, vcov(obj))
      B <- cbind(0, matrix(B, ncol=(length(obj$lev)-1)))
    }else{
      B <- cbind(0, t(coef(obj)))
    }
    xb <- X %*% B
    prob <- prop.table(exp(xb), 1)
    tmp <- sapply(1:nrow(prob), function(i)rmultinom(1, 1, prob[i,])  )
    out <- factor(apply(tmp, 2, which.max),
                  levels=1:length(obj$lev),
                  labels=obj$lev)

  }
  if(inherits(obj, "polr")){
    cdf <- switch(obj$method,
                  logistic = plogis,
                  probit = pnorm,
                  cloglog = pGumbel,
                  cauchit = pcauchy)

    X <- X[,-1]
    b <- c(coef(obj), obj$zeta)
    if(incl_b_var){
      B <- mvrnorm(1, b, vcov(obj))
      Z <- B[grep("\\|", names(B))]
      B <- B[-grep("\\|", names(B))]
    }else{
      B <- coef(obj)
      Z <- obj$zeta
    }
    xb <- X %*% B
    Z <- c(-Inf, Z, Inf)
    qs <- sapply(Z, function(z)cdf(z-xb))
    prob <- sapply(2:ncol(qs), function(i)qs[,i] - qs[,(i-1)])
    tmp <- sapply(1:nrow(prob), function(i)rmultinom(1, 1, prob[i,])  )
    out <- factor(apply(tmp, 2, which.max),
                  levels=1:length(obj$lev),
                  labels=obj$lev)
  }
  if(inherits(obj, "glm")){
    if(family(obj)$family == "binomial"){
      if(incl_b_var){
        b <- mvrnorm(1, coef(obj),vcov(obj))
      }else{
        b <- coef(obj)
      }
      prob <- family(obj)$linkinv(X %*% b)
      out <- rbinom(length(prob), 1, prob)
    }else{
      stop(paste0("GLM family ", family(obj)$family, " not currently supported.\n"))
    }
  }
  if(inherits(obj, "lm") & !inherits(obj, "glm")){
    if(incl_b_var){
      b <- mvrnorm(1, coef(obj),vcov(obj))
    }else{
      b <- coef(obj)
    }
    prob <- X %*% b
    out <- rnorm(nrow(X), prob, summary(obj)$sigma)
  }
  if(ret == "draw"){
    out
  }else{
    prob
  }
}

#' Find type of variable
#'
#' Find the variable type to identify the appropriate model for the variable
#'
#' @param x a vector of values to be evaluated
#' @param ... other arguments to be passed down, currently not implemented.
#'
#' @export
find_type <- function(x, ...){
  if(is.numeric(x) & length(unique(na.omit(x))) > 2){
    type <- "lm"
  }
  if(length(unique(na.omit(x))) == 2 & !is.character(x)){
    type <- "glm"
  }
  if(is.factor(x)){
    if(is.ordered(x)){
      type <- "polr"
    }else{
      type <- "multinom"
    }
  }
  type
}


#' Specify a Vote Path Analysis
#'
#' The vote path model estimates independent models that would otherwise use a block-recursive
#' structure.  Each model in each block is fit independently using the appropriate error distribution.
#'
#' @param blocks A list in which each element contains a vector of variable names for variables that will represent
#' that block's effects.  The list may be named or not.
#' @param models If `NULL`, the algorithm will use `find_type()` to find the appropriate model.  This assumes that variables
#' that should be estimated using ordinal regression models are ordered factors.  Otherwise, models can be a named vector of
#' model types for some or all of the variables included in the model.
#' @param data A data frame containing all of the variables in the `blocks` list.
#' @param ... Other arguments to be passed down.
#'
#' @importFrom stats lm glm family binomial
#' @importFrom MASS polr
#' @importFrom nnet multinom
#'
#' @export
vote_path <- function(blocks, models = NULL, data, ...){

  types <- sapply(c(unlist(blocks)), function(nm)find_type(data[[nm]]))
  allvars <- c(unlist(blocks))
  dv <- allvars[length(allvars)]
  if(!is.null(models)){
    types[match(names(models), names(types))] <- models
  }
  res <- vector(mode="list", length=length(blocks)-1)
  for(i in 2:(length(blocks)-1)){
    mlist <- vector(mode="list", length=length(blocks[[i]]))
    for(j in 1:length(blocks[[i]])){
      form <- reformulate(blocks[[(i-1)]], response=blocks[[i]][j])
      arglist <- list(formula = form, data=data)
      if(types[blocks[[i]][j]] %in% c("polr", "multinom")){
        arglist$Hess <- TRUE
      }
      if(types[blocks[[i]][j]] == "glm"){
        arglist$family <- binomial
      }
      mlist[[j]] <- do.call(types[blocks[[i]][j]], arglist)
      }
    res[[(i-1)]] <- mlist
  }
  fullform <- reformulate(allvars[-length(allvars)],
                          response=dv)
  arglist <- list(formula = fullform,
                  data=data)
  if(types[dv] %in% c("polr", "multinom")){
    arglist$Hess <- TRUE
  }
  if(types[dv] == "glm"){
    arglist$family <- binomial
  }
  res[[(length(blocks)-1)]] <- do.call(types[dv], arglist)
  out <- list(models = res, blocks = blocks)
  class(out) <- "votepath"
  out
}


#' Simulate Variable Effect
#'
#' Simulate a variable's effect in a vote path analysis.
#'
#' @param obj An object of class `votepath`
#' @param data A data frame that contains all of the variables form the analysis.
#' @param varname The name of a variable whose effect will be evaluated.
#' @param diffchange The amount to change `varname`.  Changes will be `x-.5*diffchange` and `x+.5*diffchange`.
#' For categorical variables, the first and last categories will be chosen.  The `diffchange` parameter
#' defines the change so long as `vals=NULL`.  If `vals` is not `NULL`, then those values will be used for everyone.
#' @param vals A vector of length 2 giving the values that will be used to evaluate the effect size.  This will override
#' `diffchange`.  The values must be of the same class as the variable being changed.  For example, if the variable being
#' changed is a factor, the `vals` vector also has to be a factor with the same levels as the variable in `varname`.
#' @param b_var Logical indicating whether sampling variability on the coefficients should be incorporated in the simulation.
#' @param R Number of simulations to be conducted.
#' @param ... Other arguments to be passed down.
#'
#' @importFrom progress progress_bar
#' @export
#'
#'
#'
sim_effect <- function(obj, data, varname, diffchange=c("unit", "sd"), vals=NULL, b_var = TRUE, R=100, ...){
  mods <- obj$models
  blocks <- obj$blocks
  out_i <- out_d <- out_t <- NULL
  pb <- progress_bar$new(total = R)
  if(!is.null(vals) & length(vals) != 2)stop("vals must be a vector of length 2\n")
  if(!is.null(vals) & inherits(data[[varname]], "factor") & !inherits(vals, "factor"))stop("vals must have the same class as varname\n")
  which_block <- min(which(sapply(blocks, function(x)max(varname == x)) == 1))
  if(is.factor(data[[varname]]) & is.null(vals)){
    levs <- levels(data[[varname]])
    vals <- factor(c(1,length(levs)),
                   levels=1:length(levs),
                   labels=levs)
  }
  for(k in 1:R){
  new_0 <- new_1 <- data
  if(is.null(vals)){
    delta <- ifelse(diffchange == "sd", sd(data[[varname]], na.rm=TRUE), 1)
    new_0[[varname]] <- new_0[[varname]] - .5*delta
    new_1[[varname]] <- new_1[[varname]] + .5*delta
  }else{
    new_0[[varname]] <- vals[1]
    new_1[[varname]] <- vals[2]
  }
  for(i in (which_block+1):(length(blocks)-1)){
    for(j in 1:length(blocks[[i]])){
      new_0[[blocks[[i]][j]]] <- draw_val(mods[[(i-1)]][[j]], new_0, incl_b_var = b_var)
      new_1[[blocks[[i]][j]]] <- draw_val(mods[[(i-1)]][[j]], new_1, incl_b_var = b_var)
    }
  }

  res0 <- draw_val(mods[[length(mods)]], new_0, ret="expected", incl_b_var = b_var)
  res1 <- draw_val(mods[[length(mods)]], new_1, ret="expected", incl_b_var = b_var)
  total_effect <- res0-res1

  d0 <- draw_val(mods[[length(mods)]], data, ret="expected", incl_b_var = b_var)
  d1 <- draw_val(mods[[length(mods)]], data, ret="expected", incl_b_var = b_var)
  direct_effect <- d1-d0

  indirect_effect <- total_effect - direct_effect

  out_t <- rbind(out_t, colMeans(total_effect))
  out_d <- rbind(out_d, colMeans(direct_effect))
  out_i <- rbind(out_i, colMeans(indirect_effect))
  pb$tick()
  }
  res <- list(total = out_t, direct= out_d, indirect=out_i)
}

