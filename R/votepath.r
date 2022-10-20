pGumbel <- function (q, mu = 0, sigma = 1){
  stopifnot(sigma > 0)
  exp(-exp(-((q - mu)/sigma)))
}


#' Generate draws from the distribution of coefficients
#' @param x model object whose coefficients are to be drawn
#' @param ... Other arguments to be passed down 
#' @param n number of draws to take
#' @importFrom MASS mvrnorm
#' 
#' @export
draw_coef <- function(x, ..., n=1000){
  UseMethod("draw_coef")
}

#' @method draw_coef lm
#' @export
draw_coef.lm <- function(x, ..., n=1000){
  args <- list(...)
  if("var" %in% names(args)){
    v <- args$var
  }else{
    v <- vcov(x)
  }
  B <- mvrnorm(n, coef(x), v)
  list(B=array(t(B), dim=c(length(coef(x)), 1, n)))
}

#' @method draw_coef glm
#' @export
draw_coef.glm <- function(x, ..., n){
  draw_coef.lm(x, ..., n)
}

#' @method draw_coef multinom
#' @export
draw_coef.multinom <- function(x, ..., n){
  args <- list(...)
  if("var" %in% names(args)){
    v <- args$var
  }else{
    v <- vcov(x)
  }
  b <- c(t(coef(x)))
  b0 <- matrix(0, ncol = ncol(coef(x)), nrow=n)
  B <- mvrnorm(n, b, v)
  B <- cbind(b0, B)
  B <- lapply(1:nrow(B), \(i)matrix(B[i,], ncol=ncol(coef(x)), byrow=TRUE))
  B_arr <- array(dim=c(nrow(coef(x))+1, ncol(coef(x)), n))
  for(i in 1:length(B)){
    B_arr[,,i] <- B[[i]]    
  }  
  list(B = B_arr)
}

#' @method draw_coef polr
#' @export
draw_coef.polr <- function(x, ..., n=1000){
  b <- c(coef(x), x$zeta)
  B <- mvrnorm(n, b, vcov(x))
  Z <- B[,grep("\\|", colnames(B)), drop=FALSE]
  B <- B[,-grep("\\|", colnames(B)), drop=FALSE]
  
  B <- array(c(t(B)), dim=c(ncol(B), 1, n))
  Z <- array(c(t(Z)), dim=c(ncol(Z), 1, n))
  list(B=B, Z=Z)
}


#' Take a random draw from a variable variable
#' @param obj Model object of class `lm`, `glm`, `polr` or `multinom`
#' @param data Data frame giving variables used to estimate model
#' @param n Array of coefficients to be used in the draw. 
#' @param ... Other arguments to be passed down, currently not implemented.
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats model.matrix coef vcov rmultinom rbinom rnorm plogis family formula pcauchy pnorm
#' @export
draw_val <- function(obj,
                     data,
                     ...){
  UseMethod("draw_val")
}

#' @method draw_val lm
#' @export
draw_val.lm <- function(obj, 
                        data, 
                        ...){
  X <- model.matrix(formula(obj), data)
  xb <- X %*% coef(obj)
  s2e <- sum(obj$residuals^2)/obj$df.residual
  rnorm(length(xb), xb, sqrt(s2e))
}

#' @method draw_val glm
#' @export
draw_val.glm <- function(obj, 
                         data, 
                         ...){
  X <- model.matrix(formula(obj), data)
  xb <- X %*% coef(obj)
  ev <- obj$family$linkinv(xb)
  if(obj$family$family == "binomial"){
    rfun <- function(n,x)rbinom(n,1,x)
  }
  if(obj$family$family == "poisson"){
    rfun <- function(n,x)rpois(n, x)
  }
  if(!obj$family$family %in% c("binomial", "poisson"))stop("Currently only binomial and poisson families are implemented.\n")
  rfun(length(ev), ev)
}

#' @method draw_val multinom
#' @export
draw_val.multinom <- function(obj, 
                              data, 
                              ...){
  X <- model.matrix(obj, data)
  b <- coef(obj)
  b <- rbind(0, b)
  xb <- X %*% t(b)
  ev <- exp(xb)
  ev <- prop.table(ev, 1)
  draws <- t(apply(ev, 1, \(x)rmultinom(1, 1, x)))
  draws <- apply(draws, 1, which.max)
  factor(draws, levels=1:length(obj$lev), labels=obj$lev)
}
                              
#' @method draw_val polr
#' @export
draw_val.polr <- function(obj, 
                              data, 
                              ...){
  cdf <- switch(obj$method, logistic = plogis, probit = pnorm)
  X <- model.matrix(obj, data)[,-1, drop=FALSE]
  b <- coef(obj)
  xb <- X %*% b
  z <- c(-Inf, obj$zeta, Inf)
  qval <- sapply(z, \(tau)cdf(tau - xb))
  ev <- sapply(2:ncol(qval), function(i)qval[,i] - qval[,(i-1)])
  draws <- t(apply(ev, 1, \(x)rmultinom(1, 1, x)))
  draws <- apply(draws, 1, which.max)
  factor(draws, levels=1:length(obj$lev), labels=obj$lev)
}


#' Find type of variable
#'
#' Find the variable type to identify the appropriate model for the variable
#'
#' @param x a vector of values to be evaluated
#' @param ... other arguments to be passed down, currently not implemented.
#'
#' @importFrom stats na.omit
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
#' @importFrom stats lm glm family binomial reformulate
#' @importFrom MASS polr
#' @importFrom nnet multinom
#'
#' @export
vote_path <- function(blocks,
                      models = NULL,
                      data, ...){

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
      if(types[blocks[[i]][j]] == "multinom"){
        arglist$maxit <- 250
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
  if(types[dv] == "multinom"){
    arglist$maxit <- 250
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
#' @param lastMod Should the prediction from the last model be a draw or a the expected value?
#' @param ... Other arguments to be passed down.
#'
#' @importFrom progress progress_bar
#' @importFrom stats sd
#' @export
#'
#'
#'
sim_effect <- function(obj,
                       data,
                       varname,
                       diffchange=c("unit", "sd"),
                       vals=NULL,
                       b_var = TRUE,
                       R=100,
                       lastMod = c("expected", "draw"),
                       ...){
  mods <- obj$models
  blocks <- obj$blocks
  lastMod <- match.arg(lastMod)
  dv <- blocks[[length(blocks)]]
  out_i <- out_d <- out_t <- out_br <-  NULL
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
  med <- vector(mode="list", length=length((which_block+1):(length(blocks)-1)))
  brvars <- c(unlist(blocks[1:which_block]))
  br_form <- reformulate(brvars, response=dv)
  dv_type <- find_type(data[[dv]])
  br_args <- list(formula = br_form, data=data)
  if(dv_type %in% c("polr", "multinom")){
    br_args$Hess <- TRUE
  }
  if(dv_type == "multinom"){
    br_args$maxit <- 250
  }
  if(dv_type == "glm"){
    br_args$family <- binomial
  }
  br_mod <- do.call(dv_type, br_args)

  for(r in 1:R){
  new_0 <- new_1 <- br_0 <- br_1 <- data
  if(is.null(vals)){
    delta <- ifelse(diffchange == "sd", sd(data[[varname]], na.rm=TRUE), 1)
    new_0[[varname]] <- br_0[[varname]] <- new_0[[varname]] - .5*delta
    new_1[[varname]] <- br_1[[varname]] <- new_1[[varname]] + .5*delta
  }else{
    new_0[[varname]] <- br_0[[varname]] <-vals[1]
    new_1[[varname]] <- br_1[[varname]] <-vals[2]
  }
  md0 <- md1 <- data
  k <- 1
  for(i in (which_block+1):(length(blocks)-1)){
    for(j in 1:length(blocks[[i]])){
      new_0[[blocks[[i]][j]]] <- md0[[blocks[[i]][j]]] <- draw_val(mods[[(i-1)]][[j]], new_0)
      new_1[[blocks[[i]][j]]] <- md1[[blocks[[i]][j]]] <- draw_val(mods[[(i-1)]][[j]], new_1)
    }
    tmp0 <- draw_val(mods[[length(mods)]], md0)
    tmp1 <- draw_val(mods[[length(mods)]], md1)
    if(!is.matrix(tmp0))tmp0 <- matrix(tmp0, ncol=1)
    if(!is.matrix(tmp1))tmp1 <- matrix(tmp1, ncol=1)
    med[[k]] <- rbind(med[[k]], colMeans(tmp1-tmp0))
    k <- k+1
  }

  res0 <- draw_val(mods[[length(mods)]], new_0)
  res1 <- draw_val(mods[[length(mods)]], new_1)
  if(!is.matrix(res0))res0 <- matrix(res0, ncol=1)
  if(!is.matrix(res1))res1 <- matrix(res1, ncol=1)
  total_effect <- res1-res0

  d0 <- draw_val(mods[[length(mods)]], br_0)
  d1 <- draw_val(mods[[length(mods)]], br_1)
  if(!is.matrix(d0))d0 <- matrix(d0, ncol=1)
  if(!is.matrix(d1))d1 <- matrix(d1, ncol=1)

  direct_effect <- d1-d0

  indirect_effect <- total_effect - direct_effect

  be0 <- draw_val(br_mod, br_0)
  be1 <- draw_val(br_mod, br_1)
  if(!is.matrix(be0))be0 <- matrix(be0, ncol=1)
  if(!is.matrix(be1))be1 <- matrix(be1, ncol=1)

  out_br <- rbind(out_br, colMeans(be1-be0))

  out_t <- rbind(out_t, colMeans(total_effect))
  out_d <- rbind(out_d, colMeans(direct_effect))
  out_i <- rbind(out_i, colMeans(indirect_effect))
  pb$tick()
  }

  cnms <- obj$models[[length(obj$models)]]$lev
  colnames(out_t) <- colnames(out_d) <- colnames(out_i) <- colnames(out_br) <- cnms
  med <- lapply(med, function(x){
    colnames(x) <- cnms
    x
  })

  res <- list(total = out_t, direct= out_d, indirect=out_i, br = out_br, mediated=med)
  class(res) <- "simeff"
  return(res)
}


##' Summary method for VotePath Simulated Effects
##'
##' @description summary method for objects of class \code{simeff}
##' @param object Object of class \code{simeff}
##' @param ... Other arguments, currently unimplemented
##' @param conf.level Level at which to make the confidence interval
##'
##' @importFrom tibble as_tibble
##' @export
##'
##' @method summary simeff
summary.simeff <- function(object, ..., conf.level=.95){
  a1 <- (1-conf.level)/2
  a2 <- 1-a1
  total_sum = apply(object$total, 2, function(x)c(mean(x), unname(quantile(x, c(a1, a2)))))
  total_sum <- as_tibble(t(total_sum), rownames="DV") %>% setNames(c("DV", "Mean", "Lower", "Upper"))

  direct_sum = apply(object$direct, 2, function(x)c(mean(x), unname(quantile(x, c(a1, a2)))))
  direct_sum <- as_tibble(t(direct_sum), rownames="DV") %>% setNames(c("DV", "Mean", "Lower", "Upper"))

  indirect_sum = apply(object$indirect, 2, function(x)c(mean(x), unname(quantile(x, c(a1, a2)))))
  indirect_sum <- as_tibble(t(indirect_sum), rownames="DV") %>% setNames(c("DV", "Mean", "Lower", "Upper"))

  cat("Total Effects:\n")
  print(total_sum)
  cat("\nIndirect Effects:\n")
  print(indirect_sum)
  cat("\nDirect Effects:\n")
  print(direct_sum)
}




