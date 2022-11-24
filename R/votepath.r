## TODO: Think about mediated effect - should the effects be calculated as they are now
##       where the main variable of interest doesn't change but all the subsequent variables
##       do?  Should it be calculated using the final model, or a model regressing the dependent
##       variable on everything from the block of interest forward, or maybe all past blocks? 

pGumbel <- function (q, mu = 0, sigma = 1){
  stopifnot(sigma > 0)
  exp(-exp(-((q - mu)/sigma)))
}

#' Draw coefficients from posterior distribution 
#' 
#' Draw coefficients from their posterior distributions implied by the model.  
#' @param obj A model object - currently supported models are `lm`, `glm`, `polr`, `multinom`, `svyglm`, `svyolr` and `svymultinom`. 
#' @param ... Other arguments to be passed down, currently not implemented.
#' 
#' @export
draw_coefs <- function(obj, ...){
  UseMethod("draw_coefs")  
}

#' @method draw_coefs lm
#' @export
draw_coefs.lm <- function(obj, ...){
  out <- MASS::mvrnorm(1, coef(obj), vcov(obj))
  out <- matrix(out, nrow=1)
  out
}

#' @method draw_coefs glm
#' @export
draw_coefs.glm <- function(obj, ...){
  out <- MASS::mvrnorm(1, coef(obj), vcov(obj))
  out <- matrix(out, nrow=1)
  out
}

#' @method draw_coefs svyglm
#' @export
draw_coefs.svyglm <- function(obj, ...){
  out <- MASS::mvrnorm(1, coef(obj), vcov(obj))
  out <- matrix(out, nrow=1)
  out
}

#' @method draw_coefs multinom
#' @export
draw_coefs.multinom <- function(obj, ...){
  b <- coef(obj)
  ncb <- ncol(b)
  v <- vcov(obj)
  v <- cbind(matrix(0, nrow=nrow(v), ncol=ncol(b)), v)
  v <- rbind(matrix(0, nrow=ncol(b), ncol = ncol(v)), v)
  b <- rbind(0, b)
  b <- c(t(b))
  out <- MASS::mvrnorm(1, b, v)
  out <- matrix(out, nrow=1)
  matrix(out[1,], ncol=ncb, byrow=TRUE)
}

#' @method draw_coefs polr
#' @export
draw_coefs.polr <- function(obj, ...){
  nb <- length(coef(obj))
  nz <- length(obj$zeta) + 2
  b <- c(-coef(obj), obj$zeta)
  v <- vcov(obj)
  out <- MASS::mvrnorm(1, b, v)
  out <- matrix(out, nrow=1)
  mat <- matrix(rep(out[1, 1:nb], nz), ncol=nz, nrow=nb)
  mat <- rbind(c(-Inf, out[1,(nb+1):ncol(out)], Inf), mat)
  t(mat)
}

#' @method draw_coefs svyolr
#' @export
draw_coefs.svyolr <- function(obj, ...){
  b <- coef(obj)
  zetas <- which(grepl("\\|", names(b)))
  betas <- setdiff(1:length(b), zetas)
  nb <- length(betas)
  nz <- length(zetas) + 2
  v <- vcov(obj)
  out <- MASS::mvrnorm(1, b, v)
  out <- matrix(out, nrow=1)
  mat <- matrix(rep(out[1, 1:nb], nz), ncol=nz, nrow=nb)
  mat <- rbind(c(-Inf, out[1,(nb+1):ncol(out)], Inf), mat)
  t(mat)
}

#' @method draw_coefs svrepstatmisc
#' @export
draw_coefs.svrepstatmisc <- function(obj, ...){
  b <- coef(obj)
  ncb <- length(b)/length(grep("\\(Intercept\\)", names(b)))
  v <- vcov(obj)
  v <- cbind(matrix(0, nrow=nrow(v), ncol=ncb), v)
  v <- rbind(matrix(0, nrow=ncb, ncol = ncol(v)), v)
  b <- c(rep(0, ncb), b)
  out <- MASS::mvrnorm(1, b, v)
  out <- matrix(out, nrow=1)
  matrix(out[1,], ncol=ncb, byrow=TRUE)
}

#' Set-up for Drawing Values
#' 
#' The `prob` function sets up values that will be past `draw_vals`.  For linear models and 
#' gaussian GLMs, these are simply predicted values.  For non-linear models, these are inverse-link 
#' transformed linear predictor values (e.g., predicted probabilities). 
#' @param obj A model object - currently supported models are `lm`, `glm`, `polr`, `multinom`, `svyglm`, `svyolr` and `svymultinom`. 
#' @param b A list of model coefficients generated by `draw_coefs()`.  
#' @param data A data frame that will be used to generate predictions. 
#' @param ... Other arguments to be passed down, currently none implemented. 
#' 
#' @importFrom stats plogis pnorm
#' @export
prob <- function(obj, b, data, ...){
  UseMethod("prob")
}

#' @method prob polr
#' @export
prob.polr <- function(obj, b, data, ...){
  pfun <- switch(obj$method, logistic = plogis, probit = pnorm)
  X <- model.matrix(obj, data=data)
  q <- pfun(X %*% t(b))
  D <- matrix(0, ncol = ncol(q)-1, nrow=ncol(q))
  for(j in 1:ncol(D)){
    D[j,j] <- -1
    D[(j+1),j] <- 1
  }
  q %*% D
}

#' @method prob svyolr
#' @export
prob.svyolr <- function(obj, b, data, ...){
  pfun <- switch(obj$method, logistic = plogis, probit = pnorm)
  X <- model.matrix(formula(obj), data=data$variables)
  q <- pfun(X %*% t(b))
  D <- matrix(0, ncol = ncol(q)-1, nrow=ncol(q))
  for(j in 1:ncol(D)){
    D[j,j] <- -1
    D[(j+1),j] <- 1
  }
  q %*% D
}

#' @method prob lm
#' @export
prob.lm <- function(obj, b, data, ...){
  X <- model.matrix(obj, data=data)
  X %*% t(b)
}

#' @method prob glm
#' @export
prob.glm <- function(obj, b, data, ...){
  X <- model.matrix(obj, data=data)
  obj$family$linkinv(X %*% t(b))
}

#' @method prob svyglm
#' @export
prob.svyglm <- function(obj, b, data, ...){
  X <- model.matrix(obj, data=data$variables)
  obj$family$linkinv(X %*% t(b))
}

#' @method prob multinom
#' @export
prob.multinom <- function(obj, b, data, ...){
  X <- model.matrix(formula(obj), data=data) 
  q <- exp(X %*% t(b))
  prop.table(q, 1)
}

#' @method prob svrepstatmisc
#' @export
prob.svrepstatmisc <- function(obj, b, data, ...){
  X <- model.matrix(attr(obj, "formula"), data=data$variables) 
  q <- exp(X %*% t(b))
  prop.table(q, 1)
}

## TODO: finish draw_val for polr, svyolr, multinom, svymultinom and svyglm
#' Take a random draw from a variable variable
#' @param obj Model object of class `lm`, `glm`, `polr` or `multinom`
#' @param probs A list that was produced by the `prob()` function.  
#' @param e matrix or vector of errors to add to the predicted values to generate the draw in an `lm` or gaussian `glm`.
#' @param ... Other arguments to be passed down, currently not implemented.
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats model.matrix coef vcov rmultinom rbinom rnorm plogis family formula pcauchy pnorm
#' @export
draw_val <- function(obj,
                     probs, 
                     e=NULL, 
                     ...){
  UseMethod("draw_val")
}

#' @method draw_val lm
#' @export
draw_val.lm <- function(obj, 
                        probs,
                        e = NULL, 
                        ...){
  if(is.null(e)){
    sigma <- summary(obj)$sigma
    rnorm(length(probs), probs, sigma)
  }else{
    probs + e    
  }
}

#' @method draw_val glm
#' @export
draw_val.glm <- function(obj, 
                         probs, 
                         e=NULL,
                         ...){
  if(!obj$family$family %in% c("binomial", "gaussian"))stop("Currently only gaussian and binomial GLMs are implemented.\n")
  if(family(obj)$family == "gaussian"){
    if(is.null(e)){
      sigma <- summary(obj)$dispersion
      rnorm(length(probs), probs, sigma)
    }else{
      probs + e
    }
  }else{
    rbinom(length(probs), 1, probs)
  }
}

#'@method draw_val svyglm
#'@export
draw_val.svyglm <- function(obj, 
                            probs, 
                            e=NULL,
                            ...){
  if(!obj$family$family %in% c("binomial", "gaussian"))stop("Currently only gaussian and binomial GLMs are implemented.\n")
  if(obj$family$family == "gaussian"){
    if(is.null(e)){
      sigma <- c(summary(obj)$dispersion)
      rnorm(length(probs), probs, sigma)
    }else{
      probs + e
    }
  }else{
    rbinom(length(probs), 1, probs)
  }
}

#'@method draw_val polr
#'@export
draw_val.polr <- function(obj, probs, e=NULL, ...){
   fac <- apply(probs, 1, function(z)which.max(rmultinom(1, 1, z)))
   fac <- factor(fac, levels=1:length(obj$lev), labels=obj$lev)
   unname(fac)
}

#'@method draw_val svyolr
#'@export
draw_val.svyolr <- function(obj, probs, e=NULL, ...){
  draw_val.polr(obj, probs, ...)  
}

#'@method draw_val multinom
#'@export
draw_val.multinom <- function(obj, probs, e=NULL, ...){
  fac <- apply(probs, 1, function(z)which.max(rmultinom(1, 1, z)))
  fac <- factor(fac, levels=1:length(obj$lev), labels=obj$lev)
  unname(fac)
}

#'@method draw_val svrepstatmisc
#'@export
draw_val.svrepstatmisc <- function(obj, probs, e=NULL, ...){
  fac <- apply(probs, 1, function(z)which.max(rmultinom(1, 1, z)))
  fac <- factor(fac, levels=1:length(attr(obj, "lev")), labels=attr(obj, "lev"))
  unname(fac)
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
#' @importFrom survey svyglm svyolr svydesign as.svrepdesign
#'
#' @export
vote_path <- function(blocks,
                      models = NULL,
                      data, ...){
  orig_bl <- blocks
  blocks <- unname(blocks)
  if(!inherits(data, "survey.design")){
    types <- sapply(c(unlist(unname(blocks))), function(nm)find_type(data[[nm]]))
  }else{
    types <- sapply(c(unlist(unname(blocks))), function(nm)find_type(data$variables[[nm]]))
  }
  if(inherits(data, "survey.design")){
    if(any(types == "lm")){
      types[which(types == "lm")] <- "svyglm_lm"
    }
    if(any(types == "glm")){
      types[which(types == "glm")] <- "svyglm_bin"
    }
    if(any(types == "polr")){
      types[which(types == "polr")] <- "svyolr"
    }
    if(any(types == "multinom")){
      types[which(types == "multinom")] <- "svymultinom"
    }
  }
  if(any(types == "svymultinom")){
    if(!requireNamespace("svrepmisc", quietly=TRUE)){
      stop("Must install and load svrepmisc with e.g., remotes::install_github('carlganz/svrepmisc')\n")
    }
  }
  allvars <- c(unlist(unname(blocks)))
  dv <- allvars[length(allvars)]
  if(!is.null(models)){
    types[match(names(models), names(types))] <- models
  }
  res <- vector(mode="list", length=length(blocks)-1)
  for(i in 2:(length(blocks)-1)){
    mlist <- vector(mode="list", length=length(blocks[[i]]))
    for(j in 1:length(blocks[[i]])){
      form <- reformulate(blocks[[(i-1)]], response=blocks[[i]][j])
      if(!inherits(data, "survey.design")){
        arglist <- list(formula = form, data=data)
      }else{
        arglist <- list(formula = form, design=data)
      }
      if(types[blocks[[i]][j]] %in% c("polr")){
        arglist$Hess <- TRUE
      }
      if(types[blocks[[i]][j]] %in% c("multinom", "svymultinom")){
        arglist$maxit <- 250
        arglist$Hess <- TRUE
        arglist$trace <- FALSE
        if(types[blocks[[i]][j]] == "svymultinom"){
          arglist$design <- as.svrepdesign(arglist$design)
        }
        
      }
      if(types[blocks[[i]][j]] == "glm"){
        arglist$family <- binomial
      }
      if(types[blocks[[i]][j]] == "svyglm_lm"){
        types[blocks[[i]][j]] <- "svyglm"
      }
      if(types[blocks[[i]][j]] == "svyglm_bin"){
        types[blocks[[i]][j]] <- "svyglm"
        arglist$family <- binomial
      }
      mlist[[j]] <- do.call(types[blocks[[i]][j]], arglist)
      if(types[blocks[[i]][j]] == "svymultinom"){
        attr(mlist[[j]], "formula") <- arglist$formula
        attr(mlist[[j]], "lev") <- levels(data$variables[[blocks[[i]][j]]])
      }
      }
    res[[(i-1)]] <- mlist
  }
  fullform <- reformulate(allvars[-length(allvars)],
                          response=dv)
  if(!inherits(data, "survey.design")){
    arglist <- list(formula = fullform,
                    data=data)
  }else{
    arglist <- list(formula = fullform,
                    design=data)
  }
  if(types[dv] %in% c("polr")){
    arglist$Hess <- TRUE
  }
  if(types[dv] %in% c("multinom", "svymultinom")){
    arglist$maxit <- 250
    arglist$Hess <- TRUE
    arglist$trace <- FALSE
    if(types[dv] == "svymultinom"){
      arglist$design <- as.svrepdesign(arglist$design)
    }
  }
  if(types[dv] == "glm"){
    arglist$family <- binomial
  }
  if(types[dv] == "svyglm_lm"){
    types[dv] <- "svyglm"
  }
  if(types[dv] == "svyglm_bin"){
    types[dv] <- "svyglm"
    arglist$family <- binomial
  }
  res[[(length(blocks)-1)]] <- do.call(types[dv], arglist)
  if(types[dv] == "svymultinom"){
    attr(res[[(length(blocks)-1)]], "formula") <- fullform
    attr(res[[(length(blocks)-1)]], "lev") <- levels(data$variables[[dv]])
  }
  out <- list(models = res, blocks = orig_bl)
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
#' @param R Number of simulations to be conducted.
#' @param ... Other arguments to be passed down.
#'
#' @importFrom progress progress_bar
#' @importFrom stats sd
#' @export
sim_effect <- function(obj,
                       data,
                       varname,
                       diffchange=c("unit", "sd"),
                       vals=NULL,
                       R=100,
                       ...){
  mods <- obj$models
  blocks <- obj$blocks
  dv <- blocks[[length(blocks)]]
  if(!is.null(names(blocks))){
    nm_blk <- names(blocks)
  }else{
    nm_blk <- paste0("Block ", 1:length(blocks))
  }
  out_i <- out_d <- out_t <- out_br <-  NULL
  pb <- progress_bar$new(total = R)
  if(!is.null(vals) & length(vals) != 2)stop("vals must be a vector of length 2\n")
  if(!is.null(vals) & inherits(data[[varname]], "factor") & !inherits(vals, "factor"))stop("vals must have the same class as varname\n")
  which_block <- min(which(sapply(blocks, function(x)varname %in% x)))
  #which_block <- min(which(sapply(blocks, function(x)max(varname == x)) == 1))
  if(is.factor(data[[varname]]) & is.null(vals)){
    levs <- levels(data[[varname]])
    vals <- factor(c(1,length(levs)),
                   levels=1:length(levs),
                   labels=levs)
  }
  med <- vector(mode="list", length=length((which_block+1):(length(blocks)-1)))
  te <- de <- ie <- NULL
  for(r in 1:R){
    new_0 <- new_1 <- iedat_0 <- iedat_1 <- data
    if(is.null(vals)){
      delta <- ifelse(diffchange == "sd", sd(data[[varname]], na.rm=TRUE), 1)
      new_0[[varname]] <- new_0[[varname]] - .5*delta
      new_1[[varname]] <- new_1[[varname]] + .5*delta
    }else{
      new_0[[varname]] <- vals[1]
      new_1[[varname]] <- vals[2]
    }
    b_direct <- draw_coefs(mods[[length(mods)]], R=1)
    p0_direct <- prob(mods[[length(mods)]], b_direct, new_0)
    p1_direct <- prob(mods[[length(mods)]], b_direct, new_1)
    if(inherits(mods[[length(mods)]], "lm") & !inherits(mods[[length(mods)]], "glm")){
      e_final <- rnorm(nrow(new_0), 0, summary(mods[[length(mods)]])$sigma)
    }else{
      e_final <- NULL
    }
    if(inherits(mods[[length(mods)]], "glm")){
      if(mods[[length(mods)]]$family$family == "gaussian"){
        e_final <- rnorm(nrow(new_0), 0, summary(mods[[length(mods)]])$dispersion)
      }
    }
    d0_direct <- draw_val(mods[[length(mods)]], p0_direct, e=e_final)
    d1_direct <- draw_val(mods[[length(mods)]], p1_direct, e=e_final)
    if(!is.matrix(d0_direct) & !is.factor(d0_direct)){
      d0_direct <- matrix(d0_direct, ncol=1)
      d1_direct <- matrix(d1_direct, ncol=1)
    }  
    if(is.factor(d0_direct)){
      tab0 <- table(d0_direct)/sum(table(d0_direct))
      tab1 <- table(d1_direct)/sum(table(d1_direct))
    }
    if(is.matrix(d0_direct)){
      de <- rbind(de, colMeans(d1_direct-d0_direct))
    }else{
      de <- rbind(de, tab1-tab0)
    }
    k <- 1
    for(i in (which_block+1):(length(blocks)-1)){
      for(j in 1:length(blocks[[i]])){
        b <- draw_coefs(mods[[(i-1)]][[j]], R=1)
        p0 <- prob(mods[[(i-1)]][[j]], b, new_0)
        p1 <- prob(mods[[(i-1)]][[j]], b, new_1)
        if(inherits(mods[[(i-1)]][[j]], "lm") & !inherits(mods[[(i-1)]][[j]], "glm")){
          e <- rnorm(nrow(new_0), 0, summary(mods[[(i-1)]][[j]])$sigma)
        }else{
          e <- NULL
        }
        if(inherits(mods[[(i-1)]][[j]], "glm")){
          if(mods[[(i-1)]][[j]]$family$family == "gaussian"){
            e <- rnorm(nrow(new_0), 0, summary(mods[[(i-1)]][[j]])$dispersion)
          }
        }
        new_0[[blocks[[i]][j]]] <- iedat_0[[blocks[[i]][j]]] <-draw_val(mods[[(i-1)]][[j]], p0, e=e)
        new_1[[blocks[[i]][j]]] <- iedat_1[[blocks[[i]][j]]] <-draw_val(mods[[(i-1)]][[j]], p1, e=e)
      }
      p0_med <- prob(mods[[length(mods)]], b_direct, iedat_0)
      p1_med <- prob(mods[[length(mods)]], b_direct, iedat_1)
      tmp0_med <- draw_val(mods[[length(mods)]], p0_med, e=e_final)
      tmp1_med <- draw_val(mods[[length(mods)]], p1_med, e=e_final)
      if(!is.matrix(tmp0_med) & !is.factor(tmp0_med)){
        tmp0_med <- matrix(tmp0_med, ncol=1)
        tmp1_med <- matrix(tmp1_med, ncol=1)
      }
      if(is.factor(tmp0_med)){
        tab0_med <- table(tmp0_med)/sum(table(tmp0_med))
        tab1_med <- table(tmp1_med)/sum(table(tmp1_med))
      }
      if(is.matrix(tmp0_med)){
        med[[k]] <- rbind(med[[k]], colMeans(tmp1_med-tmp0_med))
      }else{
        med[[k]] <- rbind(med[[k]], tab1_med-tab0_med)
      }
      names(med)[k] <- nm_blk[i]
      k <- k+1
    }
    

    ## Total Effect
    p0_final <- prob(mods[[length(mods)]], b_direct, new_0)
    p1_final <- prob(mods[[length(mods)]], b_direct, new_1)
    tmp0 <- draw_val(mods[[length(mods)]], p0_final, e=e_final)
    tmp1 <- draw_val(mods[[length(mods)]], p1_final, e=e_final)
    if(!is.matrix(tmp0) & !is.factor(tmp0)){
      tmp0 <- matrix(tmp0, ncol=1)
      tmp1 <- matrix(tmp1, ncol=1)
    }
    if(is.factor(tmp0)){
      tab0 <- table(tmp0)/sum(table(tmp0))
      tab1 <- table(tmp1)/sum(table(tmp1))
    }
    if(is.matrix(tmp0)){
      te <- rbind(te, colMeans(tmp1-tmp0))
    }else{
      te <- rbind(te, tab1-tab0)
    }

    ## Calculated Indirect effect
    p0i_final <- prob(mods[[length(mods)]], b_direct, iedat_0)
    p1i_final <- prob(mods[[length(mods)]], b_direct, iedat_1)
    tmp0i <- draw_val(mods[[length(mods)]], p0i_final, e=e_final)
    tmp1i <- draw_val(mods[[length(mods)]], p1i_final, e=e_final)
    if(!is.matrix(tmp0i) & !is.factor(tmp0i)){
      tmp0i <- matrix(tmp0i, ncol=1)
      tmp1i <- matrix(tmp1i, ncol=1)
    }
    if(is.factor(tmp0i)){
      tab0i <- table(tmp0i)/sum(table(tmp0i))
      tab1i <- table(tmp1i)/sum(table(tmp1i))
    }
    if(is.matrix(tmp0i)){
      ie <- rbind(ie, colMeans(tmp1i-tmp0i))
    }else{
      ie <- rbind(ie, tab1i-tab0i)
    }
    pb$tick()
  }
  ie_infer <- te-de
  de_infer <- te-ie
  te_infer <- ie+de
  if("lev" %in% names(obj$models[[length(obj$models)]])){
    cnms <- obj$models[[length(obj$models)]]$lev
    colnames(de) <- colnames(te) <- colnames(ie) <- 
      colnames(de_infer) <- colnames(te_infer) <- colnames(ie_infer) <- cnms
    med <- lapply(med, function(x){
      colnames(x) <- cnms
      x
    })
  }
  
  res <- list(total = te, direct= de, indirect=ie, total_infer = te_infer, direct_infer = de_infer, indirect_infer=ie_infer, med=med)
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
##' @importFrom magrittr `%>%`
##' @importFrom dplyr mutate bind_rows
##' @importFrom stats quantile
##' 
##' @export
##'
##' @method summary simeff
summary.simeff <- function(object, ..., conf.level=.95){
  a1 <- (1-conf.level)/2
  a2 <- 1-a1
  total_sum = apply(object$total, 2, function(x)c(mean(x), unname(quantile(x, c(a1, a2)))))
  rownames(total_sum) <- c("Mean", "Lower", "Upper")
  total_sum <- as_tibble(t(total_sum), rownames="DV")  %>% 
    mutate(effect = "Total (calculated)")

  direct_sum = apply(object$direct, 2, function(x)c(mean(x), unname(quantile(x, c(a1, a2)))))
  rownames(direct_sum) <- c("Mean", "Lower", "Upper")
  direct_sum <- as_tibble(t(direct_sum), rownames="DV") %>% 
    mutate(effect = "Direct (calculated)")

  indirect_sum = apply(object$indirect, 2, function(x)c(mean(x), unname(quantile(x, c(a1, a2)))))
  rownames(indirect_sum) <- c("Mean", "Lower", "Upper")
  indirect_sum <- as_tibble(t(indirect_sum), rownames="DV") %>% 
    mutate(effect = "Indirect (calculated)")

  total_sum_infer = apply(object$total_infer, 2, function(x)c(mean(x), unname(quantile(x, c(a1, a2)))))
  rownames(total_sum_infer) <- c("Mean", "Lower", "Upper")
  total_sum_infer <- as_tibble(t(total_sum_infer), rownames="DV") %>% 
    mutate(effect = "Total (inferred)")
  
  direct_sum_infer = apply(object$direct_infer, 2, function(x)c(mean(x), unname(quantile(x, c(a1, a2)))))
  rownames(direct_sum_infer) <- c("Mean", "Lower", "Upper")
  direct_sum_infer <- as_tibble(t(direct_sum_infer), rownames="DV") %>% 
    mutate(effect = "Direct (inferred)")
  
  indirect_sum_infer = apply(object$indirect_infer, 2, function(x)c(mean(x), unname(quantile(x, c(a1, a2)))))
  rownames(indirect_sum_infer) <- c("Mean", "Lower", "Upper")
  indirect_sum_infer <- as_tibble(t(indirect_sum_infer), rownames="DV") %>% 
    mutate(effect = "Indirect (inferred)")

  med_sum <- lapply(seq_along(object$med), function(m){
    z = apply(object$med[[m]], 2, function(x)c(mean(x), unname(quantile(x, c(a1, a2)))))
    rownames(z) <- c("Mean", "Lower", "Upper")
    z <- as_tibble(t(z), rownames="DV") %>% 
      mutate(effect = paste0("Mediated (", names(object$med)[m], ")"))
    z})  
  names(med_sum) <- names(object$med)
  
  res <- list(
    total_c = total_sum, 
    direct_c = direct_sum, 
    indirect_c = indirect_sum, 
    total_i = total_sum_infer, 
    direct_i = direct_sum_infer, 
    indirect_i = indirect_sum_infer, 
    mediated = med_sum)
  class(res) <- "summary.simeff"
  res
}

#' Print method for summary.simeff objects. 
#' 
#' Print method for the output from summary.simeff. 
#' 
#' @param obj Object of class `summary.simeff`. 
#' @param ... Other arguments to be passed down, currently unimplemeneted. 
#' 
#' @method print summary.simeff
#' @export
print.summary.simeff <- function(obj, 
                                 ..., 
                                 total_effect = c("calculated", "inferred"), 
                                 direct_effect = c("calculated", "inferred"), 
                                 indirect_effect = c("inferred", "calculated"), 
                                 print_mediated = FALSE){
  te <- match.arg(total_effect)
  de <- match.arg(direct_effect)
  ie <- match.arg(indirect_effect)
  te_type <- ifelse(te == "calculated", "total_c", "total_i")
  de_type <- ifelse(de == "calculated", "direct_c", "direct_i")
  ie_type <- ifelse(ie == "calculated", "indirect_c", "indirect_i")
  t_eff <- obj[[te_type]]
  d_eff <- obj[[de_type]]
  i_eff <- obj[[ie_type]]
  cat("Total Effect\n")
  print(t_eff)
  cat("Direct Effect\n")
  print(d_eff)
  cat("Indirect Effect\n")
  print(i_eff)
  if(print_mediated){
    cat("Mediated Effects\n")
    print(obj$mediated)
  }
}


