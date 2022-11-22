sim_effect <- function(obj,
                       data,
                       varname,
                       diffchange=c("unit", "sd"),
                       vals=NULL,
                       b_var = TRUE,
                       R=100,
                       ...){
  mods <- obj$models
  blocks <- obj$blocks
  dv <- blocks[[length(blocks)]]
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
  dv_type <- find_type(data[[dv]])
  #  md0 <- md1 <- data
  for(r in 1:R){
    new_0 <- new_1 <-  data
    if(is.null(vals)){
      delta <- ifelse(diffchange == "sd", sd(data[[varname]], na.rm=TRUE), 1)
      new_0[[varname]] <- new_0[[varname]] - .5*delta
      new_1[[varname]] <- new_1[[varname]] + .5*delta
    }else{
      new_0[[varname]] <- vals[1]
      new_1[[varname]] <- vals[2]
    }
    k <- 1
    for(i in (which_block+1):(length(blocks)-1)){
      for(j in 1:length(blocks[[i]])){
        b <- draw_coefs(models[[i]][[j]], R=1)
        p0 <- prob(mods[[i]][[j]], b, new0)
        p1 <- prob(mods[[i]][[j]], b, new1)
        if(inherits(mods[[i]][[j]], "lm")){
          e <- rnorm(nrow(new0), 0, summary(mods[[i]][[j]]$sigma))
        }else{
          e <- NULL
        }
        if(inherits(mods[[i]][[j]], "glm")){
          if(mods[[i]][[j]]$family$family == "gaussian"){
            e <- rnorm(nrow(new0), 0, summary(mods[[i]][[j]])$dispersion)
          }
        }
        new_0[[blocks[[i]][j]]] <- draw_val(models[[i]][[j]], p0, e=e)
        new_1[[blocks[[i]][j]]] <- draw_val(models[[i]][[j]], p1, e=e)
      }
    }
      ## Should we use md0 and md1 or new0 and new1 in this part of the prediction?
    b_final <- draw_coefs(mods[[legnth(mods)]], R=1)
    p0_final <- prob(mods[[length(mods)]], b_final, new0)
    p1_final <- prob(mods[[length(mods)]], b_final, new1)
    if(inherits(mods[[length(mods)]], "lm")){
      e <- rnorm(nrow(new0), 0, summary(mods[[length(mods)]]$sigma))
    }else{
      e <- NULL
    }
    if(inherits(mods[[length(mods)]], "glm")){
      if(mods[[length(mods)]]$family$family == "gaussian"){
        e <- rnorm(nrow(new0), 0, summary(mods[[length(mods)]])$dispersion)
      }
    }
    tmp0 <- draw_val(mods[[length(mods)]], p0_final, e=e)
    tmp1 <- draw_val(mods[[length(mods)]], p1_final, e=e)
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
