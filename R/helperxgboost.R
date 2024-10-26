#' Helper functions for xgboost
#'
#' This script contains the utils functions that are used in ReSurv xgboost.
#' In particular, it finds for each contribution the individuals at risk
#' in the tie.
#'
#' @import data.table

risks_in_the_tie <- function(starts_i,
                             stops,
                             stops_i){

  nstops <- length(stops)
  risksets <- vector(mode="list", length=nstops)

  for(i in 1:nstops){

    stop <- stops[i]

    ivec<- ( (starts_i < stop) & (stops_i >= stop))

    risksets[i] <- list(which(ivec))

  }


  return(risksets)
}

#find for each contribution the events in the tie
events_in_the_tie <- function(starts_i,
                              stops,
                              stops_i){

  nstops <- length(stops)
  eventsset<- vector(mode="list", length=nstops)

  for(i in 1:nstops){

    stop <- stops[i]

    ivec<- ( stops_i == stop)

    eventsset[i] <- list(which(ivec))

  }

  return(eventsset)
}


exp_sum_computer <- function(x,ypred){

  sum(exp(ypred[x]))
}


cox_evaluation_metrix <- function(preds,
                                  dtrain){

  risk_sets <- attr(dtrain, 'risk_sets')
  event_sets <- attr(dtrain, 'event_sets')
  efron_c<- attr(dtrain, 'efron_c')
  tieid<- attr(dtrain, 'tieid')

  exp_p_sum <- rep(sapply(risk_sets,FUN=exp_sum_computer, ypred=preds),tieid)
  exp_p_tie <- rep(sapply(event_sets,FUN=exp_sum_computer, ypred=preds),tieid)
  exp_p <- exp(preds)

  r_k <- exp_p_sum-efron_c*exp_p_tie

  lkh<-(exp_p/r_k)

  value <- -sum(log(lkh))
  return(list(metric = "log-partial likelihood", value = value/length(preds) ))
}

##

cox_loss_objective <- function(preds,dtrain){

  Ti <- attr(dtrain, 'truncation')
  Ei <- attr(dtrain, 'claim_arrival')

  risk_sets <- attr(dtrain, 'risk_sets')
  event_sets <- attr(dtrain, 'event_sets')

  efron_c<-attr(dtrain, 'efron_c')
  tieid<- attr(dtrain, 'tieid')
  # obs_tie<- attr(dtrain, 'groups')

  exp_p_sum <- sapply(risk_sets,
                      FUN=exp_sum_computer,
                      ypred=preds)

  exp_p_tie <- sapply(event_sets,
                      FUN=exp_sum_computer,
                      ypred=preds)

  tmp1 <- data.table(risks_s = rep(exp_p_sum, tieid),
                     events_s = rep(exp_p_tie, tieid),
                     efron_c=efron_c,
                     ties = Ei)
  tmp_alpha_i=tmp1[, .(alpha_i = sum(1/(risks_s-efron_c*events_s))), by = ties]$alpha_i
  alpha_i = vector("numeric",length = max(Ei))
  alpha_i[unique(Ei)] = tmp_alpha_i

  tmp_beta_i=tmp1[, .(beta_i = sum(efron_c/(risks_s-efron_c*events_s))), by = ties]$beta_i
  beta_i = vector("numeric",length = max(Ei))
  beta_i[unique(Ei)] = tmp_beta_i

  tmp_gamma_i=tmp1[, .(gamma_i = sum(1/(risks_s-efron_c*events_s)^2)), by = ties]$gamma_i
  gamma_i = vector("numeric",length = max(Ei))
  gamma_i[unique(Ei)] = tmp_gamma_i

  tmp_omega_i=tmp1[, .(omega_i = sum((1-(1-efron_c)^2)/(risks_s-efron_c*events_s)^2)), by = ties]$omega_i
  omega_i = vector("numeric",length = max(Ei))
  omega_i[unique(Ei)] = tmp_omega_i



  exp_p <- exp(preds)
  n <- length(exp_p)
  J <- length(alpha_i)
  alpha_i_lt = vector("numeric",n)
  gamma_i_lt = vector("numeric",n)
  beta_i_lt = vector("numeric",n)
  omega_i_lt = vector("numeric",n)

  for(i in 1:n){

    alpha_i_lt[i]= sum(alpha_i[(Ti[i]+1):Ei[i]])
    beta_i_lt[i] = beta_i[Ei[i]]
    gamma_i_lt[i]= sum(gamma_i[(Ti[i]+1):Ei[i]])
    omega_i_lt[i] = omega_i[Ei[i]]
  }

  #we consider the nll
  grad <- exp_p*(alpha_i_lt-beta_i_lt)-1

  hess <- grad-(exp_p^2)*(gamma_i_lt-omega_i_lt)+1
  return(list(grad=grad,hess=hess))
}


utils::globalVariables(c("risks_in_the_tie","starts_i","stops","stops_i","risks_s","events_s","ties"))



