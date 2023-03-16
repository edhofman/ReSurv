pkg.env <- new.env()

pkg.env$check.all.present <- function(x,check.on){

  "
  This function checks that you have all the periods in the data,
  from the minimum record to the maximum record.

  "

  tmp <- x

  v <- diff(as.integer(sort(unique(x))))

  if(sum(v>1)>0){

    warning(paste0("Some", check.on, "are missing in the data"))

  }

}

pkg.env$check.time.units <- function(input_time_unit,
                             output_time_unit){

  "
  This function checks that the input output time transformation is consistent.
  E.g. you cannot turn quarters into semesters, you can instead turn trimesters into semesters.

  "

  if((1/input_time_unit)%%(1/output_time_unit) != 0){

    stop('The provided time intervals are not subsettable.')

  }

}

pkg.env$encode.variables <- function(x){
  "
  This function encodes the periods.
  We impose that the indexization starts from 1.

  "
  seq <- min(x):max(x)

  dim = length(seq)

  v <- 1:length(seq)

  v[match(x,seq)]


}


pkg.env$formula.editor <- function(continuous_features,
                           categorical_features,
                           continuous_features_spline,
                           degrees_of_freedom,
                           input_output='i'){
  "
  This util edits creates the string that is used for model fitting in a compact way.
  @param continuous_features: string vector of continuous features to be included in the linear predictor.
  @param categorical_features: string vector of categorical features to be included in the linear predictor.
  @param continuous_features_spline: logical value, T if a spline is added to model the continuous features.
  @param degrees_of_freedom: degrees of freedom of the spline.
  @param input_output: set to input ('i') or output ('o') depending on the formula that we require.
  "


  tmp.cat <- switch(!is.null(categorical_features), paste(categorical_features, collapse='+'), NULL)
  tmp.cont <- switch(!is.null(continuous_features), paste(continuous_features, collapse='+'), NULL)
  tmp.splines <- switch((!is.null(continuous_features) & continuous_features_spline),paste0("pspline(",continuous_features, ",degree=3,df=",degrees_of_freedom,")"),NULL)

  string_formula<- paste(paste0("survival::Surv","(TR_",input_output,", DP_rev_",input_output,", I) ~ "),paste(c(tmp.cat,tmp.cont,tmp.splines), collapse='+'))
  string_formula


}


"This is a vectorized version of the grepl function.
See the grepl function documentation."
pkg.env$vgrepl <- Vectorize(grepl, vectorize.args = "pattern")


pkg.env$model.matrix.creator <- function(data,
                                 select_columns){
  "
  This function encodes the matrices that we need for model fitting.

  "

  #individual_data$training.data
  X <- data %>%
    fastDummies::dummy_cols(select_columns = select_columns, #individual_data$categorical_features
                            remove_selected_columns = T)

  tmp.cond=as.logical(apply(pkg.env$vgrepl(pattern=select_columns,
                                   x=colnames(X)), #individual_data$categorical_features
                            MARGIN=1,
                            sum))

  X <- X[,tmp.cond]

  return(X)

}


pkg.env$model.matrix.extract.hazard.names <- function(X,
                                                      string_formula,
                                                      data){

  formula_ct <- as.formula(string_formula)
  Y<-model.extract(model.frame(formula_ct, data=data),"response")

  enter <- Y[, 1]
  exit <- Y[, 2]
  event <- Y[, 3] != 0
  sco <- exp(rep(0, nrow(Y)))

  time <- sort(seq(1,max(exit[event]), by=1)) #might be times with no events

  X_unique <- unique(X)

  names_hazard <- (data.frame(X_unique) %>%
                     rowwise() %>%
                     mutate(name = paste0(names(.)[c_across() == 1], collapse = ',')))$name

  return(list(enter=enter,
              exit=exit,
              event=event,
              sco=sco,
              time=time,
              names_hazard=names_hazard))


}


# Fitting routines

pkg.env$fit_cox_model <- function(data,
                          formula_ct,
                          X,
                          X_ams){
  "This function is the fitting routine for the cox model."

  cox <- coxph(formula_ct, data=data)
  beta_2 <- cox$coef
  beta<-c(0,beta_2)


  Xb <- as.matrix(X)%*%beta

  X_ams <- cbind(X_ams, Xb)

  beta_ams = unique(round(X_ams,10) )[,ncol(X_ams)] #if no round some systems has too high precision.


  list(
    cox=cox,
    beta=beta,
    beta_2=beta_2,
    Xb=Xb,
    X_ams=X_ams,
    beta_ams=beta_ams
  )
}



# Handling the baseline

pkg.env$hazard_baseline_model <- function(data,
                                  cox,
                                  hazard=NULL,
                                  baseline,
                                  conversion_factor,
                                  nk=50,
                                  nbin=48,
                                  phi=1){


  if(baseline == "breslow"){

    bs_hazard <- basehaz(cox, centered=FALSE) %>%
      mutate(hazard = hazard-lag(hazard,default=0))

    bs_hazard2 = tibble(DP_rev_i = bs_hazard$time,
                        hazard=bs_hazard$hazard) %>%
      mutate(hazard = hazard-lag(hazard, default=0))
  }

  if(baseline == "spline"){
    bs_hazard=bshazard::bshazard(pkg.env$formula.editor(continuous_features=NULL,
                                                categorical_features="1",
                                                continuous_features_spline=F),
                                 data=data[(data$AP_i-1)%%(conversion_factor^-1)==0 & data$claim_type==0,],
                                 nk=nk,
                                 nbin=nbin,
                                 phi=phi)
    bs_hazard2 <- tibble(DP_rev_i = bs_hazard$time,
                         hazard = bs_hazard$hazard)
  }

  return(list(bs_hazard=bs_hazard,
              bs_hazard2=bs_hazard2))

}

# Hazard computation

pkg.env$hazard_f<-function(i,
                   enter,
                   time,
                   exit,
                   event){
  "
  This function computes the hazard over a matrix.

  "
  rs <- (enter < time[i] & exit >= time[i])

  O<- sum(event[exit == time[i]])
  E <- sum(rs,-1/2*event[exit == time[i]] )

  c(O/E)}

pkg.env$dissect_hazard_name <- function(names_hazard, name = "AP"){

  "
  Get accident period and covariates to use for later grouping.

  "

  names_hazard <- as.character(names_hazard)
  start_position <- gregexpr('_', names_hazard)[[1]][2] +1
  end_position <- gregexpr(',', names_hazard)[[1]][1] -1

  AP <- substr(names_hazard, start_position, end_position )

  covariate <- substr(names_hazard, end_position+2, nchar(names_hazard) )
  if(name == "AP"){
    return(
      as.numeric(AP)
    )
  }
  if(name=="covariate"){
    return(
      covariate
    )
  }
}

pkg.env$hazard_data_frame <- function(hazard, conversion_factor){

  "
  Convert hazard matrix to dataframe and add grouping variables.

  "

  hazard <- as.data.frame(hazard)

  hazard$DP_rev_i <- 1:nrow(hazard)

  hazard_t <- reshape2::melt(hazard, id.vars = "DP_rev_i")

  hazard_t$AP_i <- sapply(hazard_t$variable, pkg.env$dissect_hazard_name)
  hazard_t$covariate <- mapply(pkg.env$dissect_hazard_name, hazard_t$variable,  MoreArgs = list(name= "covariate"))
  hazard_t$AP_o  <- ceiling(hazard_t$AP_i*conversion_factor)

  #Group is pr. covariate, output accident period
  groups <- unique(data.frame(AP_o = hazard_t$AP_o, covariate = hazard_t$covariate)) %>%
    mutate(group = row_number())

  hazard_group <- hazard_t %>%  left_join(groups, by=c("AP_o", "covariate")) %>%
    select(DP_rev_i, AP_i, value, group)

  return(
    hazard_group
  )
}

pkg.env$m_to_q_hazard <- function(i,
                          hazard,
                          frame_tmp,
                          conversion_factor){
  "
  Group monthly hazard to quarterly.

  "
  inv_conv_factor <- 1/conversion_factor #Inverse conversion factor corresponds to width of output development period measured in input development periods

  h_tmp <- data.frame(hazard[,((i-1)*inv_conv_factor+2 ):(i*inv_conv_factor+1)]) %>% # extract hazard to be grouped
    mutate(DP_rev_i = row_number()) %>%
    filter(DP_rev_i<=max(frame_tmp$time_w)) %>% #Filter to only have hazard values that are used in output development period
    filter(DP_rev_i>=min(frame_tmp$time_w))

  mask_frame <- is.na(frame_tmp2[,c(2:ncol(frame_tmp2),1)]) #create a frame holding 1/0 for relevant values in said parallelogram. Rearrange some of the columns to get correct ordering
  h_tmp[mask_frame] <- 0 #replace values, that doesn't appear in parallelogram with 0

  h_tmp[,1:inv_conv_factor] <- cumsum(h_tmp[,1:inv_conv_factor])

  h_tmp[mask_frame] <- 0 #replace values, that doesn't appear in parallelogram with 0

  #Joint the weights to the hazard and apply the weighing to get a single estimate for the entire development period
  h_tmp2 <- h_tmp %>%  left_join(frame_tmp2, by=c("DP_rev_i"="time_w"))

  hazard2<- sum(h_tmp2[,1:inv_conv_factor]*h_tmp2[,(inv_conv_factor+2):ncol(h_tmp2)], na.rm=T)/sum(h_tmp2[,(inv_conv_factor+2):ncol(h_tmp2)],na.rm=T)

  return(hazard2)
}

pkg.env$m_to_q_hazard_2 <- function(i,
                                  hazard_data_frame,
                                  frame_tmp,
                                  frame_tmp2,
                                  conversion_factor){
  "
  Group monthly hazard to quarterly.

  "

  #select relevant hazard values, and create p_month variable
  h_tmp_2 <- hazard_data_frame %>%  filter(group==i) %>%
    filter(DP_rev_i<=max(frame_tmp$time_w)) %>% #Filter to only have hazard values that are used in output development period
    filter(DP_rev_i>=min(frame_tmp$time_w)) %>%
    mutate(p_month = (AP_i-1)%%(1/(conversion_factor))+1)

  #Join the weights, which also give an indication of which development periods, pr. accident period is relevant in the parallelogram
  h_tmp_2 <- h_tmp_2 %>%  left_join(frame_tmp, by = c("p_month" = "p_month", "DP_rev_i" = "time_w")) %>%
    filter(!is.na(weight)) %>%
    group_by(AP_i) %>%
    mutate(value_cum = cumsum(value)) %>%
    ungroup()

  hazard2 <- sum(h_tmp_2$value_cum*h_tmp_2$weight) / sum(h_tmp_2$weight)

  return(hazard2)
}









