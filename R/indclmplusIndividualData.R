#' Fit individual chain ladder plus models.
#'
#' This function fits and computes the reserves for the individual clmplus models.
#'
#' @param model: model to be chosen for the fitting.
#' @param tie: ties handling, default is efron approach.
#' @param baseline: handling the baseline hazard. Default is a spline.
#'
#' @return individual chain ladder plus fit.
#' @export
ReSurv <- function(IndividualData,
                   hazard.model="cox",
                   tie='efron',
                   baseline="spline",
                   percentage_data_training=.8,
                   hparameters=list(),
                   seed=1
){

  UseMethod("ReSurv")

}

#' Fit individual chain ladder plus models.
#'
#' This function fits and computes the reserves for the individual clmplus models.
#'
#' @param model: model to be chosen for the fitting.
#' @param tie: ties handling. Default is efron approach.
#' @param baseline: handling the baseline hazard. Default is a spline.
#'
#' @return individual chain ladder plus fit.
#' @export
ReSurv.default <- function(IndividualData,
                               hazard.model="cox",
                               tie='efron',
                               baseline="spline",
                           hparameters=list(),
                           seed=1){

  message('The object provided must be of class IndividualData')

}



#' Fit chain-ladder+ to reverse time triangles.
#'
#' This function fits and computes the reserves for the invidual clmplus models.
#'
#' @param model: model to be chosen for the fitting.
#' @param tie: ties handling, default is efron approach.
#' @param baseline: handling the baseline hazard. Default is a spline.
#'
#' @return individual chain ladder plus fit.
#' @export
ReSurv.IndividualData <- function(IndividualData,
                               hazard.model="cox",
                               tie='efron',
                               baseline="spline",
                               continuous.features.scaling.method="minmax",
                               random_seed=1,
                               hparameters=list(),
                               percentage_data_training=.8
                               ){


  set.seed(random_seed)

  formula_ct <- as.formula(IndividualData$string_formula_i)

  X_i <- pkg.env$model.matrix.creator(data= IndividualData$training.data,
                                        select_columns = c('AP_i',IndividualData$categorical_features))

  hz_names_i = pkg.env$model.matrix.extract.hazard.names(X=X_i,
                                                         string_formula=IndividualData$string_formula_i,
                                                         data=IndividualData$training.data)

  ##################################################################################
  # By now I created a separate function that runs the code. We need to think if there is a smart
  # way of creating only one matrix. We do the same procudure twice.
  # Create X matrix

  X <- pkg.env$model.matrix.creator(data= IndividualData$training.data,
                            select_columns = IndividualData$categorical_features)

  scaler <- pkg.env$scaler(continuous.features.scaling.method=continuous.features.scaling.method)
  Xc <- IndividualData$training.data %>%
    summarise(across(all_of(IndividualData$continuous_features),
                     scaler))

  ##################################################################################
  # Here we place the differen fitting routines

  if(hazard.model=="cox"){model.out <- pkg.env$fit_cox_model(data=IndividualData$training,
                                                     formula_ct,
                                                     X,
                                                     X_i)}

  if(hazard.model=="deepsurv"){

    training_test_split = pkg.env$check.traintestsplit(percentage_data_training)

    datads_pp = pkg.env$deep_surv_pp(X=cbind(X,Xc),
                           Y=IndividualData$training.data[,c("DP_rev_i", "I", "TR_i")],
                           training_test_split = training_test_split)



    model.out <- pkg.env$fit_deep_surv(datads_pp,
                                       hparameters=hparameters)

    return(model.out)


    }

  ##################################################################################

  ##################################################################################
  # The following steps are data specific.
  # They need to be generalized.

  baseline_out <- pkg.env$hazard_baseline_model(data=IndividualData$training.data,
                                        cox=model.out$cox,
                                        hazard=NULL,
                                        baseline=baseline,
                                        conversion_factor=IndividualData$conversion_factor,
                                        nk=50,
                                        nbin=48,
                                        phi=1)

  ##################################################################################

  l = length(model.out$beta_ams)

  #need placeholder for latest i mirror cl behaviour
  tmp= sapply(1:l,function(x) c(baseline_out$bs_hazard$hazard)*exp(model.out$beta_ams[x]))

  hazard_cl <- (sapply(seq_along(hz_names_i$time),
                       pkg.env$hazard_f,
                       enter= hz_names_i$enter,
                       time=hz_names_i$time,
                       exit=hz_names_i$exit,
                       event=hz_names_i$event))

  hazard = as.matrix(cbind(hazard_cl,
                           tmp ))


  max_DP <- max(IndividualData$training$DP_rev_o)

  ############################################################
  #check
  hazard_q <- matrix(nrow=max_DP, ncol=(ncol(hazard)-1)*IndividualData$conversion_factor)
  ############################################################


  #group to quarters, this is relatively time consuming
  for( i in 1:max_DP){
    frame_tmp <- data.table(IndividualData$training) %>% filter(TR_o<i) %>%
      filter(DP_rev_o>=i) %>%
      mutate(time_w = round(ifelse(DP_rev_o==i, DP_rev_i, AP_i-1+3*(i-AP_o+1) ),10) ) %>%
      mutate(weight = ifelse(DP_rev_o==i, (DP_rev_i-1)%%3+1, 3 )) %>%
      mutate(DP_o_tmp = i) %>%
      mutate(p_month = (AP_i-1)%%3+1) %>%
      group_by(p_month, time_w) %>%
      summarize(weight=sum(weight*I) )

    frame_tmp2 <- frame_tmp %>%  reshape2::dcast(time_w ~p_month, value.var="weight")

    hazard_q[i,] <- mapply(pkg.env$m_to_q_hazard,
                           1:ncol(hazard_q),
                           MoreArgs=list(hazard=hazard,
                           frame_tmp=frame_tmp,
                           frame_tmp2=frame_tmp2))}

  #create monthly development factors, though not outputted later on.
  colnames(hazard) <- c("CL",hz_names_i$names_hazard )

  df_monthly <- (2+hazard)/(2-hazard)

  #Remove last row,since doesn't make sense
  df_monthly <- as.data.frame(df_monthly[1:(nrow(df_monthly)-1),]) %>%
    map_df(rev) %>%
    mutate(DM=row_number())

  X_o <- pkg.env$model.matrix.creator(data= IndividualData$training.data,
                                        select_columns = c('AP_o',IndividualData$categorical_features))

  hz_names_o = pkg.env$model.matrix.extract.hazard.names(X=X_o,
                                                         string_formula=IndividualData$string_formula_o,
                                                         data=IndividualData$training.data)

hazard_cl <- (sapply(seq_along(hz_names_o$time),
                     pkg.env$hazard_f,
                     enter= hz_names_o$enter,
                     time= hz_names_o$time,
                     exit= hz_names_o$exit,
                     event= hz_names_o$event))

q_hazard = as.matrix(cbind(hazard_cl, hazard_q )   )

colnames(q_hazard) <- c("CL",hz_names_o$names_hazard )


df_quarterly <- (2+q_hazard)/(2-q_hazard)


df_quarterly <- as.data.frame(df_quarterly[1:(nrow(df_quarterly)-1),]) %>%
  map_df(rev) %>%
  mutate(DP_i=row_number())

out=list(df = df_quarterly,
         hazard= q_hazard,
         IndividualData=IndividualData)

class(out) <- c('ReSurvFit')

out

}




