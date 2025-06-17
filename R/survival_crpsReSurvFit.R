#' Survival continuously ranked probability score.
#'
#' Return the Survival Continuously Ranked Probability Score (SCRPS) of a \code{ReSurv} model.
#'
#' The model fit uses the theoretical framework of Hiabu et al. (2023), that relies on the
#'
#' @param ReSurvFit ReSurvFit object to use for the score computation.
#' @param user_data_set data.frame provided from the user to compute the survival CRPS, optional.
#'
#' @return Survival CRPS, \code{data.table} that contains the CRPS (\code{crps}) for each observation (\code{id}).
#'
#' @import data.table
#' @importFrom stats complete.cases
#'
#' @references
#' Pittarello, G., Hiabu, M., & Villegas, A. M. (2023). Chain Ladder Plus: a versatile approach for claims reserving. arXiv preprint arXiv:2301.03858.
#'
#' Therneau, T. M., & Lumley, T. (2015). Package ‘survival’. R Top Doc, 128(10), 28-33.
#'
#' Katzman, J. L., Shaham, U., Cloninger, A., Bates, J., Jiang, T., & Kluger, Y. (2018). DeepSurv: personalized treatment recommender system using a Cox proportional hazards deep neural network. BMC medical research methodology, 18(1), 1-12.
#'
#' Chen, T., He, T., Benesty, M., & Khotilovich, V. (2019). Package ‘xgboost’. R version, 90, 1-66.
#'
#' @export
survival_crps <- function(ReSurvFit,
                          user_hazard_frame=NULL,
                          user_data_set = NULL){

  UseMethod("survival_crps")

}


#' Survival continuously ranked probability score.
#'
#' Return the Survival Continuously Ranked Probability Score (SCRPS) of a \code{ReSurv} model.
#'
#' The model fit uses the theoretical framework of Hiabu et al. (2023), that relies on the
#'
#' @param ReSurvFit ReSurvFit object to use for the score computation.
#' @param user_data_set data.frame provided from the user to compute the survival CRPS, optional.
#'
#' @return Survival CRPS, \code{data.table} that contains the CRPS (\code{crps}) for each observation (\code{id}).
#'
#' @import data.table
#'
#' @references
#' Pittarello, G., Hiabu, M., & Villegas, A. M. (2023). Chain Ladder Plus: a versatile approach for claims reserving. arXiv preprint arXiv:2301.03858.
#'
#' Therneau, T. M., & Lumley, T. (2015). Package ‘survival’. R Top Doc, 128(10), 28-33.
#'
#' Katzman, J. L., Shaham, U., Cloninger, A., Bates, J., Jiang, T., & Kluger, Y. (2018). DeepSurv: personalized treatment recommender system using a Cox proportional hazards deep neural network. BMC medical research methodology, 18(1), 1-12.
#'
#' Chen, T., He, T., Benesty, M., & Khotilovich, V. (2019). Package ‘xgboost’. R version, 90, 1-66.
#'
#' @export
survival_crps.default <- function(ReSurvFit,
                                  user_hazard_frame=NULL,
                                  user_data_set = NULL){

  message('The object provided must be of class ReSurvFit')

}

#' Survival continuously ranked probability score.
#'
#' Return the Survival Continuously Ranked Probability Score (SCRPS) of a \code{ReSurv} model.
#'
#' The model fit uses the theoretical framework of Hiabu et al. (2023), that relies on the
#'
#' @param ReSurvFit ReSurvFit object to use for the score computation.
#' @param user_data_set data.frame provided from the user to compute the survival CRPS, optional.
#'
#' @return Survival CRPS, \code{data.table} that contains the CRPS (\code{crps}) for each observation (\code{id}).
#'
#' @import data.table
#'
#' @references
#' Pittarello, G., Hiabu, M., & Villegas, A. M. (2023). Chain Ladder Plus: a versatile approach for claims reserving. arXiv preprint arXiv:2301.03858.
#'
#' Therneau, T. M., & Lumley, T. (2015). Package ‘survival’. R Top Doc, 128(10), 28-33.
#'
#' Katzman, J. L., Shaham, U., Cloninger, A., Bates, J., Jiang, T., & Kluger, Y. (2018). DeepSurv: personalized treatment recommender system using a Cox proportional hazards deep neural network. BMC medical research methodology, 18(1), 1-12.
#'
#' Chen, T., He, T., Benesty, M., & Khotilovich, V. (2019). Package ‘xgboost’. R version, 90, 1-66.
#'
#' @export
survival_crps.ReSurvFit <- function(ReSurvFit,
                                    user_hazard_frame=NULL,
                                    user_data_set = NULL){



  if (is.null(user_hazard_frame)) {
    hazard_frame <- ReSurvFit$hazard_frame %>%
      select(-DP_i) %>%
      rename(dev_f_i = f_i, cum_dev_f_i = cum_f_i)
  } else{
    # hazard_frame<- user_hazard_frame%>%
    #   select(-DP_i) %>%
    #   rename(dev_f_i = f_i, cum_dev_f_i = cum_f_i)
    # browser()
    tmp <- IndividualDataPP(user_hazard_frame,
                            id=NULL,
                            categorical_features = ReSurvFit$IndividualDataPP$categorical_features,
                            continuous_features =   ReSurvFit$IndividualDataPP$continuous_features,
                            accident_period = ReSurvFit$IndividualDataPP$accident_period,
                            calendar_period = ReSurvFit$IndividualDataPP$calendar_period,
                            input_time_granularity = ReSurvFit$IndividualDataPP$input_time_granularity,
                            output_time_granularity = ReSurvFit$IndividualDataPP$output_time_granularity,
                            years = ReSurvFit$IndividualDataPP$years)


    columns_for_grouping <- unique(c(ReSurvFit$IndividualDataPP$categorical_features,ReSurvFit$IndividualDataPP$continuous_features,"AP_i"))

    out <- as.data.table(tmp$full.data)

    out <- out[DP_rev_i <= TR_i,.(.N),by=columns_for_grouping][,..columns_for_grouping]

    l4 <- list()

    l4$DP_rev_i <- min(tmp$training.data[,'DP_rev_i']):max(tmp$training.data[,'DP_rev_i'])

    l4<-do.call(CJ, c(l4, sorted = FALSE))

    newdata<-as.data.frame(setkey(out[,c(k=1,.SD)],k)[l4[,c(k=1,.SD)],allow.cartesian=TRUE][,k:=NULL])


    if(ReSurvFit$hazard_model=="COX"){

      Xc_tmp_bsln <- copy(out)

      for(cft in tmp$continuous_features){


        mnv <- min(tmp$training.data[cft])
        mxv <- max(tmp$training.data[cft])

        Xc_tmp_bsln[[cft]] <-2*(Xc_tmp_bsln[[cft]]-mnv)/(mxv-mnv)-1

      }


      X_tmp_bsln <- pkg.env$model.matrix.creator(data= out,
                                                 select_columns = tmp$categorical_features,
                                                 remove_first_dummy=T)


      X_tmp_bsln=cbind(X_tmp_bsln,Xc_tmp_bsln)



    newdata.bs <- pkg.env$df.2.fcst.nn.pp(data=tmp$training.data,
                                          newdata=newdata,
                                          continuous_features=tmp$continuous_features,
                                          categorical_features=tmp$categorical_features)


    Y =as.data.table(tmp$full.data[,c("DP_rev_i", "I", "TR_i")])
    Y <- Y[DP_rev_i <= TR_i,]


    model.out <- ReSurvFit$model.out$model.out
    coxlp <-  predict(model.out$cox,
                      newdata=newdata,
                      'lp')

    expg <- exp(coxlp)

    hazard_frame <- cbind(newdata, expg)
    colnames(hazard_frame)[dim(hazard_frame)[2]]="expg"


    max_dp_i =  pkg.env$maximum.time(tmp$years,tmp$input_time_granularity)

    unique_baseline <- ReSurvFit$hazard_frame %>%
      group_by(DP_i) %>%
      summarise(baseline = unique(baseline)) %>%
      mutate(DP_rev_i = max_dp_i - DP_i+1)


    hazard_frame <- hazard_frame %>%
      left_join(unique_baseline,by=c("DP_rev_i"))

    hazard_frame[,'hazard'] <- hazard_frame[,'baseline']*hazard_frame[,'expg']


    hazard_frame_updated <- pkg.env$hazard_data_frame(hazard=hazard_frame,
                                                      # Om.df=Om.df,
                                                      eta_old=.5,
                                                      categorical_features = tmp$categorical_features,
                                                      continuous_features = tmp$continuous_features,
                                                      calendar_period_extrapolation = tmp$calendar_period_extrapolation)


    hazard_frame <-  hazard_frame_updated %>%
      mutate(DP_i=max_dp_i-DP_rev_i+1) %>%
      relocate(DP_i, .after =  AP_i) %>%
      rename(f_i=dev_f_i,
             cum_f_i=cum_dev_f_i)


    }

    if(ReSurvFit$hazard_model=="XGB"){

      # browser()


      Y =as.data.table(tmp$full.data[,c("DP_rev_i", "I", "TR_i")])


      Xc <- copy(out)

      for(cft in tmp$continuous_features){


        mnv <- min(tmp$training.data[cft])
        mxv <- max(tmp$training.data[cft])

        Xc[[cft]] <-2*(Xc[[cft]]-mnv)/(mxv-mnv)-1

      }

      # Xc <- IndividualDataPP$training.data %>%
      #   reframe(across(all_of(IndividualDataPP$continuous_features),
      #                  scaler))



      # browser()
      if(!is.null(tmp$categorical_features)){

        X <- pkg.env$model.matrix.creator(data= out,
                                                   select_columns = tmp$categorical_features,
                                                   remove_first_dummy=T)

        # X <- pkg.env$model.matrix.creator(data= IndividualDataPP$training.data,
        #                                   select_columns = IndividualDataPP$categorical_features,
        #                                   remove_first_dummy=T)

        X=cbind(X,Xc)
      }else{


        X <- Xc


      }


      newdata.mx <- pkg.env$df.2.fcst.xgboost.pp(data=tmp$training.data,
                                                 newdata=newdata,
                                                 continuous_features=tmp$continuous_features,
                                                 categorical_features=tmp$categorical_features)


      pred <- predict(ReSurvFit$model.out$model.out,newdata.mx)

      expg <- exp(pred)

      hazard_frame <- cbind(newdata, expg)
      colnames(hazard_frame)[dim(hazard_frame)[2]]="expg"


      max_dp_i =  pkg.env$maximum.time(tmp$years,tmp$input_time_granularity)

      unique_baseline <- ReSurvFit$hazard_frame %>%
        group_by(DP_i) %>%
        summarise(baseline = unique(baseline)) %>%
        mutate(DP_rev_i = max_dp_i - DP_i+1)


      hazard_frame <- hazard_frame %>%
        left_join(unique_baseline,by=c("DP_rev_i"))

      hazard_frame[,'hazard'] <- hazard_frame[,'baseline']*hazard_frame[,'expg']

      hazard_frame_updated <- pkg.env$hazard_data_frame(hazard=hazard_frame,
                                                        # Om.df=Om.df,
                                                        eta_old=.5,
                                                        categorical_features = tmp$categorical_features,
                                                        continuous_features = tmp$continuous_features,
                                                        calendar_period_extrapolation = tmp$calendar_period_extrapolation)


      hazard_frame <-  hazard_frame_updated %>%
        mutate(DP_i=max_dp_i-DP_rev_i+1) %>%
        relocate(DP_i, .after =  AP_i) %>%
        rename(f_i=dev_f_i,
               cum_f_i=cum_dev_f_i)


    }

    if(ReSurvFit$hazard_model=="NN"){

      Y =as.data.table(tmp$full.data[,c("DP_rev_i", "I", "TR_i")])


      Xc <- copy(out)

      for(cft in tmp$continuous_features){


        mnv <- min(tmp$training.data[cft])
        mxv <- max(tmp$training.data[cft])

        Xc[[cft]] <-2*(Xc[[cft]]-mnv)/(mxv-mnv)-1

      }


      if(!is.null(tmp$categorical_features)){

        X <- pkg.env$model.matrix.creator(data= out,
                                          select_columns = tmp$categorical_features,
                                          remove_first_dummy=T)

        # X <- pkg.env$model.matrix.creator(data= IndividualDataPP$training.data,
        #                                   select_columns = IndividualDataPP$categorical_features,
        #                                   remove_first_dummy=T)

        X=cbind(X,Xc)
      }else{


        X <- Xc


      }
      # browser()

      newdata.mx <- pkg.env$df.2.fcst.nn.pp(data=tmp$training.data,
                                                 newdata=newdata,
                                                 continuous_features=tmp$continuous_features,
                                                 categorical_features=tmp$categorical_features)


      x_fc= reticulate::np_array(as.matrix(newdata.mx), dtype = "float32")

      beta_ams <- ReSurvFit$model.out$model.out$predict(input=x_fc)


      expg <- exp(beta_ams)

      hazard_frame <- cbind(newdata, expg)
      colnames(hazard_frame)[dim(hazard_frame)[2]]="expg"


      max_dp_i =  pkg.env$maximum.time(tmp$years,tmp$input_time_granularity)

      unique_baseline <- ReSurvFit$hazard_frame %>%
        group_by(DP_i) %>%
        summarise(baseline = unique(baseline)) %>%
        mutate(DP_rev_i = max_dp_i - DP_i+1)


      hazard_frame <- hazard_frame %>%
        left_join(unique_baseline,by=c("DP_rev_i"))

      hazard_frame[,'hazard'] <- hazard_frame[,'baseline']*hazard_frame[,'expg']

      hazard_frame_updated <- pkg.env$hazard_data_frame(hazard=hazard_frame,
                                                        # Om.df=Om.df,
                                                        eta_old=.5,
                                                        categorical_features = tmp$categorical_features,
                                                        continuous_features = tmp$continuous_features,
                                                        calendar_period_extrapolation = tmp$calendar_period_extrapolation)


      hazard_frame <-  hazard_frame_updated %>%
        mutate(DP_i=max_dp_i-DP_rev_i+1) %>%
        relocate(DP_i, .after =  AP_i) %>%
        rename(f_i=dev_f_i,
               cum_f_i=cum_dev_f_i)






    }



    hazard_frame <- data.table(hazard_frame)

    # Simplify the code and save useful attributes
    categorical_features <- ReSurvFit$IndividualDataPP$categorical_features
    continuous_features <- ReSurvFit$IndividualDataPP$continuous_features
    max_dp_i =  pkg.env$maximum.time(ReSurvFit$IndividualDataPP$years,
                                     ReSurvFit$IndividualDataPP$input_time_granularity)
    conversion_factor =ReSurvFit$IndividualDataPP$conversion_factor
    calendar_period_extrapolation = ReSurvFit$IndividualDataPP$calendar_period_extrapolation

    # find groups
    hazard_frame <- hazard_frame[,
                                 ix_group:=.GRP,
                                 by=c(categorical_features,
                                      continuous_features)]

    # find unique combinations of groups
    tmp_unq_grp <- data.table(unique(data.frame(hazard_frame)[c(categorical_features,
                                                                continuous_features,
                                                                'ix_group')]))

    # find the test set

    test_for_crps = ReSurvFit$IndividualDataPP$full.data %>%
      filter(DP_rev_i <= TR_i)

    # Save the different curves

    hazard_list<- split(hazard_frame, hazard_frame$ix_group)

    hazard_list<-lapply(hazard_list, function(x) x[order(DP_rev_i),
                                                   .(DP_rev_i,
                                                     cdf2_i = (1-S_i)^2,
                                                     S2_i=S_i^2,
                                                     x.vals=diff(c(0,DP_rev_i)))])


    test_for_crps = merge(test_for_crps,
                          tmp_unq_grp,
                          by=c(categorical_features, continuous_features),
                          all.x=T,
                          all.y=F)

    test_for_crps[['id']] = 1:dim(test_for_crps)[1]

    test_for_crps<-test_for_crps[complete.cases(test_for_crps),]

    setDT(test_for_crps)

    tmp <- test_for_crps[,.(crps=survival_information(x=DP_rev_i,
                                                      group=ix_group,
                                                      hazard_list = hazard_list)),
                         by=id]


    return(tmp)


  }





  # Elaborate features on the test set
  if(is.null(user_data_set)){


    # Simplify the code and save useful attributes
    categorical_features <- ReSurvFit$IndividualDataPP$categorical_features
    continuous_features <- ReSurvFit$IndividualDataPP$continuous_features
    max_dp_i =  pkg.env$maximum.time(ReSurvFit$IndividualDataPP$years,
                                     ReSurvFit$IndividualDataPP$input_time_granularity)
    conversion_factor =ReSurvFit$IndividualDataPP$conversion_factor
    calendar_period_extrapolation = ReSurvFit$IndividualDataPP$calendar_period_extrapolation

    test_for_crps=ReSurvFit$IndividualDataPP$full.data %>%
      filter(DP_rev_i <= TR_i)%>%
    mutate(
      DP_rev_o = floor(max_dp_i*conversion_factor)-ceiling(DP_i*conversion_factor+((AP_i-1)%%(1/conversion_factor))*conversion_factor) +1,
      AP_o = ceiling(AP_i*conversion_factor)
    ) %>%
    mutate(TR_o= AP_o-1) %>%
    mutate(across(all_of(categorical_features),
                  as.factor)) %>%
    select(all_of(categorical_features),
           all_of(continuous_features),
           all_of(switch(calendar_period_extrapolation, 'RP_i', NULL)),
           AP_i,
           AP_o,
           DP_i,
           DP_rev_i,
           DP_rev_o,
           TR_i,
           TR_o,
           I) %>%
    as.data.table()

    # in case the test set is not available we do not compute it and return NULL
    if(dim(test_for_crps)[1]==0){
      warning('No test set available in your data to compute the Survival CRPS')
      return(NULL)
    }

    # find groups
    setDT(hazard_frame)
    hazard_frame <- hazard_frame[,
                                 ix_group:=.GRP,
                                 by=c(categorical_features,
                                      continuous_features)]

    tmp_unq_grp <- data.table(unique(data.frame(hazard_frame)[c(categorical_features,
                                                                continuous_features,
                                                                'ix_group')]))



    }else{


      tmp_cond= colnames(ReSurvFit$IndividualDataPP$starting.data) %in% colnames(user_data_set)

      tmp_training_set = as.data.table(ReSurvFit$IndividualDataPP$starting.data)[,..tmp_cond]

      # Simple rbind (full starting data and new data to compute CRPS)
      tmp_fdata = rbind(tmp_training_set,
                        user_data_set)



      conversion_factor= ReSurvFit$IndividualDataPP$conversion_factor
      continuous_features=ReSurvFit$IndividualDataPP$continuous_features
      categorical_features=ReSurvFit$IndividualDataPP$categorical_features

      # find groups
      setDT(hazard_frame)
      hazard_frame <- hazard_frame[,
                                   ix_group:=.GRP,
                                   by=c(categorical_features,
                                        continuous_features)]

      # find unique combinations of groups
      tmp_unq_grp <- data.table(unique(data.frame(hazard_frame)[c(categorical_features,
                                                                  continuous_features,
                                                                  'ix_group')]))


      # Process the data
      tmp_idata = IndividualDataPP(tmp_fdata,
                                 continuous_features=continuous_features,
                                 categorical_features=categorical_features,
                                 accident_period=ReSurvFit$IndividualDataPP$accident_period,
                                 calendar_period=ReSurvFit$IndividualDataPP$calendar_period,
                                 input_time_granularity=ReSurvFit$IndividualDataPP$input_time_granularity,
                                 output_time_granularity=ReSurvFit$IndividualDataPP$output_time_granularity,
                                 years=ReSurvFit$IndividualDataPP$years,
                                 calendar_period_extrapolation=ReSurvFit$IndividualDataPP$calendar_period_extrapolation,
                                 continuous_features_spline=NULL)

      test_for_crps = tmp_idata$full.data %>%
        filter(DP_rev_i <= TR_i)

      max_dp_i =  pkg.env$maximum.time(ReSurvFit$IndividualDataPP$years,ReSurvFit$IndividualDataPP$input_time_granularity)


      calendar_period_extrapolation=ReSurvFit$IndividualDataPP$calendar_period_extrapolation

      test_for_crps=test_for_crps%>%
        mutate(
          DP_rev_o = floor(max_dp_i*conversion_factor)-ceiling(DP_i*conversion_factor+((AP_i-1)%%(1/conversion_factor))*conversion_factor) +1,
          AP_o = ceiling(AP_i*conversion_factor)
        ) %>%
        mutate(TR_o= AP_o-1) %>%
        mutate(across(all_of(categorical_features),
                      as.factor)) %>%
        select(all_of(categorical_features),
               all_of(continuous_features),
               all_of(switch(calendar_period_extrapolation, 'RP_i', NULL)),
               AP_i,
               AP_o,
               DP_i,
               DP_rev_i,
               DP_rev_o,
               TR_i,
               TR_o,
               I) %>%
        as.data.table()


    }


  # Save the different curves

  hazard_list<- split(hazard_frame, hazard_frame$ix_group)

  hazard_list<-lapply(hazard_list, function(x) x[order(DP_rev_i),
                                                 .(DP_rev_i,
                                                   cdf2_i = (1-S_i)^2,
                                                   S2_i=S_i^2,
                                                   x.vals=diff(c(0,DP_rev_i)))])


  test_for_crps = merge(test_for_crps,
                        tmp_unq_grp,
                        by=c(categorical_features, continuous_features),
                        all.x=T,
                        all.y=F)

  test_for_crps[['id']] = 1:dim(test_for_crps)[1]

  test_for_crps<-test_for_crps[complete.cases(test_for_crps),]

  tmp <- test_for_crps[,.(crps=survival_information(x=DP_rev_i,
                                                   group=ix_group,
                                                   hazard_list = hazard_list)),
                       by=id]


  return(tmp)

}





