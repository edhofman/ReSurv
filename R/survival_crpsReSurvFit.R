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
survival_crps <- function(ReSurvFit,
                          user_data_set= NULL){

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
                                  user_data_set= NULL){

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
                                    user_data_set= NULL){


  hazard_frame <- ReSurvFit$hazard_frame

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



  # Elaborate features on the test set
  if(is.null(user_data_set)){

    # in case the test set is not available we do not compute it and return NULL
    if(dim(test_for_crps)[1]==0){
      warning('No test set available in your data to compute the Survival CRPS')
      return(NULL)
    }

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
    as.data.table()}else{

      tmp_cond= colnames(ReSurvFit$IndividualDataPP$starting.data) %in% colnames(user_data_set)

      tmp_training_set = ReSurvFit$IndividualDataPP$starting.data[,tmp_cond]

      # Simple rbind (full starting data and new data to compute CRPS)
      tmp_fdata = rbind(tmp_training_set,
                        user_data_set)

      # Process the data
      tmp_idata = IndividualDataPP(tmp_fdata,
                                 continuous_features=continuous_features,
                                 categorical_features=categorical_features,
                                 accident_period=ReSurvFit$IndividualDataPP$accident_period,
                                 calendar_period="CM",
                                 input_time_granularity=ReSurvFit$IndividualDataPP$input_time_granularity,
                                 output_time_granularity=ReSurvFit$IndividualDataPP$output_time_granularity,
                                 years=ReSurvFit$IndividualDataPP$years,
                                 calendar_period_extrapolation=ReSurvFit$IndividualDataPP$calendar_period_extrapolation,
                                 continuous_features_spline=NULL)

      test_for_crps = tmp_idata$full.data %>%
        filter(DP_rev_i <= TR_i)

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


  # browser()
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







