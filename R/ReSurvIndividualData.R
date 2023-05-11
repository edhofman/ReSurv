#' Fit \code{ReSurv} models on the individual data.
#'
#' This function fits and computes the reserves for the \code{ReSurv} models
#'
#' The model fit uses the theoretical framework of Hiabu et al. (2023), that relies on the
#' correspondence between hazard models and development factors:
#'
#' To be completed with final notation of the paper.
#'
#' The \code{ReSurv} package assumes proportional hazard models.
#' Given an i.i.d. sample \eqn{\left\{y_i,x_i\right\}_{i=1, \ldots, n}} the individual hazard at time \eqn{t} is:
#'
#' \eqn{\lambda_i(t)=\lambda_0(t)e^{y_i(x_i)}}
#'
#' Composed of a baseline \eqn{\lambda_0(t)} and a proportional effect \eqn{e^{y_i(x_i)}}.
#'
#' Currently, the implementation allows to optimize the partial likelihood (concerning the proportional effects) using one of the following statistical learning approaches:
#' \itemize{
#' \item{\href{https://cran.r-project.org/web/packages/survival/vignettes/survival.pdf}{Cox}}
#' \item{\href{https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-018-0482-1}{Neural Networks}}
#' \item{\href{https://xgboost.readthedocs.io/en/stable/}{eXtreme Gradient Boosting}}
#' \item{\href{https://cran.r-project.org/web/packages/LTRCtrees/LTRCtrees.pdf}{Left-truncated Right-Censored Trees}}
#' }
#'
#'
#' @param IndividualData IndividualData object to use for the \code{ReSurv} fit.
#' @param hazard_model \code{character}, hazard model supported from our package, must be provided as a string. The model can be chosen from:
#' \itemize{
#' \item{\code{"cox"}: Standard Cox model for the hazard.}
#' \item{\code{"deep_surv"}: Deep Survival Neural Network.}
#' \item{\code{"xgboost"}: eXtreme Gradient Boosting.}
#' }
#' @param tie ties handling, default is the Efron approach.
#' @param baseline handling the baseline hazard. Default is a spline.
#' @param continuous_features_scaling_method method to preprocess the features
#' @param random_seed \code{integer}, random seed set for reproducibility
#' @param hparameters \code{list}, hyperparameters for the machine learning models. It will be disregarded for the cox approach.
#' @param percentage_data_training \code{numeric}, percentage of data used for training on the upper triangle.
#' @param grouping_method \code{character}, use probability or exposure approach to group from input to output development factors. Choice between:
#' \itemize{
#' \item{\code{"exposure"}}
#' \item{\code{"probability"}}
#' }
#' Default is \code{"exposure"}.
#' @param check_value \code{numeric}, check hazard value on initial granularity, if above threshold we increase granularity to try and adjust the development factor.
#'
#'
#' @return ReSurv fit.
#'
#' @import reticulate
#' @import tidyverse
#' @import xgboost
#' @import rpart
#' @import LTRCtrees
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
ReSurv <- function(IndividualData,
                   hazard_model="cox",
                   tie='efron',
                   baseline="spline",
                   continuous_features_scaling_method="minmax",
                   random_seed=1964,
                   hparameters=list(),
                   percentage_data_training=.8,
                   grouping_method = "exposure",
                   check_value = 1.85
){

  UseMethod("ReSurv")

}
#' Fit \code{ReSurv} models on the individual data.
#'
#' This function fits and computes the reserves for the \code{ReSurv} models
#'
#' The model fit uses the theoretical framework of Hiabu et al. (2023), that relies on the
#' correspondence between hazard models and development factors:
#'
#' To be completed with final notation of the paper.
#'
#' The \code{ReSurv} package assumes proportional hazard models.
#' Given an i.i.d. sample \eqn{\left\{y_i,x_i\right\}_{i=1, \ldots, n}} the individual hazard at time \eqn{t} is:
#'
#' \eqn{\lambda_i(t)=\lambda_0(t)e^{y_i(x_i)}}
#'
#' Composed of a baseline \eqn{\lambda_0(t)} and a proportional effect \eqn{e^{y_i(x_i)}}.
#'
#' Currently, the implementation allows to optimize the partial likelihood (concerning the proportional effects) using one of the following statistical learning approaches:
#' \itemize{
#' \item{\href{https://cran.r-project.org/web/packages/survival/vignettes/survival.pdf}{Cox}}
#' \item{\href{https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-018-0482-1}{Neural Networks}}
#' \item{\href{https://xgboost.readthedocs.io/en/stable/}{eXtreme Gradient Boosting}}
#' \item{\href{https://cran.r-project.org/web/packages/LTRCtrees/LTRCtrees.pdf}{Left-truncated Right-Censored Trees}}
#' }
#'
#'
#' @param IndividualData IndividualData object to use for the \code{ReSurv} fit.
#' @param hazard_model \code{character}, hazard model supported from our package, must be provided as a string. The model can be chosen from:
#' \itemize{
#' \item{\code{"cox"}: Standard Cox model for the hazard.}
#' \item{\code{"deep_surv"}: Deep Survival Neural Network.}
#' \item{\code{"xgboost"}: eXtreme Gradient Boosting.}
#' }
#' @param tie ties handling, default is the Efron approach.
#' @param baseline handling the baseline hazard. Default is a spline.
#' @param continuous_features_scaling_method method to preprocess the features
#' @param random_seed \code{integer}, random seed set for reproducibility
#' @param hparameters \code{list}, hyperparameters for the machine learning models. It will be disregarded for the cox approach.
#' @param percentage_data_training \code{numeric}, percentage of data used for training on the upper triangle.
#' @param grouping_method \code{character}, use probability or exposure approach to group from input to output development factors. Choice between:
#' \itemize{
#' \item{\code{"exposure"}}
#' \item{\code{"probability"}}
#' }
#' Default is \code{"exposure"}.
#' @param check_value \code{numeric}, check hazard value on initial granularity, if above threshold we increase granularity to try and adjust the development factor.
#'
#'
#' @return ReSurv fit.
#'
#' @import reticulate
#' @import tidyverse
#' @import xgboost
#' @import rpart
#'
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
ReSurv.default <- function(IndividualData,
                           hazard_model="cox",
                           tie='efron',
                           baseline="spline",
                           continuous_features_scaling_method="minmax",
                           random_seed=1,
                           hparameters=list(),
                           percentage_data_training=.8,
                           grouping_method = "exposure",
                           check_value = 1.85){

  message('The object provided must be of class IndividualData')

}



#' Fit \code{ReSurv} models on the individual data.
#'
#' This function fits and computes the reserves for the \code{ReSurv} models
#'
#' The model fit uses the theoretical framework of Hiabu et al. (2023), that relies on the
#' correspondence between hazard models and development factors:
#'
#' To be completed with final notation of the paper.
#'
#' The \code{ReSurv} package assumes proportional hazard models.
#' Given an i.i.d. sample \eqn{\left\{y_i,x_i\right\}_{i=1, \ldots, n}} the individual hazard at time \eqn{t} is:
#'
#' \eqn{\lambda_i(t)=\lambda_0(t)e^{y_i(x_i)}}
#'
#' Composed of a baseline \eqn{\lambda_0(t)} and a proportional effect \eqn{e^{y_i(x_i)}}.
#'
#' Currently, the implementation allows to optimize the partial likelihood (concerning the proportional effects) using one of the following statistical learning approaches:
#' \itemize{
#' \item{\href{https://cran.r-project.org/web/packages/survival/vignettes/survival.pdf}{Cox}}
#' \item{\href{https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-018-0482-1}{Neural Networks}}
#' \item{\href{https://xgboost.readthedocs.io/en/stable/}{eXtreme Gradient Boosting}}
#' \item{\href{https://cran.r-project.org/web/packages/LTRCtrees/LTRCtrees.pdf}{Left-truncated Right-Censored Trees}}
#' }
#'
#'
#' @param IndividualData IndividualData object to use for the \code{ReSurv} fit.
#' @param hazard_model \code{character}, hazard model supported from our package, must be provided as a string. The model can be chosen from:
#' \itemize{
#' \item{\code{"cox"}: Standard Cox model for the hazard.}
#' \item{\code{"deep_surv"}: Deep Survival Neural Network.}
#' \item{\code{"xgboost"}: eXtreme Gradient Boosting.}
#' }
#' @param tie ties handling, default is the Efron approach.
#' @param baseline handling the baseline hazard. Default is a spline.
#' @param continuous_features_scaling_method method to preprocess the features
#' @param random_seed \code{integer}, random seed set for reproducibility
#' @param hparameters \code{list}, hyperparameters for the machine learning models. It will be disregarded for the cox approach.
#' @param percentage_data_training \code{numeric}, percentage of data used for training on the upper triangle.
#' @param grouping_method \code{character}, use probability or exposure approach to group from input to output development factors. Choice between:
#' \itemize{
#' \item{\code{"exposure"}}
#' \item{\code{"probability"}}
#' }
#' Default is \code{"exposure"}.
#' @param check_value \code{numeric}, check hazard value on initial granularity, if above threshold we increase granularity to try and adjust the development factor.
#'
#' @return ReSurv fit.
#'
#' @import reticulate
#' @import tidyverse
#' @import xgboost
#' @import rpart
#'
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
ReSurv.IndividualData <- function(IndividualData,
                               hazard_model="cox",
                               tie='efron',
                               baseline="spline",
                               continuous_features_scaling_method="minmax",
                               random_seed=1,
                               hparameters=list(),
                               percentage_data_training=.8,
                               grouping_method = "exposure",
                               check_value = 1.85
                               ){


  set.seed(random_seed)

  formula_ct <- as.formula(IndividualData$string_formula_i)

  newdata <- create.df.2.fcst(IndividualData=IndividualData,
                                      hazard_model=hazard_model)


  # create data frame of occurrencies to weight development factors
  Om.df <-   pkg.env$create.om.df(training.data=IndividualData$training.data,
                         input_time_granularity=IndividualData$input_time_granularity,
                         years=IndividualData$years)

  if(hazard_model=="cox"){

    data=IndividualData$training.data
    model.out <- pkg.env$fit_cox_model(data=data,
                                       formula_ct=formula_ct,
                                       newdata=newdata)

    # tmp <- pkg.env$spline_hp(hparameters,IndividualData)


    # data <- IndividualData$training.data

    bs_hazard <- basehaz( model.out$cox, centered=FALSE) %>%
      mutate(hazard = hazard-lag(hazard,default=0))

    # baseline_out <- pkg.env$hazard_baseline_model(data=IndividualData$training.data,
    #                                               cox=cox.model,
    #                                               hazard=NULL,
    #                                               baseline=baseline,
    #                                               conversion_factor=IndividualData$conversion_factor,
    #                                               nk=tmp$nk,
    #                                               nbin=tmp$nbin,
    #                                               phi=tmp$phi)


    bsln <- data.frame(baseline=bs_hazard$hazard,
                       DP_rev_i=ceiling(bs_hazard$time))  #$hazard

    hazard_frame <- cbind(newdata, model.out$expg)
    colnames(hazard_frame)[dim(hazard_frame)[2]]="expg"



    }

  if(hazard_model=="deepsurv"){

    X <- pkg.env$model.matrix.creator(data= IndividualData$training.data,
                                      select_columns = IndividualData$categorical_features)

    scaler <- pkg.env$scaler(continuous_features_scaling_method=continuous_features_scaling_method)

    Xc <- IndividualData$training.data %>%
      summarize(across(all_of(IndividualData$continuous_features),
                       scaler))

    training_test_split = pkg.env$check.traintestsplit(percentage_data_training)

    X = cbind(X,Xc)

    Y=IndividualData$training.data[,c("DP_rev_i", "I", "TR_i")]

    datads_pp = pkg.env$deep_surv_pp(X=X,
                           Y=Y,
                           training_test_split = training_test_split)

    hparameters <- pkg.env$nn_hparameter_nodes_grid(hparameters)

    hparameters <- list(params=as.list.data.frame(hparameters),
                        verbose=hparameters$verbose,
                        epochs = hparameters$epochs,
                        num_workers = hparameters$num_workers)


    model.out <- pkg.env$fit_deep_surv(datads_pp,
                                       params=hparameters$params,
                                       verbose = hparameters$verbose,
                                       epochs = hparameters$epochs,
                                       num_workers = hparameters$num_workers,
                                       seed = random_seed)

    # bsln <- model.out$compute_baseline_hazards(
    #   input = datads_pp$x_train,
    #   target = datads_pp$y_train,
    #   batch_size = hparameters$batch_size)

    bsln <- pkg.env$baseline.calc(hazard_model = hazard_model,
                                  model.out = model.out,
                                  X=X,
                                  Y=Y)

    newdata.mx <- pkg.env$df.2.fcst.nn.pp(data=IndividualData$training.data,
                                          newdata=newdata,
                                          continuous_features=IndividualData$continuous_features,
                                          categorical_features=IndividualData$categorical_features)

    x_fc= reticulate::np_array(as.matrix(newdata.mx), dtype = "float32")

    beta_ams <- model.out$predict(input=x_fc,
                                  batch_size=hparameters$batch_size,
                                  num_workers=hparameters$num_workers)
    #make to hazard relative to initial model, to have similiar interpretation as standard cox

    benchmark_id <- pkg.env$benchmark_id(X = X,
                                           Y =Y ,
                                           newdata.mx = newdata.mx
    )

    pred_relative <- beta_ams - beta_ams[benchmark_id]

    expg <- exp(pred_relative)

    hazard_frame <- cbind(newdata,expg)
    bsln <- data.frame(baseline=bsln,
                       DP_rev_i=sort(as.integer(unique(IndividualData$training.data$DP_rev_i))))

  }

  if(hazard_model == "xgboost"){

    X <- pkg.env$model.matrix.creator(data= IndividualData$training.data,
                                  select_columns = IndividualData$categorical_features,
                                  remove_first_dummy=T)

    scaler <- pkg.env$scaler(continuous_features_scaling_method = continuous_features_scaling_method)

    Xc <- IndividualData$training.data %>%
      summarize(across(all_of(IndividualData$continuous_features),
                       scaler))

    training_test_split = pkg.env$check.traintestsplit(percentage_data_training)

    X=cbind(X,Xc)

    Y=IndividualData$training.data[,c("DP_rev_i", "I", "TR_i")]

    datads_pp <- pkg.env$xgboost_pp(X=X,
                                    Y=Y,
                                    training_test_split=training_test_split)

    model.out <- pkg.env$fit_xgboost(datads_pp,
                                     hparameters=hparameters)
    return(list(data=X,
                model.out = model.out))

    bsln <- pkg.env$baseline.calc(hazard_model = hazard_model,
                                  model.out = model.out,
                                  X=X,
                                  Y=Y)

    newdata.mx <- pkg.env$df.2.fcst.xgboost.pp(data=IndividualData$training.data,
                                          newdata=newdata,
                                          continuous_features=IndividualData$continuous_features,
                                          categorical_features=IndividualData$categorical_features)

    pred <- predict(model.out,newdata.mx)

    #make to hazard relative to initial model, to have similiar interpretation as standard cox
    newdata.bs <- pkg.env$df.2.fcst.nn.pp(data=IndividualData$training.data,
                                          newdata=newdata,
                                          continuous_features=IndividualData$continuous_features,
                                          categorical_features=IndividualData$categorical_features)
    benchmark_id <- pkg.env$benchmark_id(X = X,
                                         Y =Y ,
                                         newdata.mx = newdata.bs,
                                         remove_first_dummy=T)


    pred_relative <- pred - pred[benchmark_id]

    expg <- exp(pred_relative)

    hazard_frame <- cbind(newdata,expg)

    bsln <- data.frame(baseline=bsln,
                       DP_rev_i=sort(as.integer(unique(IndividualData$training.data$DP_rev_i))))

  }

  if(hazard_model == "LTRCtrees"){

    X <- pkg.env$model.matrix.creator(data= IndividualData$training.data,
                                      select_columns = IndividualData$categorical_features,
                                      remove_first_dummy=T)

    scaler <- pkg.env$scaler(continuous_features_scaling_method = continuous_features_scaling_method)

    Xc <- IndividualData$training.data %>%
      summarize(across(all_of(IndividualData$continuous_features),
                       scaler))

    training_test_split = pkg.env$check.traintestsplit(percentage_data_training)

    X=cbind(X,Xc)

    Y=IndividualData$training.data[,c("DP_rev_i", "I", "TR_i")]

    control.pars <- do.call(rpart.control, hparameters)

    model.out <- pkg.env$fit_LTRCtrees(data=IndividualData$training.data,
                                       formula_ct=formula_ct,
                                       newdata=newdata,
                                       control.pars)




    bsln <- pkg.env$baseline.calc(hazard_model = hazard_model,
                                  model.out = model.out$cox,
                                  X=X,
                                  Y=Y,
                                  training_df=IndividualData$training.data)


    #make to hazard relative to initial model, to have similiar interpretation as standard cox


    #make to hazard relative to initial model, to have similiar interpretation as standard cox
    # benchmark_id <- pkg.env$benchmark_id(X = X,
    #                                      Y =Y ,
    #                                      newdata.mx = newdata)

    pred <- predict(model.out$cox,newdata)

    benchmark_id <- 1
    # pred_relative <- model.out$expg/model.out$expg[benchmark_id]

    pred_relative <- exp(pred - pred[benchmark_id])

    # exp(pred_relative)

    hazard_frame <- cbind(newdata, expg=pred_relative)
    # colnames(hazard_frame)[dim(hazard_frame)[2]]="expg"

    bsln <- data.frame(baseline=bsln,
                       DP_rev_i=sort(as.integer(unique(IndividualData$training.data$DP_rev_i))))



  }

  ##################################################################################


  hazard_frame <- hazard_frame %>%
    full_join(bsln,
              by="DP_rev_i") %>%
    as.data.frame() %>%
    replace_na(list(baseline=0))

  hazard_frame[,'hazard'] <- hazard_frame[,'baseline']*hazard_frame[,'expg']



  #return(hazard_frame)

  #need placeholder for latest i mirror cl behaviour

  # hazard_cl <- (sapply(seq_along(hz_names_i$time),
  #                      pkg.env$hazard_f,
  #                      enter= hz_names_i$enter,
  #                      time=hz_names_i$time,
  #                      exit=hz_names_i$exit,
  #                      event=hz_names_i$event))




  ############################################################
  #check

  #hazard_q <- matrix(nrow=max_DP, ncol=(ncol(hazard)-1)*IndividualData$conversion_factor)
  #eta_o <- c()


  ############################################################

  #Add development and relevant survival values to the hazard_frame
  hazard_frame_updated <- pkg.env$hazard_data_frame(hazard=hazard_frame,
                                                    Om.df=Om.df,
                                                    categorical_features = IndividualData$categorical_features,
                                                    continuous_features = IndividualData$continuous_features,
                                                    calendar_period_extrapolation = IndividualData$calendar_period_extrapolation)


  hazard_frame_grouped <- pkg.env$covariate_mapping(
    hazard_frame = hazard_frame_updated,
    categorical_features = IndividualData$categorical_features,
    continuous_features = IndividualData$continuous_features,
    conversion_factor = IndividualData$conversion_factor,
    calendar_period_extrapolation = IndividualData$calendar_period_extrapolation
    )

  # I add missing observations starting from ALL the combinations of AP_i and DP_i
  # in the training data for EVERY combination of features.
  # In case some development (accident) periods are missing I fill the holes
  # (e.g. for a given I do not observe DP 11, 12 in the middle of the triangle but my DP go from 1 to 25 and I add them)


  missing.obsevations <- pkg.env$fill_data_frame(data=IndividualData$full.data,
                                                 continuous_features=IndividualData$continuous_features,
                                                 categorical_features=IndividualData$categorical_features,
                                                 years=IndividualData$years,
                                                 input_time_granularity=IndividualData$input_time_granularity,
                                                 conversion_factor=IndividualData$conversion_factor)


  latest_observed <- pkg.env$latest_observed_values_i(
    data_reserve= bind_rows(IndividualData$training.data, missing.obsevations),
    groups = hazard_frame_grouped$groups,
    categorical_features = IndividualData$categorical_features,
    continuous_features = IndividualData$continuous_features,
    calendar_period_extrapolation = IndividualData$calendar_period_extrapolation
  )

  max_DP <- max(bind_rows(IndividualData$training.data, missing.obsevations)$DP_rev_o)


  expected_i <- pkg.env$predict_i(
    hazard_data_frame = hazard_frame_grouped$hazard_group,
    latest_cumulative = latest_observed$latest_cumulative,
    grouping_method = "exposure"
  )

  df_i <- pkg.env$retrieve_df_i(
    hazard_data_frame = hazard_frame_grouped$hazard_group,
    groups = hazard_frame_grouped$groups
  )

  hazard_frame_input <- pkg.env$input_hazard_frame(
    hazard_frame = hazard_frame_grouped$hazard_group,
    expected_i = expected_i ,
    categorical_features = IndividualData$categorical_features,
    continuous_features = IndividualData$continuous_features,
    df_i = df_i,
    groups = hazard_frame_grouped$groups)

  ############################################################
  # If input dimension not equal to output dimension, we perform a grouping #

  if(IndividualData$conversion_factor != 1){

  development_periods <- distinct(select(data.frame(IndividualData$training), AP_i, AP_o))

  # Calculate the minimum and maximum development periods for each row in development_periods for each DP_rev_o
  dp_ranges <- t(lapply(1:max_DP, function(DP_rev_o) {
    cbind(DP_rev_o, development_periods,
          min_dp = with(development_periods, AP_i+1/(IndividualData$conversion_factor)*(DP_rev_o-AP_o)),
          max_dp = with(development_periods, AP_i-1+1/(IndividualData$conversion_factor)*(DP_rev_o-AP_o+1)))
  }
  ))

  dp_ranges <- do.call(rbind, dp_ranges)


  #check_input_hazard <- pkg.env$check_input_hazard(hazard_frame_input,
  #                           check_value=check_value)

  #If we exeed the check value, we calculate on ouput granularity, predict on output granularity, and distribute evenly in the relevant input-periods.
  #From here we do simple chain-ladder to calculate new development factor.
  if(#check_input_hazard
    FALSE
     ){
    development_factor_o <- mapply(pkg.env$i_to_o_development_factor,
                                   1:max(hazard_frame_grouped$groups$group_o),
                                   MoreArgs=list(hazard_data_frame=hazard_frame_grouped$hazard_group,
                                                 expected_i = expected_i,
                                                 dp_ranges = dp_ranges,
                                                 groups = hazard_frame_grouped$groups,
                                                 observed_pr_dp = latest_observed$observed_pr_dp,
                                                 latest_cumulative = latest_observed$latest_cumulative,
                                                 conversion_factor = IndividualData$conversion_factor,
                                                 grouping_method = "probability"))

    if(ncol(hazard_frame_grouped$groups) == 5){
      colnames(development_factor_o) <- unique(c(paste0("AP_o_",hazard_frame_grouped$groups$AP_o,",", hazard_frame_grouped$groups$covariate )))
    }
    else{
      colnames(development_factor_o) <- c(hazard_frame_grouped$groups$covariate )
    }

    df_o <- as.data.frame(development_factor_o[1:(nrow(development_factor_o)-1),]) %>%
      map_df(rev) %>%
      mutate(DP_o=row_number())

    #We only update for relevant periods, hence for example for accident periods, where we have already seen the development, we just put to 1.
    hazard_frame_grouped$hazard_group <- pkg.env$update_hazard_frame(
      hazard_frame_input=hazard_frame_input,
      hazard_frame_grouped=hazard_frame_grouped$hazard_group,
      df_o=df_o,
      latest_observed_i = latest_observed$observed_pr_dp,
      groups = hazard_frame_grouped$groups,
      categorical_features=IndividualData$categorical_features,
      continuous_features = IndividualData$continuous_features,
      conversion_factor = IndividualData$conversion_factor,
      check_value = check_value
    )

    expected_i <- pkg.env$predict_i(
      hazard_data_frame = hazard_frame_grouped$hazard_group,
      latest_cumulative = latest_observed$latest_cumulative,
      grouping_method = "exposure"
    )

    df_i <- pkg.env$retrieve_df_i(
      hazard_data_frame = hazard_frame_grouped$hazard_group,
      groups = hazard_frame_grouped$groups,
      adjusted=T
    )

    hazard_frame_input <- pkg.env$input_hazard_frame(
      hazard_frame = hazard_frame_grouped$hazard_group,
      expected_i = expected_i ,
      categorical_features = IndividualData$categorical_features,
      continuous_features = IndividualData$continuous_features,
      df_i = df_i,
      groups = hazard_frame_grouped$groups,
      adjusted=T)



  }

  expected_o <-pkg.env$predict_o(expected_i = expected_i,
                                 groups = hazard_frame_grouped$groups,
                                 conversion_factor = IndividualData$conversion_factor)



  development_factor_o <- mapply(pkg.env$i_to_o_development_factor,
                           1:max(hazard_frame_grouped$groups$group_o),
                           MoreArgs=list(hazard_data_frame=hazard_frame_grouped$hazard_group,
                                         expected_i = expected_i,
                                         dp_ranges = dp_ranges,
                                         groups = hazard_frame_grouped$groups,
                                         observed_pr_dp = latest_observed$observed_pr_dp,
                                         latest_cumulative = latest_observed$latest_cumulative,
                                         conversion_factor = IndividualData$conversion_factor,
                                         grouping_method = grouping_method))



#We only have 5 groups if AP is included as a covariate
  if(ncol(hazard_frame_grouped$groups) == 5){
    colnames(development_factor_o) <- unique(c(paste0("AP_o_",hazard_frame_grouped$groups$AP_o,",", hazard_frame_grouped$groups$covariate )))
  }
  else{
    colnames(development_factor_o) <- c(hazard_frame_grouped$groups$covariate )
  }

  df_o <- as.data.frame(development_factor_o[1:(nrow(development_factor_o)-1),]) %>%
    map_df(rev) %>%
    mutate(DP_o=row_number())


  hazard_frame_output <- pkg.env$output_hazard_frame(
    hazard_frame_input=hazard_frame_input,
    expected_o=expected_o,
    categorical_features=IndividualData$categorical_features,
    continuous_features=IndividualData$continuous_features,
    df_o=df_o,
    groups = hazard_frame_grouped$groups
  )

  out=list(df_output = df_o,
           df_input = df_i,
           hazard_frame_input = hazard_frame_input,
           hazard_frame_output = hazard_frame_output,
           IndividualData=IndividualData)

  class(out) <- c('ReSurvFit')

  return(out)
  }

  out=list(model.out=model.out,
           df_input = df_i,
           hazard_frame_input = hazard_frame_input,
           IndividualData=IndividualData)

  class(out) <- c('ReSurvFit')

  return(out)
}

#' Draft for plot of \code{ReSurvFit} models for simulated data.

# plot.ReSurvFit <- function(ReSurv){
#
#   "
#   Plot different development patterns. Currently only works for simulations with no accident-period dependency.
#   "
#   warning("Plotting functionality is not implemented yet.")
  # df_output<- ReSurv$df_output
  # df_input <- ReSurv$df_input

  #Theoretichal,assuming equal claim reporting throughout year.(kinda unrealistic)
  # df_output$DP_i <-(df_output$DP_o * 1/resurv.fit.xgboost$IndividualData$conversion_factor)-1/2*1/individual_data$conversion_factor+1/2
  #
  #
  # graph <-  df_output %>%  mutate(type = "Output") %>% select(-c(DP_o)) %>%
  #   reshape2::melt(id.vars =c("DP_i", "type"))  %>%
  #   rbind(
  #     df_input %>%  mutate(type = "Input") %>%
  #       reshape2::melt(id.vars =c("DP_i", "type"))
  #   ) %>%  group_by(type,variable) %>%  arrange(DP_i) %>% mutate(
  #     value_cum= 1/rev(cumprod(rev(value)))
  #   )

  # #Observed based
  # input_expected <- ReSurv$hazard_frame_input %>%
  #   left_join(ReSurv$IndividualData$training.data[,c("AP_i","DP_rev_i","claim_type", "I")] %>% group_by(AP_i, DP_rev_i, claim_type) %>%  summarize(I=sum(I)), by=c("AP_i","DP_rev_i","claim_type")) %>%
  #   mutate(DP_i = max(ReSurv$hazard_frame_input$DP_rev_i)-DP_rev_i + 1) %>%
  #   filter(DP_i <= max(ReSurv$hazard_frame_input$DP_rev_i)-AP_i+1) %>%
  #   mutate(DP_i_o = ceiling(DP_i*resurv.fit.xgboost$IndividualData$conversion_factor) ) %>%
  #   group_by(DP_i_o, claim_type) %>% summarize(I_E = sum(I_expected, na.rm=T)) %>%  ungroup() %>%
  #   group_by(claim_type) %>%
  #   arrange(DP_i_o) %>%  mutate(I_E_i=cumsum(I_E))
  #
  # output_expected <-ReSurv$hazard_frame_output %>%
  #   left_join(ReSurv$IndividualData$training.data[,c("AP_o","DP_rev_o","claim_type", "I")] %>% group_by(AP_o, DP_rev_o, claim_type) %>%  summarize(I=sum(I)), by=c("AP_o","DP_rev_o","claim_type")) %>%
  #   mutate(DP_o = max(ReSurv$hazard_frame_output$DP_rev_o)-DP_rev_o + 1) %>%
  #   filter(DP_o <= max(ReSurv$hazard_frame_output$DP_rev_o)-AP_o+1) %>%
  #   group_by(DP_o,claim_type) %>% summarize(I_E = sum(I_expected, na.rm=T)) %>%  ungroup() %>%
  #   group_by(claim_type) %>%
  #   arrange(DP_o) %>%  mutate(I_E_o=cumsum(I_E))
  #
  # DP_correction <- output_expected %>% left_join(input_expected,
  #                                            by=c("DP_o" = "DP_i_o",
  #                                                 "claim_type")) %>%
  #   mutate(correction_factor = I_E_o / I_E_i) %>%
  #   select(DP_o, claim_type, correction_factor) %>%
  #   mutate(DP_o = DP_o * 1/ReSurv$IndividualData$conversion_factor)
  #
  # df_output$DP_i <-(df_output$DP_o * 1/ReSurv$IndividualData$conversion_factor)
  #
  #
  #
  # graph <-  df_output %>%  mutate(type = "Output") %>% select(-c(DP_o)) %>%
  #   reshape2::melt(id.vars =c("DP_i", "type"))  %>%
  #   rbind(
  #     df_input %>%  mutate(type = "Input") %>%
  #       reshape2::melt(id.vars =c("DP_i", "type"))
  #   ) %>%  group_by(type,variable) %>%  arrange(DP_i) %>% mutate(
  #     value_cum= 1/rev(cumprod(rev(value)))
  #   ) %>%
  #   mutate(claim_type = substr(variable, 12,12)) %>%
  #   left_join(
  #     DP_correction, by =c("DP_i"= "DP_o", "claim_type")
  #   ) %>%
  #   mutate(DP_i = case_when(
  #     type == "Output" ~ DP_i *correction_factor,
  #     TRUE ~ DP_i))


  # plot.df <- ggplot(data=graph) +geom_line(aes(x=DP_i, y=value_cum, col=type)) + facet_grid(~variable) +
  #   theme_bw() +  labs(title = "Estimated claim development", x = "Input development period",
  #                      y = "Development percentage", color = "Granularity") +
  #   scale_y_continuous(labels = scales::percent_format(accuracy = 1))
  #
  #
  #
  # true <- ReSurv$IndividualData$full.data
  # conversion_factor <- ReSurv$IndividualData$conversion_factor
  #
  # true_output <- true %>%
  #   filter(DP_rev_i <= TR_i) %>%
  #   mutate(
  #     DP_rev_o = floor(max(DP_i)*conversion_factor)-ceiling(DP_i*conversion_factor+((AP_i-1)%%(1/conversion_factor))*conversion_factor) +1,
  #     AP_o = ceiling(AP_i*conversion_factor)
  #   ) %>%
  #   group_by(claim_type, AP_o, DP_rev_o) %>%
  #   mutate(claim_type = as.character(claim_type)) %>%
  #   summarize(I=sum(I), .groups = "drop")
  #
  # true_input <- true %>%
  #   filter(DP_rev_i <= TR_i) %>%
  #   group_by(claim_type, AP_i, DP_rev_i) %>%
  #   mutate(claim_type = as.character(claim_type)) %>%
  #   summarize(I=sum(I), .groups = "drop")
  #
  # output_triangle_plot <-ReSurv$hazard_frame_output[,c("claim_type","AP_o", "DP_rev_o", "I_expected")] %>%
  #   inner_join(true_output, by =c("claim_type","AP_o", "DP_rev_o")) %>%
  #   mutate(expected_ratio = (I-I_expected)/I) %>%
  #   mutate(DP_o = max(ReSurv$hazard_frame_output$DP_rev_o)-DP_rev_o + 1)
  #
  # input_triangle_plot <-ReSurv$hazard_frame_input[,c("claim_type","AP_i", "DP_rev_i", "I_expected")] %>%
  #   inner_join(true_input, by =c("claim_type","AP_i", "DP_rev_i")) %>%
  #   mutate(expected_ratio = (I-I_expected)/I) %>%
  #   mutate(DP_i= max(ReSurv$hazard_frame_input$DP_rev_i)-DP_rev_i + 1)
  #
  # plot.triangle.ratio.output <- ggplot(data = output_triangle_plot, aes(AP_o, DP_o, fill = expected_ratio))+
  #   facet_grid(~claim_type) +
  #   scale_y_reverse() +
  #   geom_tile(color = "white")+
  #   scale_fill_gradient2(low = "blue", high = "red", mid = "green",
  #                        midpoint = 0, limit = c(-1,1), space = "Lab", oob=scales::squish,
  #                        name="(Actual-Expected)/Actual") +
  #   theme_bw()+
  #   theme(axis.text.x = element_text(angle = 45, vjust = 1,
  #                                    size = 12, hjust = 1))+
  #   coord_fixed() + ggtitle("Expected ratio on unseen data") +
  #   xlab("Development Output Period") + ylab("Accident Output Period")
  #
  # plot.triangle.ratio.input <- ggplot(data = input_triangle_plot, aes(AP_i, DP_i, fill = expected_ratio))+
  #   facet_grid(~claim_type) +
  #   scale_y_reverse() +
  #   geom_tile(color = "white")+
  #   scale_fill_gradient2(low = "blue", high = "red", mid = "green",
  #                        midpoint = 0, limit = c(-2,2), space = "Lab", oob=scales::squish,
  #                        name="(Actual-Expected)/Actual") +
  #   theme_bw()+
  #   theme(axis.text.x = element_text(angle = 45, vjust = 1,
  #                                    size = 12, hjust = 1))+
  #   coord_fixed() + ggtitle("Expected ratio on unseen data") +
  #   xlab("Development Input Period") + ylab("Accident Input Period")
  #
  # return(list(plot.df,plot.triangle.ratio))

# }




