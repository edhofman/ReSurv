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
                   percentage_data_training=.8){

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
  # Om.df <-   pkg.env$create.om.df(training.data=IndividualData$training.data,
                                  # input_time_granularity=IndividualData$input_time_granularity,
                                  # years=IndividualData$years)



  if(hazard_model=="cox"){


    data=IndividualData$training.data

    X=data %>%
      select(c(IndividualData$continuous_features,IndividualData$categorical_features))

    Y=data[,c("DP_rev_i", "I", "TR_i")]

    model.out <- pkg.env$fit_cox_model(data=data,
                                       formula_ct=formula_ct,
                                       newdata=newdata)

    # tmp <- pkg.env$spline_hp(hparameters,IndividualData)


    ## OLD BASELINE COMPUTATION (BRESLOW)
    # bs_hazard <- basehaz( model.out$cox, centered=FALSE) %>%
    #   mutate(hazard = hazard-lag(hazard,default=0))
    #
    #
    # bsln <- data.frame(baseline=bs_hazard$hazard,
    #                    DP_rev_i=ceiling(bs_hazard$time))  #$hazard

    ## NEW BASELINE COMPUTATION (RESURV)

    X_tmp_bsln <- pkg.env$model.matrix.creator(data= IndividualData$training.data,
                                      select_columns = IndividualData$categorical_features,
                                      remove_first_dummy=T)

    scaler <- pkg.env$scaler(continuous_features_scaling_method = continuous_features_scaling_method)

    Xc_tmp_bsln <- IndividualData$training.data %>%
      reframe(across(all_of(IndividualData$continuous_features),
                     scaler))

    # training_test_split = pkg.env$check.traintestsplit(percentage_data_training)

    X_tmp_bsln=cbind(X_tmp_bsln,Xc_tmp_bsln)

    bsln <- pkg.env$baseline.calc(hazard_model = hazard_model,
                                  model.out = model.out,
                                  X=X_tmp_bsln,
                                  Y=Y)


    bsln <- data.frame(baseline=bsln,
                       DP_rev_i=sort(as.integer(unique(IndividualData$training.data$DP_rev_i))))

    ### make it relative
    newdata.bs <- ReSurv:::pkg.env$df.2.fcst.nn.pp(data=IndividualData$training.data,
                                                   newdata=newdata,
                                                   continuous_features=IndividualData$continuous_features,
                                                   categorical_features=IndividualData$categorical_features)

    benchmark_id <- ReSurv:::pkg.env$benchmark_id(X = X_tmp_bsln,
                                                  Y =Y ,
                                                  newdata.mx = newdata.bs,
                                                  remove_first_dummy=T)


    pred_relative <- model.out$cox_lp-model.out$cox_lp[benchmark_id]

    ###

    hazard_frame <- cbind(newdata, exp(pred_relative))
    colnames(hazard_frame)[dim(hazard_frame)[2]]="expg"



    is_lkh <- pkg.env$evaluate_lkh_cox(X_train=X,
                                    Y_train=Y,
                                    model=model.out)


    os_lkh <- NULL


  }

  if(hazard_model=="deepsurv"){

    X <- pkg.env$model.matrix.creator(data= IndividualData$training.data,
                                      select_columns = IndividualData$categorical_features)

    scaler <- pkg.env$scaler(continuous_features_scaling_method=continuous_features_scaling_method)

    Xc <- IndividualData$training.data %>%
      reframe(across(all_of(IndividualData$continuous_features),
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



    is_lkh <- pkg.env$evaluate_lkh_nn(X_train=datads_pp$lkh_eval_data$data_train,
                                   Y_train=datads_pp$lkh_eval_data$y_train,
                                   model=model.out)

    os_lkh <- pkg.env$evaluate_lkh_nn(X_train=datads_pp$lkh_eval_data$data_val,
                                      Y_train=datads_pp$lkh_eval_data$y_val,
                                      model=model.out)

  }

  if(hazard_model == "xgboost"){

    X <- pkg.env$model.matrix.creator(data= IndividualData$training.data,
                                      select_columns = IndividualData$categorical_features,
                                      remove_first_dummy=T)

    scaler <- pkg.env$scaler(continuous_features_scaling_method = continuous_features_scaling_method)

    Xc <- IndividualData$training.data %>%
      reframe(across(all_of(IndividualData$continuous_features),
                       scaler))

    training_test_split = pkg.env$check.traintestsplit(percentage_data_training)

    X=cbind(X,Xc)

    Y=IndividualData$training.data[,c("DP_rev_i", "I", "TR_i")]

    datads_pp <- pkg.env$xgboost_pp(X=X,
                                    Y=Y,
                                    training_test_split=training_test_split)

    model.out <- pkg.env$fit_xgboost(datads_pp,
                                     hparameters=hparameters)

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

    # compute the likelihood of the fitted model (upper triangle)
    is_lkh <- pkg.env$evaluate_lkh_xgb(X_train=X,
                                       Y_train=Y,
                                       dset='is',
                                       samples_cn=datads_pp$samples_cn,
                                       model=model.out)

    os_lkh <- pkg.env$evaluate_lkh_xgb(X_train=X,
                                       Y_train=Y,
                                       dset='os',
                                       samples_cn=datads_pp$samples_cn,
                                       model=model.out)

  }

  if(hazard_model == "LTRCtrees"){

    X <- pkg.env$model.matrix.creator(data= IndividualData$training.data,
                                      select_columns = IndividualData$categorical_features,
                                      remove_first_dummy=T)

    scaler <- pkg.env$scaler(continuous_features_scaling_method = continuous_features_scaling_method)

    Xc <- IndividualData$training.data %>%
      reframe(across(all_of(IndividualData$continuous_features),
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


    pred <- predict(model.out$cox,newdata)

    benchmark_id <- 1
    # pred_relative <- model.out$expg/model.out$expg[benchmark_id]

    pred_relative <- exp(pred - pred[benchmark_id])

    # exp(pred_relative)

    hazard_frame <- cbind(newdata, expg=pred_relative)
    # colnames(hazard_frame)[dim(hazard_frame)[2]]="expg"

    bsln <- data.frame(baseline=bsln,
                       DP_rev_i=sort(as.integer(unique(IndividualData$training.data$DP_rev_i))))



    is_lkh <- pkg.env$evaluate_lkh_LTRCtrees(X_train=IndividualData$training.data %>% select(c(IndividualData$categorical_features,IndividualData$continuous_features)),
                                    Y_train=Y,
                                    model=model.out)

    os_lkh <- NULL

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
                                                    # Om.df=Om.df,
                                                    categorical_features = IndividualData$categorical_features,
                                                    continuous_features = IndividualData$continuous_features,
                                                    calendar_period_extrapolation = IndividualData$calendar_period_extrapolation)


  out=list(model.out=list(data=X,
                          model.out=model.out),
           # Om.df=Om.df,
           is_lkh=is_lkh,
           os_lkh=os_lkh,
           hazard_frame = hazard_frame_updated,
           hazard_model = hazard_model,
           IndividualData = IndividualData)

  class(out) <- c('ReSurvFit')

  return(out)
}




