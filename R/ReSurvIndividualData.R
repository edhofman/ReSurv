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
#' \item{\href{https://cran.r-project.org/web/packages/survival/vignettes/survival.pdf}{COX}}
#' \item{\href{https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-018-0482-1}{Neural Networks}}
#' \item{\href{https://xgboost.readthedocs.io/en/stable/}{eXtreme Gradient Boosting}}
#' }
#'
#'
#' @param IndividualDataPP IndividualDataPP object to use for the \code{ReSurv} fit.
#' @param hazard_model \code{character}, hazard model supported from our package, must be provided as a string. The model can be chosen from:
#' \itemize{
#' \item{\code{"COX"}: Standard Cox model for the hazard.}
#' \item{\code{"NN"}: Deep Survival Neural Network.}
#' \item{\code{"XGB"}: eXtreme Gradient Boosting.}
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
#' @return \code{ReSurv} fit. A list containing
#' \itemize{
#' \item{\code{model.out}: \code{list} containing the pre-processed covariates data for the fit (\code{data}) and the basic model output (\code{model.out};COX, XGB or NN).}
#' \item{\code{is_lkh}: \code{numeric} Training negative log likelihood.}
#' \item{\code{os_lkh}:  \code{numeric} Validation  negative log likelihood. Not available for COX.}
#' \item{\code{hazard_frame}: \code{data.frame} containing the fitted hazard model with the corresponding covariates. It contains:}
#'    \itemize{
#'    \item{\code{expg}: fitted risk score.}
#'    \item{\code{baseline}: fitted baseline.}
#'    \item{\code{hazard}: fitted hazard rate (\code{expg}*\code{baseline}).}
#'    \item{\code{f_i}: fitted development factors.}
#'    \item{\code{cum_f_i}: fitted cumulative development factors.}
#'    \item{\code{S_i}:fitted survival function.}
#'    \item{\code{S_i_lag}:fitted survival function (lag version, for further information see \code{?dplyr::lag}).}
#'    \item{\code{S_i_lead}:fitted survival function (lead version, for further information see \code{?dplyr::lead}).}
#'    }
#' \item{\code{hazard_model}: \code{string} chosen hazard model (COX, NN or XGB)}
#' \item{\code{IndividualDataPP}: starting \code{IndividualDataPP} object.}
#' }
#'
#' @import reticulate
#' @import tidyverse
#' @import xgboost
#' @import rpart
#' @import data.table
#' @importFrom dplyr reframe full_join
#' @importFrom tidyr replace_na
#'
#' @references
#' Munir, H., Emil, H., & Gabriele, P. (2023). A machine learning approach based on survival analysis for IBNR frequencies in non-life reserving. arXiv preprint arXiv:2312.14549.
#'
#' Therneau, T. M., & Lumley, T. (2015). Package ‘survival’. R Top Doc, 128(10), 28-33.
#'
#' Katzman, J. L., Shaham, U., Cloninger, A., Bates, J., Jiang, T., & Kluger, Y. (2018). DeepSurv: personalized treatment recommender system using a Cox proportional hazards deep neural network. BMC medical research methodology, 18(1), 1-12.
#'
#' Chen, T., He, T., Benesty, M., & Khotilovich, V. (2019). Package ‘xgboost’. R version, 90, 1-66.
#'
#' @export
ReSurv <- function(IndividualDataPP,
                   hazard_model="COX",
                   tie='efron',
                   baseline="spline",
                   continuous_features_scaling_method="minmax",
                   random_seed=1,
                   hparameters=list(),
                   percentage_data_training=.8,
                   grouping_method = "exposure",
                   check_value = 1.85){

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
#' \item{\href{https://cran.r-project.org/web/packages/survival/vignettes/survival.pdf}{COX}}
#' \item{\href{https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-018-0482-1}{Neural Networks}}
#' \item{\href{https://xgboost.readthedocs.io/en/stable/}{eXtreme Gradient Boosting}}
#' }
#'
#'
#' @param IndividualDataPP IndividualDataPP object to use for the \code{ReSurv} fit.
#' @param hazard_model \code{character}, hazard model supported from our package, must be provided as a string. The model can be chosen from:
#' \itemize{
#' \item{\code{"COX"}: Standard Cox model for the hazard.}
#' \item{\code{"NN"}: Deep Survival Neural Network.}
#' \item{\code{"XGB"}: eXtreme Gradient Boosting.}
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
#' @return \code{ReSurv} fit. A list containing
#' \itemize{
#' \item{\code{model.out}: \code{list} containing the pre-processed covariates data for the fit (\code{data}) and the basic model output (\code{model.out};COX, XGB or NN).}
#' \item{\code{is_lkh}: \code{numeric} Training negative log likelihood.}
#' \item{\code{os_lkh}:  \code{numeric} Validation  negative log likelihood. Not available for COX.}
#' \item{\code{hazard_frame}: \code{data.frame} containing the fitted hazard model with the corresponding covariates. It contains:}
#'    \itemize{
#'    \item{\code{expg}: fitted risk score.}
#'    \item{\code{baseline}: fitted baseline.}
#'    \item{\code{hazard}: fitted hazard rate (\code{expg}*\code{baseline}).}
#'    \item{\code{f_i}: fitted development factors.}
#'    \item{\code{cum_f_i}: fitted cumulative development factors.}
#'    \item{\code{S_i}:fitted survival function.}
#'    \item{\code{S_i_lag}:fitted survival function (lag version, for further information see \code{?dplyr::lag}).}
#'    \item{\code{S_i_lead}:fitted survival function (lead version, for further information see \code{?dplyr::lead}).}
#'    }
#' \item{\code{hazard_model}: \code{string} chosen hazard model (COX, NN or XGB)}
#' \item{\code{IndividualDataPP}: starting \code{IndividualDataPP} object.}
#' }
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
ReSurv.default <- function(IndividualDataPP,
                           hazard_model="COX",
                           tie='efron',
                           baseline="spline",
                           continuous_features_scaling_method="minmax",
                           random_seed=1,
                           hparameters=list(),
                           percentage_data_training=.8,
                           grouping_method = "exposure",
                           check_value = 1.85){

  message('The object provided must be of class IndividualDataPP')

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
#' \item{\href{https://cran.r-project.org/web/packages/survival/vignettes/survival.pdf}{COX}}
#' \item{\href{https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-018-0482-1}{Neural Networks}}
#' \item{\href{https://xgboost.readthedocs.io/en/stable/}{eXtreme Gradient Boosting}}
#' }
#'
#'
#' @param IndividualDataPP IndividualDataPP object to use for the \code{ReSurv} fit.
#' @param hazard_model \code{character}, hazard model supported from our package, must be provided as a string. The model can be chosen from:
#' \itemize{
#' \item{\code{"COX"}: Standard Cox model for the hazard.}
#' \item{\code{"NN"}: Deep Survival Neural Network.}
#' \item{\code{"XGB"}: eXtreme Gradient Boosting.}
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
#' @return \code{ReSurv} fit. A list containing
#' \itemize{
#' \item{\code{model.out}: \code{list} containing the pre-processed covariates data for the fit (\code{data}) and the basic model output (\code{model.out};COX, XGB or NN).}
#' \item{\code{is_lkh}: \code{numeric} Training negative log likelihood.}
#' \item{\code{os_lkh}:  \code{numeric} Validation  negative log likelihood. Not available for COX.}
#' \item{\code{hazard_frame}: \code{data.frame} containing the fitted hazard model with the corresponding covariates. It contains:}
#'    \itemize{
#'    \item{\code{expg}: fitted risk score.}
#'    \item{\code{baseline}: fitted baseline.}
#'    \item{\code{hazard}: fitted hazard rate (\code{expg}*\code{baseline}).}
#'    \item{\code{f_i}: fitted development factors.}
#'    \item{\code{cum_f_i}: fitted cumulative development factors.}
#'    \item{\code{S_i}:fitted survival function.}
#'    \item{\code{S_i_lag}:fitted survival function (lag version, for further information see \code{?dplyr::lag}).}
#'    \item{\code{S_i_lead}:fitted survival function (lead version, for further information see \code{?dplyr::lead}).}
#'    }
#' \item{\code{hazard_model}: \code{string} chosen hazard model (COX, NN or XGB)}
#' \item{\code{IndividualDataPP}: starting \code{IndividualDataPP} object.}
#' }
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
ReSurv.IndividualDataPP <- function(IndividualDataPP,
                                  hazard_model="COX",
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

  formula_ct <- as.formula(IndividualDataPP$string_formula_i)

  newdata <- create.df.2.fcst(IndividualDataPP=IndividualDataPP,
                              hazard_model=hazard_model)


  # logical: check if we work with a baseline model
  is_baseline_model = is.null(c(IndividualDataPP$categorical_features,
                                IndividualDataPP$continuous_features))


  # create data frame of occurrencies to weight development factors
  # Om.df <-   pkg.env$create.om.df(training.data=IndividualDataPP$training.data,
                                  # input_time_granularity=IndividualDataPP$input_time_granularity,
                                  # years=IndividualDataPP$years)



  if(hazard_model=="COX"){

    data=IndividualDataPP$training.data

    X=data %>%
      select(c(IndividualDataPP$continuous_features,IndividualDataPP$categorical_features))

    Y=data[,c("DP_rev_i", "I", "TR_i")]

    model.out <- pkg.env$fit_cox_model(data=data,
                                       formula_ct=formula_ct,
                                       newdata=newdata)

    # tmp <- pkg.env$spline_hp(hparameters,IndividualDataPP)


    ## OLD BASELINE COMPUTATION (BRESLOW)
    # bs_hazard <- basehaz( model.out$cox, centered=FALSE) %>%
    #   mutate(hazard = hazard-lag(hazard,default=0))
    #
    #
    # bsln <- data.frame(baseline=bs_hazard$hazard,
    #                    DP_rev_i=ceiling(bs_hazard$time))  #$hazard

    ## NEW BASELINE COMPUTATION (RESURV)

    if(is_baseline_model){

      X_tmp_bsln = data.frame(rep(1,dim(Y)[1]))

    }else{

      scaler <- pkg.env$scaler(continuous_features_scaling_method = continuous_features_scaling_method)

      Xc_tmp_bsln <- IndividualDataPP$training.data %>%
      reframe(across(all_of(IndividualDataPP$continuous_features),
                     scaler))


    if(!is.null(IndividualDataPP$categorical_features)){


      X_tmp_bsln <- pkg.env$model.matrix.creator(data= IndividualDataPP$training.data,
                                      select_columns = IndividualDataPP$categorical_features,
                                      remove_first_dummy=T)


      X_tmp_bsln=cbind(X_tmp_bsln,Xc_tmp_bsln)

    }else{

      X_tmp_bsln= Xc_tmp_bsln

    }

      }

    bsln <- pkg.env$baseline.calc(hazard_model = hazard_model,
                                  model.out = model.out,
                                  X=X_tmp_bsln,
                                  Y=Y)


    bsln <- data.frame(baseline=bsln,
                       DP_rev_i=sort(as.integer(unique(IndividualDataPP$training.data$DP_rev_i))))

    ### make it relative

    if(is_baseline_model){

      newdata.bs <- data.frame(intercept_1 = rep(1, dim(newdata)[1]))

      benchmark_id <- pkg.env$benchmark_id(X = X_tmp_bsln,
                                           Y =Y ,
                                           newdata.mx = newdata.bs,
                                           remove_first_dummy=F)

    }else{

      newdata.bs <- pkg.env$df.2.fcst.nn.pp(data=IndividualDataPP$training.data,
                                                     newdata=newdata,
                                                     continuous_features=IndividualDataPP$continuous_features,
                                                     categorical_features=IndividualDataPP$categorical_features)

      benchmark_id <- pkg.env$benchmark_id(X = X_tmp_bsln,
                                           Y =Y ,
                                           newdata.mx = newdata.bs,
                                           remove_first_dummy=T)


      }





    pred_relative <- model.out$cox_lp-model.out$cox_lp[benchmark_id]

    ###

    hazard_frame <- cbind(newdata, exp(pred_relative))
    colnames(hazard_frame)[dim(hazard_frame)[2]]="expg"



    is_lkh <- pkg.env$evaluate_lkh_cox(X_train=X,
                                    Y_train=Y,
                                    model=model.out)


    os_lkh <- NULL


  }

  if(hazard_model=="NN"){


    Y=IndividualDataPP$training.data[,c("DP_rev_i", "I", "TR_i")]

    training_test_split = pkg.env$check.traintestsplit(percentage_data_training)

    if(is_baseline_model){

      X <- data.frame(intercept_1 = rep(1,dim(Y)[1]))

    }else{

      scaler <- pkg.env$scaler(continuous_features_scaling_method=continuous_features_scaling_method)

      Xc <- IndividualDataPP$training.data %>%
        reframe(across(all_of(IndividualDataPP$continuous_features),
                       scaler))

      if(!is.null(IndividualDataPP$categorical_features)){

        X <- pkg.env$model.matrix.creator(data= IndividualDataPP$training.data,
                                      select_columns = IndividualDataPP$categorical_features)

        X = cbind(X,Xc)


        }else{


          X <- Xc

        }

    }

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

    if(is_baseline_model){

      newdata.mx <- data.frame(intercept_1= rep(1,dim(newdata)[1]))

    }else{

    newdata.mx <- pkg.env$df.2.fcst.nn.pp(data=IndividualDataPP$training.data,
                                          newdata=newdata,
                                          continuous_features=IndividualDataPP$continuous_features,
                                          categorical_features=IndividualDataPP$categorical_features)}



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
                       DP_rev_i=sort(as.integer(unique(IndividualDataPP$training.data$DP_rev_i))))



    if(!inherits(datads_pp$lkh_eval_data$data_train,"data.frame")){


      is_lkh <- pkg.env$evaluate_lkh_nn(X_train=as.data.frame(datads_pp$lkh_eval_data$data_train),
                                        Y_train=datads_pp$lkh_eval_data$y_train,
                                        model=model.out)

      os_lkh <- pkg.env$evaluate_lkh_nn(X_train=as.data.frame(datads_pp$lkh_eval_data$data_val),
                                        Y_train=datads_pp$lkh_eval_data$y_val,
                                        model=model.out)


    }else{

      is_lkh <- pkg.env$evaluate_lkh_nn(X_train=datads_pp$lkh_eval_data$data_train,
                                        Y_train=datads_pp$lkh_eval_data$y_train,
                                        model=model.out)

      os_lkh <- pkg.env$evaluate_lkh_nn(X_train=datads_pp$lkh_eval_data$data_val,
                                        Y_train=datads_pp$lkh_eval_data$y_val,
                                        model=model.out)

    }




  }

  if(hazard_model == "XGB"){

    Y=IndividualDataPP$training.data[,c("DP_rev_i", "I", "TR_i")]

    training_test_split = pkg.env$check.traintestsplit(percentage_data_training)

    if(is_baseline_model){

      X= data.frame(intercept_1 = rep(1,dim(Y)[1]))

    }else{

      scaler <- pkg.env$scaler(continuous_features_scaling_method = continuous_features_scaling_method)

      Xc <- IndividualDataPP$training.data %>%
        reframe(across(all_of(IndividualDataPP$continuous_features),
                       scaler))


      if(!is.null(IndividualDataPP$categorical_features)){

        X <- pkg.env$model.matrix.creator(data= IndividualDataPP$training.data,
                                      select_columns = IndividualDataPP$categorical_features,
                                      remove_first_dummy=T)

        X=cbind(X,Xc)
        }else{


          X <- Xc


        }

      }



    datads_pp <- pkg.env$xgboost_pp(X=X,
                                    Y=Y,
                                    training_test_split=training_test_split)

    model.out <- pkg.env$fit_xgboost(datads_pp,
                                     hparameters=hparameters)

    bsln <- pkg.env$baseline.calc(hazard_model = hazard_model,
                                  model.out = model.out,
                                  X=X,
                                  Y=Y)


    if(is_baseline_model){

      newdata.mx <- xgboost::xgb.DMatrix(as.matrix(rep(1, dim(newdata)[1])))

    }else{

      newdata.mx <- pkg.env$df.2.fcst.xgboost.pp(data=IndividualDataPP$training.data,
                                               newdata=newdata,
                                               continuous_features=IndividualDataPP$continuous_features,
                                               categorical_features=IndividualDataPP$categorical_features)
      }

    pred <- predict(model.out,newdata.mx)



    if(is_baseline_model){

      newdata.bs <- data.frame(intercept_1 = rep(1, dim(newdata)[1]))

      benchmark_id <- pkg.env$benchmark_id(X = X,
                                           Y =Y ,
                                           newdata.mx = newdata.bs,
                                           remove_first_dummy=F)

    }else{
    #make to hazard relative to initial model, to have similiar interpretation as standard cox
    newdata.bs <- pkg.env$df.2.fcst.nn.pp(data=IndividualDataPP$training.data,
                                          newdata=newdata,
                                          continuous_features=IndividualDataPP$continuous_features,
                                          categorical_features=IndividualDataPP$categorical_features)

    benchmark_id <- pkg.env$benchmark_id(X = X,
                                         Y =Y ,
                                         newdata.mx = newdata.bs,
                                         remove_first_dummy=T)}


    pred_relative <- pred - pred[benchmark_id]

    expg <- exp(pred_relative)

    hazard_frame <- cbind(newdata,expg)

    bsln <- data.frame(baseline=bsln,
                       DP_rev_i=sort(as.integer(unique(IndividualDataPP$training.data$DP_rev_i))))

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

  #hazard_q <- matrix(nrow=max_DP, ncol=(ncol(hazard)-1)*IndividualDataPP$conversion_factor)
  #eta_o <- c()


  ############################################################

  #Add development and relevant survival values to the hazard_frame
  hazard_frame_updated <- pkg.env$hazard_data_frame(hazard=hazard_frame,
                                                    # Om.df=Om.df,
                                                    categorical_features = IndividualDataPP$categorical_features,
                                                    continuous_features = IndividualDataPP$continuous_features,
                                                    calendar_period_extrapolation = IndividualDataPP$calendar_period_extrapolation)


  out_hz_frame <-  hazard_frame_updated %>%
      mutate(DP_i=pkg.env$maximum.time(IndividualDataPP$years, IndividualDataPP$input_time_granularity)-DP_rev_i+1) %>%
      relocate(DP_i, .after =  AP_i) %>%
      rename(f_i=dev_f_i,
             cum_f_i=cum_dev_f_i)

  out=list(model.out=list(data=X,
                          model.out=model.out),
           # Om.df=Om.df,
           is_lkh=is_lkh,
           os_lkh=os_lkh,
           hazard_frame = out_hz_frame,
           hazard_model = hazard_model,
           IndividualDataPP = IndividualDataPP)

  class(out) <- c('ReSurvFit')

  return(out)
}



