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
#' \item{\href{https://github.com/therneau/survival}{COX}}
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
#'
#'
#'@examples
#'
#' input_data_0 <- data_generator(
#' random_seed = 1964,
#' scenario = "alpha",
#' time_unit = 1,
#' years = 4,
#' period_exposure = 100)
#'
#' individual_data <- IndividualDataPP(data = input_data_0,
#' categorical_features = "claim_type",
#' continuous_features = "AP",
#' accident_period = "AP",
#' calendar_period = "RP",
#' input_time_granularity = "years",
#' output_time_granularity = "years",
#' years=4)
#'
#'
#' resurv_fit_cox <- ReSurv(individual_data,
#' hazard_model = "COX")
#'
#'
#'
#'

#' @import xgboost
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
                   hazard_model = "COX",
                   tie = "efron",
                   baseline = "spline",
                   continuous_features_scaling_method = "minmax",
                   random_seed = 1,
                   hparameters = list(),
                   percentage_data_training = .8,
                   grouping_method = "probability",
                   check_value = 1.85,
                   eta=0.5,
                   simplifier=FALSE){

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
#' \item{\href{https://github.com/therneau/survival}{COX}}
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


#' @import xgboost

#'
#'
#'
#'
#'@examples
#'
#' input_data_0 <- data_generator(
#' random_seed = 1964,
#' scenario = "alpha",
#' time_unit = 1,
#' years = 4,
#' period_exposure = 100)
#'
#' individual_data <- IndividualDataPP(data = input_data_0,
#' categorical_features = "claim_type",
#' continuous_features = "AP",
#' accident_period = "AP",
#' calendar_period = "RP",
#' input_time_granularity = "years",
#' output_time_granularity = "years",
#' years=4)
#'
#'
#' resurv_fit_cox <- ReSurv(individual_data,
#' hazard_model = "COX")
#'
#'
#'
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
                           hazard_model = "COX",
                           tie = "efron",
                           baseline = "spline",
                           continuous_features_scaling_method = "minmax",
                           random_seed = 1,
                           hparameters = list(),
                           percentage_data_training = .8,
                           grouping_method = "exposure",
                           check_value = 1.85,
                           eta=0.5,
                           simplifier=FALSE){

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
#' \item{\href{https://github.com/therneau/survival}{COX}}
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


#' @import xgboost

#'
#'
#'
#'
#'@examples
#'
#' input_data_0 <- data_generator(
#' random_seed = 1964,
#' scenario = "alpha",
#' time_unit = 1,
#' years = 4,
#' period_exposure = 100)
#'
#' individual_data <- IndividualDataPP(data = input_data_0,
#' categorical_features = "claim_type",
#' continuous_features = "AP",
#' accident_period = "AP",
#' calendar_period = "RP",
#' input_time_granularity = "years",
#' output_time_granularity = "years",
#' years=4)
#'
#'
#' resurv_fit_cox <- ReSurv(individual_data,
#' hazard_model = "COX")
#'
#'
#'
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
                                  hazard_model = "COX",
                                  tie = "efron",
                                  baseline = "spline",
                                  continuous_features_scaling_method = "minmax",
                                  random_seed = 1,
                                  hparameters = list(),
                                  percentage_data_training = .8,
                                  grouping_method = "exposure",
                                  check_value = 1.85,
                                  eta=0.5,
                                  simplifier=FALSE
){



  cont_f <- IndividualDataPP$data_information$continuous_features
  cat_f <- IndividualDataPP$data_information$categorical_features

  set.seed(random_seed)

  formula_ct <- as.formula(IndividualDataPP$data_information$string_formula_i)

  if(simplifier){

    columns_for_grouping <- unique(c(cont_f,cat_f,"AP_i"))

    tmp <- as.data.table(IndividualDataPP$training.data)

    out <- tmp[,.(.N),by=columns_for_grouping][,..columns_for_grouping]

    l4 <- list()

    l4$DP_rev_i <- min(IndividualDataPP$training.data[,'DP_rev_i']):max(IndividualDataPP$training.data[,'DP_rev_i'])

    l4<-do.call(CJ, c(l4, sorted = FALSE))

    newdata<-as.data.frame(setkey(out[,c(k=1,.SD)],k)[l4[,c(k=1,.SD)],allow.cartesian=TRUE][,k:=NULL])


  }else{
  newdata <- create.df.2.fcst(IndividualDataPP=IndividualDataPP,
                              hazard_model=hazard_model)}


  # logical: check if we work with a baseline model
  is_baseline_model = is.null(c(cont_f,
                                cat_f))


  if(hazard_model=="COX"){

    data=IndividualDataPP$training.data

    X=data[,.SD,
           .SDcols=c(cont_f,
                     cat_f)]

    # X=data %>%
    #   select(c(IndividualDataPP$data_information$continuous_features,IndividualDataPP$data_information$categorical_features))

    Y=data[,.SD,
           .SDcols=c("DP_rev_i", "I", "TR_i")]

    # model.out <- pkg.env$fit_cox_model(data=data,
    #                                    formula_ct=formula_ct,
    #                                    newdata=newdata)
    cox <- coxph(formula_ct, data=data, ties="efron")
    cox_lp <- predict(cox,newdata=newdata,'lp',reference='zero')

    cox_training_lp <- predict(cox,newdata=data[order(DP_rev_i)],'lp',reference='zero')

    model.out <- list(
      cox=cox,
      cox_lp=cox_lp,
      expg = exp(cox_lp),
      train_expg= cox_training_lp#exp(cox_training_lp)
    )
    ## NEW BASELINE COMPUTATION (RESURV)

    if(is_baseline_model){

      X_tmp_bsln = data.frame(rep(1,dim(Y)[1]))

    }else{




    if(!is.null(cat_f)){


      X_tmp_bsln <- pkg.env$model.matrix.creator(data= IndividualDataPP$training.data,
                                      select_columns = cat_f,
                                      remove_first_dummy=T)




      if(!is.null(cont_f)){
        scaler <- pkg.env$scaler(continuous_features_scaling_method = continuous_features_scaling_method)

        Xc_tmp_bsln<- IndividualDataPP$training.data[,lapply(.SD,scaler),.SDcols=IndividualDataPP$data_information$continuous_features]

        X_tmp_bsln=cbind(X_tmp_bsln,Xc_tmp_bsln)
      }


    }


    if(!is.null(cont_f) & is.null(cat_f)){

      scaler <- pkg.env$scaler(continuous_features_scaling_method = continuous_features_scaling_method)

      Xc_tmp_bsln<- IndividualDataPP$training.data[,lapply(.SD,scaler),.SDcols=IndividualDataPP$data_information$continuous_features]


      X_tmp_bsln= Xc_tmp_bsln

    }

      }

    bsln <- pkg.env$baseline.calc(hazard_model = hazard_model,
                                  model.out = model.out,
                                  X=X_tmp_bsln,
                                  Y=Y)


    bsln <- data.table(baseline=bsln,
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
                                                     continuous_features=IndividualDataPP$data_information$continuous_features,
                                                     categorical_features=IndividualDataPP$data_information$categorical_features)

      benchmark_id <- pkg.env$benchmark_id(X = X_tmp_bsln,
                                           Y =Y ,
                                           newdata.mx = newdata.bs,
                                           remove_first_dummy=T)


      }





    pred_relative <- model.out$cox_lp-model.out$cox_lp[benchmark_id]

    ###

    # hazard_frame <- cbind(newdata, exp(pred_relative))
    # colnames(hazard_frame)[dim(hazard_frame)[2]]="expg"

    hazard_frame = copy(newdata)
    setDT(hazard_frame)
    hazard_frame[,expg:=exp(pred_relative)]


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
        reframe(across(all_of(IndividualDataPP$data_information$continuous_features),
                       scaler))

      if(!is.null(IndividualDataPP$data_information$categorical_features)){

        X <- pkg.env$model.matrix.creator(data= IndividualDataPP$training.data,
                                      select_columns = IndividualDataPP$data_information$categorical_features)

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



    bsln <- pkg.env$baseline.calc(hazard_model = hazard_model,
                                  model.out = model.out,
                                  X=X,
                                  Y=Y)

    if(is_baseline_model){

      newdata.mx <- data.frame(intercept_1= rep(1,dim(newdata)[1]))

    }else{

    newdata.mx <- pkg.env$df.2.fcst.nn.pp(data=IndividualDataPP$training.data,
                                          newdata=newdata,
                                          continuous_features=IndividualDataPP$data_information$continuous_features,
                                          categorical_features=IndividualDataPP$data_information$categorical_features)
    }



    x_fc = as.matrix(newdata.mx)


    beta_ams <- pkg.env$predict_deepsurv(model.out$net, x_fc)

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
        reframe(across(all_of(IndividualDataPP$data_information$continuous_features),
                       scaler))


      if(!is.null(IndividualDataPP$data_information$categorical_features)){

        X <- pkg.env$model.matrix.creator(data= IndividualDataPP$training.data,
                                      select_columns = IndividualDataPP$data_information$categorical_features,
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
                                               continuous_features=IndividualDataPP$data_information$continuous_features,
                                               categorical_features=IndividualDataPP$data_information$categorical_features)
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
                                          continuous_features=IndividualDataPP$data_information$continuous_features,
                                          categorical_features=IndividualDataPP$data_information$categorical_features)

    benchmark_id <- pkg.env$benchmark_id(X = X,
                                         Y =Y ,
                                         newdata.mx = newdata.bs,
                                         remove_first_dummy=T)}


    pred_relative <- pred - pred[benchmark_id]

    expg <- exp(pred_relative)

    hazard_frame <- cbind(newdata,expg)

    bsln <- data.table(baseline=bsln,
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


  # hazard_frame <- hazard_frame %>%
  #   full_join(bsln,
  #             by="DP_rev_i") %>%
  #   as.data.frame() %>%
  #   replace_na(list(baseline=0))


  # assuming hazard_frame and bsln are already data.tables
  hazard_frame <- merge(
    hazard_frame,
    bsln,
    by = "DP_rev_i",
    all = TRUE
  )

  # replace NA in column "baseline" with 0
  hazard_frame[is.na(baseline), baseline := 0]


  # hazard_frame[,'hazard'] <- hazard_frame[,'baseline']*hazard_frame[,'expg']




  # From hazard to development factors

  hazard_frame[,hazard:=baseline*expg]

  hazard_frame[,dev_f_i := (1+(1-..eta)*hazard)/(1-..eta*hazard)]

  hazard_frame[is.na(dev_f_i), dev_f_i := 1]

  hazard_frame[dev_f_i < 0, dev_f_i := 1]

  hazard_frame <- hazard_frame[order(DP_rev_i)]

  # columns grouping
  group_cols <- c(unique(c("AP_i",IndividualDataPP$data_information$continuous_features)),IndividualDataPP$data_information$categorical_features)

  hazard_frame[, `:=`(
    cum_dev_f_i = cumprod(dev_f_i),
    S_i = fifelse(cumprod(dev_f_i) == 0, 0, 1 / cumprod(dev_f_i))
  ), by = group_cols]

  hazard_frame[, `:=`(
    S_i_lead = shift(S_i, type = "lead", fill = 0),
    S_i_lag  = shift(S_i, type = "lag",  fill = 1)
  ), by = group_cols]

  # hazard_frame[, c("expg", "baseline", "hazard") := NULL]


  cols_to_remove_na_from <- c("dev_f_i", "S_i", "S_i_lead", "S_i_lag", "cum_dev_f_i")

  hazard_frame[, (cols_to_remove_na_from) := lapply(.SD, function(x) fifelse(is.na(x), 1, x)), .SDcols = cols_to_remove_na_from]



  #Add development and relevant survival values to the hazard_frame
  # hazard_frame_updated <- pkg.env$hazard_data_frame(hazard=hazard_frame,
  #                                                   # Om.df=Om.df,
  #                                                   eta_old=eta,
  #                                                   categorical_features = IndividualDataPP$categorical_features,
  #                                                   continuous_features = IndividualDataPP$continuous_features,
  #                                                   calendar_period_extrapolation = IndividualDataPP$calendar_period_extrapolation)





 # prepare software output

  ## some columns ordering
  hazard_frame[,DP_i:=pkg.env$maximum.time(IndividualDataPP$data_information$years, IndividualDataPP$data_information$input_time_granularity)-DP_rev_i+1]

  cols <- names(hazard_frame)
  new_order <- append(cols[cols != "DP_i"], "DP_i", after = which(cols == "AP_i"))
  setcolorder(hazard_frame, new_order)

  ## some columns renaming
  setnames(hazard_frame, old = c("dev_f_i", "cum_dev_f_i"), new = c("f_i", "cum_f_i"))



  # Here we calculate the information that is needed for reserving - so we avoid carrying around copies of the main data.tables.
  ## 1. Some features combination are missing but required for prediction.
  ## 2. We need to find the range of development period required.
  #################################################################################################################
  data_information = IndividualDataPP$data_information

  max_dp_i = pkg.env$maximum.time(IndividualDataPP$data_information$years, IndividualDataPP$data_information$input_time_granularity)

  tmp.ls <- IndividualDataPP$training.data[
    (max_dp_i - DP_i + 1) > (AP_i - 1)
  ]


  setDT(tmp.ls)

  cols <- unique(c(IndividualDataPP$data_information$categorical_features,
            IndividualDataPP$data_information$continuous_features,
           "AP_i"))


  tmp.ls <- tmp.ls[,.(.N),by=cols][,.(DP_i=1:max_dp_i),by=cols] #


  tmp.existing <- unique(
    data[
      (pkg.env$maximum.time(IndividualDataPP$data_information$years, IndividualDataPP$data_information$input_time_granularity) - DP_i + 1) > (AP_i - 1),
      unique(c(IndividualDataPP$data_information$continuous_features, IndividualDataPP$data_information$categorical_features, "AP_i", "DP_i")),
      with = FALSE
    ]
  )



  test.missing <- dplyr::setdiff(x=tmp.ls,y=tmp.existing)

  if(dim(test.missing)[1]==0){
    tmp.missing <- NULL
  }else{

    tmp.missing <- copy(test.missing)

    tmp.missing[,c("DP_rev_i",
                   "TR_i",
                   "I"):=list(max_dp_i - DP_i + 1L,
                              AP_i - 1L,
                              0L)]

    tmp.missing<-tmp.missing[
      DP_rev_i > TR_i,
    ]


    tmp.missing[,
                c("DP_rev_o",
                  "AP_o"
                  ):=list(floor(max_dp_i * IndividualDataPP$data_information$conversion_factor) -
                                   ceiling(DP_i * IndividualDataPP$data_information$conversion_factor +
                                             ((AP_i - 1) %% (1 / IndividualDataPP$data_information$conversion_factor)) * IndividualDataPP$data_information$conversion_factor) + 1L,
                                 ceiling(AP_i * IndividualDataPP$data_information$conversion_factor)
                                 )

                ][,TR_o:=AP_o - 1L]


    tmp.missing[, (IndividualDataPP$data_information$categorical_features) := lapply(.SD, as.factor), .SDcols = IndividualDataPP$data_information$categorical_features]

    tmp.missing[,.SD,
                .SDcols = colnames(tmp.missing)%in%unique(c(
                  IndividualDataPP$data_information$categorical_features,
                  IndividualDataPP$data_information$continuous_features,
                  "AP_i","AP_o","DP_i","DP_rev_i","DP_rev_o","TR_i","TR_o","I"
                ))]





    }

    # the data that will be used for predictions is stored here.
    data_information$data_for_reserving <-bind_rows(IndividualDataPP$training.data, tmp.missing )

    ############################################################################

    # Also compute the ranges of development periods needed for later
    development_periods <- unique(IndividualDataPP$training.data[, .(AP_i, AP_o)])
    cf <- IndividualDataPP$data_information$conversion_factor

    # expand each (AP_i, AP_o) across all DP_rev_o = 1..max_DP and compute ranges
    max_DP <- pkg.env$maximum.time(IndividualDataPP$data_information$years,
                                   IndividualDataPP$data_information$output_time_granularity)

    dp_ranges <- development_periods[
      , .(DP_rev_o = 1:max_DP), by = .(AP_i, AP_o)  # expand rows per group
    ][
      , `:=`(
        min_dp = AP_i + (DP_rev_o - AP_o) / cf,
        max_dp = AP_i - 1 + (DP_rev_o - AP_o + 1) / cf
      )
    ]


    data_information$dp_ranges <-dp_ranges
    #################################################################################################################

  out=list(model.out=list(data=X,
                          model.out=model.out),
           hazard_frame = hazard_frame,
           data_information = data_information,
           fit_information = list(hazard_model = hazard_model,
                                  is_lkh=is_lkh,
                                  os_lkh=os_lkh))

  class(out) <- c('ReSurvFit')

  return(out)
}



