#' Survival continuously ranked probability score.
#'
#' Return the Survival Continuously Ranked Probability Score (SCRPS) of a \code{ReSurv} model.
#'
#' The model fit uses the theoretical framework of Hiabu et al. (2023), that relies on the
#'
#' @param ReSurvFit ReSurvFit object to use for the score computation.
#'
#' @return Survival CRPS.
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
survival_crps <- function(ReSurvFit){

  UseMethod("survival_crps")

}


#' Survival continuously ranked probability score.
#'
#' Return the Survival Continuously Ranked Probability Score (SCRPS) of a \code{ReSurv} model.
#'
#' The model fit uses the theoretical framework of Hiabu et al. (2023), that relies on the
#'
#' @param ReSurvFit ReSurvFit object to use for the score computation.
#'
#' @return Survival CRPS.
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
survival_crps.default <- function(ReSurvFit){

  message('The object provided must be of class ReSurvFit')

}

#' Survival continuously ranked probability score.
#'
#' Return the Survival Continuously Ranked Probability Score (SCRPS) of a \code{ReSurv} model.
#'
#' The model fit uses the theoretical framework of Hiabu et al. (2023), that relies on the
#'
#' @param ReSurvFit ReSurvFit object to use for the score computation.
#'
#' @return Survival CRPS.
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
survival_crps.ReSurvFit <- function(ReSurvFit){


  hazard_frame <- ReSurvFit$hazard_frame

  hazard_frame <- data.table(hazard_frame)

  categorical_features <- ReSurvFit$IndividualData$categorical_features
  continuous_features <- ReSurvFit$IndividualData$continuous_features

  hazard_frame <- hazard_frame[,
                               ix_group:=.GRP,
                               by=c(categorical_features,
                                    continuous_features)]

  hazard_list<- split(hazard_frame, hazard_frame$ix_group)

  hazard_list<-lapply(hazard_list, function(x) x[order(DP_rev_i),
                                                 .(DP_rev_i,
                                                   cdf2_i = (1-S_i)^2,
                                                   S2_i=S_i^2,
                                                   x.vals=c(1,diff(DP_rev_i)))])

  hazard_frame[['id']] = 1:dim(hazard_frame)[1]

  tmp <- hazard_frame[,.(crps=survival_information(x=DP_rev_i,
                                                   group=ix_group,
                                                   hazard_list = hazard_list)),by=id]


  return(tmp)


}







