#' ReSurv fit predictions
#'
#' Predict the reserve for a ReSurv fit
#'
#' @param object resurv model to use for forecasting
#' @param IndividualData individual data to preprocess for the fit.
#'
#' @return reserve predictions.
#'
#' @export
predict.ReSurvFit <- function(object,
                              IndividualData) {


  data= IndividualData$training.data
  #Same as above just for quarterly output.

  trunc_q <- ceiling(max(data$TR_i)/3)

  data_reserve <- data %>%
    filter(DP_rev_i>TR_i) %>%
    mutate(AP_o = ceiling(AP_i/3),
           DP_o = ceiling(DP_i/3+((AP_i-1)%%3)/3 ))  %>%
    mutate(across(IndividualData$categorical_features, as.factor))


  max_DP_o <- data_reserve %>% group_by(AP_o) %>%
    summarise(DP_max = trunc_q-AP_o+1) %>%
    distinct()

  data_reserve2 <- data_reserve %>%
    filter(DP_rev_i>TR_i) %>%
    #group_by_at(covariates) %>%
    select(AP_o, IndividualData$categorical_features,I) %>%
    mutate(AP_o = as.numeric(AP_o)) %>%
    left_join(max_DP_o, by="AP_o")

  data_all <- data %>%
    mutate(AP_o = ceiling(AP_i/3),
           DP_o = ceiling(DP_i/3+((AP_i-1)%%3)/3 ))  %>%
    mutate(across(IndividualData$categorical_features, as.factor))

  list_cv <- list()
  list_IBNR <- list()


  tmp_cl_ibnr <- data_reserve2 %>%  group_by(AP_o, DP_max) %>%
    summarise(I=sum(I))


  tmp_cv_ibnr <- data_reserve2 %>%  group_by(AP_o, !!sym(IndividualData$categorical_features), DP_max) %>%
    summarise(I=sum(I)) %>%  reshape2::dcast(paste0("AP_o + DP_max~",paste(IndividualData$categorical_features, collapse='+') ), value.var = 'I')


  X_tmp<-matrix(ncol=0,nrow=nrow(data_reserve))
  for(i in IndividualData$categorical_features){
    string_formula <- paste("survival::Surv(TR_i, DP_rev_i, I) ~ ",paste(i, collapse='+'))
    tmp<-model.matrix(terms(as.formula(paste(string_formula,"+0"))), model.frame(as.formula(paste(string_formula,"+0")), data=data_reserve))
    X_tmp<-cbind(X_tmp,tmp)
  }
  beta_logicals_NT<-list()

  for(k in 1:ncol(unique(X_tmp))){
    beta_logicals_NT[[k]] <- rowSums(X_tmp == unique(X_tmp)[k,][col(X_tmp)]) == ncol(X_tmp) #NoTime covariates
  }

  names_hazard <- (data.frame(unique(X_tmp)) %>%
                     rowwise() %>%
                     mutate(name = paste0(names(.)[c_across() == 1], collapse = ',')))$name


  for (i in 1:(length(beta_logicals_NT)+1) ){

    I <- matrix(nrow=nrow(tmp_cl_ibnr),  ncol=1)
    E <- matrix(nrow=nrow(tmp_cl_ibnr), ncol=nrow(tmp_cl_ibnr))
    for(j in 2:nrow(E)){

      data_all_2 <- data_all %>%
        #group_by_at(covariates) %>%
        select(AP_o, DP_o, IndividualData$categorical_features,I) %>%
        mutate(AP_o = as.numeric(AP_o)) %>%
        filter(DP_o<=(j-1)) %>%
        select(AP_o, IndividualData$categorical_features,I)

      tmp_cl <- data_all_2 %>%  group_by(AP_o) %>%
        dplyr::summarise(I=sum(I))

      tmp_cv <- data_all_2 %>%  group_by(AP_o, !!sym(IndividualData$categorical_features)) %>%
        dplyr::summarise(I=sum(I)) %>%  reshape2::dcast(paste0("AP_o~",paste(IndividualData$categorical_features, collapse='+') ), value.var = 'I')

      if(i==1){ #cl
        df_E <- object$df[j-1,i][[1]]
        E[(nrow(E)-j+2):nrow(E),j] <- as.matrix(tmp_cl[(nrow(E)-j+2):nrow(E),2]) * df_E-as.matrix(tmp_cl[(nrow(E)-j+2):nrow(E),2])

        max<-tmp_cl_ibnr[j,2][[1]]
        df_E <- prod(object$df[max:(nrow(object$df)),i])
        I[j] <- tmp_cl_ibnr[j,3][[1]] * df_E

      }
      else{
        for(k in (nrow(E)-j+2):nrow(E)){

          if((ncol(tmp_cv)) >=i){ #if no claims observed yet
            df_E <- object$df[j-1,(k+1)+(i-2)*16][[1]] #16 due to simulated structure
            E[k,j] <- sum(tmp_cv[k,i][[1]],na.rm=T) * df_E-sum(tmp_cv[k,i][[1]],na.rm=T)
            if( j>2 & j<8 & k>nrow(E)-j+2 & sum(E[k,j-1],na.rm=T)==0  ){
              E[k,j-1]<-mean(as.matrix(tmp_cv[1:(nrow(E)-j+1),i]), na.rm=T) #if we have no exposure in the first couple of months. Corresponds to a capecod approximation
            }
          }
          E[(nrow(E)-j+2):nrow(E),j][is.na(E[(nrow(E)-j+2):nrow(E),j])]<-0
        }
        max<-tmp_cv_ibnr[j,2][[1]]
        df_E <- prod(object$df[max:(nrow(object$df)),(j+1)+(i-2)*16])
        I[j] <- sum(tmp_cv_ibnr[j,1+i][[1]],na.rm=T) * df_E

      }
    }
    list_cv[[i]]<- E
    list_IBNR[[i]] <- I
  }


  list_cv[[length(list_cv)+1]]<-Reduce('+', list_cv[2:length(list_cv)])
  names(list_cv)<-c('CL',names_hazard, 'CV_combined')
  list_IBNR[[length(list_IBNR)+1]]<-Reduce('+', list_IBNR[2:length(list_IBNR)])
  names(list_IBNR)<-c('CL',names_hazard, 'CV_combined')

  return(list(list_cv=list_cv, list_IBNR = list_IBNR))

}
