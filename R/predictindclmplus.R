#' Individual clmplus fit predictions
#'
#' Predict the reserve for an individual clmplus fit
#'
#' @param indclmplus resurv model to use for forecasting
#'
#' @return reserve predictions.
#'
#' @export
predict.indclmplus <- function(ReSurvFit,
                               ...) {
  trunc_m <- 48 #manually inputted due to simulated data structure, needs to be derived from data.

  #Get relevant part of data set
  data_reserve <- indclmplus$IndividualData$training.data %>%
    filter(DP_rev_i>TR_i) %>%
    mutate(across(indclmplus$IndividualData$categorical_features, as.factor))

  #Maximum observe development month
  max_dpi <- data_reserve %>%
    group_by(AP_i) %>%
    summarise(DP_max_i = trunc_m-AP_i+1) %>%
    distinct()


  data_reserve2 <- data_reserve %>%
    filter(DP_rev_i>TR_i) %>%
    #group_by_at(covariates) %>%
    select(AP_i, indclmplus$IndividualData$categorical_features,I) %>%
    mutate(AP_i = as.numeric(AP_i)) %>%
    left_join(max_dpi, by="AP_i")


  data_all <- indclmplus$IndividualData$training.data %>%
    mutate(across(indclmplus$IndividualData$categorical_features, as.factor))

  list_cv <- list()  #Create estimates for lower part of triangle pr. covariate
  list_IBNR <- list() #Create a estimate of total missing IBNR.

  #No covariates needed for CL.
  tmp_cl_ibnr <- data_reserve2 %>%
    group_by(AP_i, DP_max_i) %>%
    summarise(I=sum(I))

  #Create triangluar structure pr.covariate.
  tmp_cv_ibnr <- data_reserve2 %>%
    group_by(AP_i, !!sym(indclmplus$IndividualData$categorical_features),
             DP_max_i) %>%
    summarise(I=sum(I)) %>%
    reshape2::dcast(paste0("AP_i + DP_max_i~",
                           paste(indclmplus$IndividualData$categorical_features, collapse='+') ),
                    value.var = 'I')


  #Relevant development factors.
  df_factor <- indclmplus$df %>%
    filter(DP_i <= max(tmp_cv_ibnr$DP)) %>% #very weird. Check this
    as.matrix()

  #Mode matrix, that doesn't remove intercepts.
  X_tmp<-matrix(ncol=0,nrow=nrow(data_reserve))
  for(i in indclmplus$IndividualData$categorical_features){
    string_formula <- paste("survival::Surv(TR_i, DP_rev_i, I) ~ ",paste(i, collapse='+'))
    tmp<-model.matrix(terms(as.formula(paste(string_formula,"+0"))),
                      model.frame(as.formula(paste(string_formula,"+0")),
                                    data=data_reserve))
    X_tmp<-cbind(X_tmp,tmp)
  }
  #Unique beta without time ( time)
  beta_logicals_NT<-list()

  #I feel like this could have been done smart, tsince the beta_logicals are only used by its length.
  for(k in 1:ncol(unique(X_tmp))){
    beta_logicals_NT[[k]] <- rowSums(X_tmp == unique(X_tmp)[k,][col(X_tmp)]) == ncol(X_tmp) #NoTime covariates
  }

  names_hazard <- (data.frame(unique(X_tmp)) %>%
                     rowwise() %>%
                     mutate(names = paste0(names(.)[c_across() == 1], collapse = ',')))$names


  #Loop through all covariates-specific development factors. SInce time is naturally included when doing IBNR estimation, we remove from loop through betas.
  for (i in 1:(length(beta_logicals_NT)+1) ){

    E <- matrix(nrow=nrow(tmp_cl_ibnr), ncol=nrow(tmp_cl_ibnr))
    I <- matrix(nrow=nrow(tmp_cl_ibnr),  ncol=1)
    for(j in 2:nrow(E)){

      data_all_2 <- data_all %>%
        #group_by_at(covariates) %>%
        select(AP_i, DP_i, indclmplus$IndividualData$categorical_features,I) %>%
        mutate(AP_i = as.numeric(AP_i)) %>%
        filter(DP_i<=(j-1)) %>%
        select(AP_i, indclmplus$IndividualData$categorical_features,I)

      tmp_cl <- data_all_2 %>%  group_by(AP_i) %>%
        dplyr::summarise(I=sum(I))

      tmp_cv <- data_all_2 %>%  group_by(AP_i, !!sym(indclmplus$IndividualData$categorical_features)) %>%
        dplyr::summarise(I=sum(I)) %>%  reshape2::dcast(paste0("AP_i~",paste(indclmplus$IndividualData$categorical_features, collapse='+') ), value.var = 'I')

      if(i==1){ #cl
        df_E <- df_factor[j-1,i]
        E[(nrow(E)-j+2):nrow(E),j] <- as.matrix(tmp_cl[(nrow(E)-j+2):nrow(E),2]) * df_E-as.matrix(tmp_cl[(nrow(E)-j+2):nrow(E),2])

        max<-tmp_cl_ibnr[j,2][[1]]
        df_E <- prod(df_factor[max:(nrow(df_factor)),i])
        I[j] <- tmp_cl_ibnr[j,3][[1]] * df_E

      }
      else{
        #Again quite manual, and could most likely be updated with a generic setup
        for(k in (nrow(E)-j+2):nrow(E)){

          if((ncol(tmp_cv)) >=i){ #if no claims observed yet
            df_E <- df_factor[j-1,(k+1)+(i-2)*48]
            E[k,j] <- sum(tmp_cv[k,i][[1]],na.rm=T) * df_E-sum(tmp_cv[k,i][[1]],na.rm=T)
            if( j>2 & j<8 & k>nrow(E)-j+2 & sum(E[k,j-1],na.rm=T)==0  ){
              E[k,j-1]<-mean(as.matrix(tmp_cv[1:(nrow(E)-j+1),i]), na.rm=T) #if we have no exposure in the first couple of months. Corresponds to a capecod approximation
            }
          }
          E[(nrow(E)-j+2):nrow(E),j][is.na(E[(nrow(E)-j+2):nrow(E),j])]<-0
        }

        max<-tmp_cv_ibnr[j,2][[1]]
        df_E <- prod(df_factor[max:(nrow(df_factor)),j+(i-2)*48])
        I[j] <- sum(tmp_cv_ibnr[j,1+i][[1]],na.rm=T) * df_E

      }
    }
    list_cv[[i]]<- E
    list_IBNR[[i]]<- I
  }

  list_cv[[length(list_cv)+1]]<-Reduce('+', list_cv[2:length(list_cv)])
  names(list_cv)<-c('CL',names_hazard, 'CV_combined')
  list_IBNR[[length(list_IBNR)+1]]<-Reduce('+', list_IBNR[2:length(list_IBNR)])
  names(list_IBNR)<-c('CL',names_hazard, 'CV_combined')

  return(list(list_cv=list_cv, list_IBNR=list_IBNR))
}
