# Hazard computation helper functions
#
# Hazard frame, covariate mapping, predictions, and development factors.
#
# @import data.table
## Hazard computation ----

hazard_f<-function(i,
                           enter,
                           time,
                           exit,
                           event){
  "
  This function computes the hazard over a matrix.

  "
  rs <- (enter < time[i] & exit >= time[i])

  O<- sum(event[exit == time[i]])
  E <- sum(rs,-1/2*event[exit == time[i]] )

  c(O/E)}


hazard_data_frame <- function(hazard,
                                      # Om.df,
                                      eta = 0.5,
                                      eta_old = NULL,
                                      categorical_features,
                                      continuous_features,
                                      calendar_period_extrapolation){
  "
  Convert hazard matrix to dataframe and add grouping variables.

  "

  if (!is.null(eta_old)) {
    eta <- eta_old
  }
  eta <- validate_eta(eta)
  continuous_features_group=unique(c("AP_i",continuous_features))

  #Calculate input development factors and corresponding survival probabilities
  hazard_frame_tmp <- hazard %>%
    # dplyr::left_join(Om.df, "DP_rev_i") %>%
    dplyr::mutate(dev_f_i = (1+(1-eta)*hazard)/(1-eta*hazard) ) %>% #Follows from the assumption that claims are distributed evenly in the input period
    # dplyr::mutate(dev_f_i = (2*Om+(Om+1)*hazard)/(2*Om-(Om-1)*hazard) ) %>%
    tidyr::replace_na(list(dev_f_i =1)) %>%
    dplyr::mutate(dev_f_i = ifelse(dev_f_i < 0, 1, dev_f_i))

  hazard_frame_tmp <- hazard_frame_tmp %>%
    dplyr::group_by(pick(dplyr::all_of(c(categorical_features, continuous_features_group)))) %>%
    dplyr::arrange(DP_rev_i) %>%
    dplyr::mutate(cum_dev_f_i = cumprod(dev_f_i)) %>%
    dplyr::mutate(S_i = ifelse(cum_dev_f_i==0,0,1/cum_dev_f_i), # to handle the ifelse statement from above
           S_i_lead = dplyr::lead(S_i, default = 0),
           S_i_lag = dplyr::lag(S_i, default = 1)) %>%
    dplyr::select(-c(expg, baseline, hazard))

  # continuous_features_group=unique(c("AP_i",continuous_features))

  hazard_frame <- hazard %>%
    dplyr::left_join(hazard_frame_tmp, c(categorical_features,
                                  continuous_features_group,
                                  "DP_rev_i")) %>%
    dplyr::mutate(dev_f_i = dplyr::coalesce(dev_f_i,1),
           S_i = dplyr::coalesce(S_i,1),
           S_i_lead = dplyr::coalesce(S_i_lead,1),
           S_i_lag = dplyr::coalesce(S_i_lag, 1),
           cum_dev_f_i = dplyr::coalesce(cum_dev_f_i,1))
  return(hazard_frame)
  }


covariate_mapping <- function(hazard_frame,
                                      categorical_features,
                                      continuous_features,
                                      conversion_factor,
                                      calendar_period_extrapolation)
{
  "
  Create a dimension table, that holds a link between inputted categorical features and the group, that is used for expected_values
  "

  #
  #Need to handle Accident/calender period effect seperatly
  if( (length(continuous_features)==1 & "AP_i" %in% continuous_features) |
      (length(continuous_features)==1 & "RP_i" %in% continuous_features) |
      (length(continuous_features)==2 & sum(c("AP_i","RP_i") %in% continuous_features))==2 ){
    continuous_features_group = NULL
  }else{
    continuous_features_group=continuous_features[!(continuous_features %in% c("AP_i","RP_i"))]
  }

  ## Generate a grouping key, used for aggregating from input periods to output periods


  setDT(hazard_frame)
  feature_cols <- unique(c(categorical_features, continuous_features_group))

  # Create feature.id efficiently
  if (is.null(feature_cols) || length(feature_cols) == 0L) {
    hazard_frame[, covariate := "0"]
  } else {
    hazard_frame[, covariate := do.call(paste, c(.SD, sep = "_")), .SDcols = feature_cols]
  }

  # hazard_frame$covariate <- name_covariates(
  #   hazard_frame,
  #   categorical_features,
  #   continuous_features_group
  # )


  #Group is pr. covariate, output accident period
  #Ongoing update to be able to handle calender-period
  if("AP_i" %in% continuous_features |
     "RP_i" %in% continuous_features){

    #Generic approach that groups by either AP_i, RP_i or both, and create the corresponding dataset.
    time_features <- continuous_features[continuous_features %in% c("AP_i","RP_i")]

    time_elements_0 <- paste(sapply(time_features, function(x){paste0(x,"=hazard_frame[['",x,"']]")}
    ), collapse=", ")
    time_elements_1 <- paste(sapply(time_features, function(x){paste0("'",x,"'")}
    ), collapse=", ")

    # expression_0 <- paste0(sprintf(
    #   "groups <- data.frame(%s, covariate = hazard_frame$covariate)",
    #   time_elements_0    ),
    #   " %>%distinct()%>%   dplyr::mutate(group_i = row_number())")


    expression_0 <- sprintf(
      "groups <- data.table(%s, covariate = hazard_frame$covariate)",
      time_elements_0    )

    eval(parse(text=expression_0))

    # Keep only unique rows and add row numbers
    groups <- unique(groups)                # distinct rows
    groups[, group_i := .I]


    setDT(hazard_frame)


    expression_1 <-(sprintf("hazard_group <-merge(
      hazard_frame,
      groups,
      by = c(%s, 'covariate'),
      all.x = TRUE
    )",time_elements_1))

    # expression_1 <- paste0(
    #   "hazard_group <- hazard_frame %>%  dplyr::left_join(groups, by=",
    #   sprintf(
    #     "c(%s, 'covariate'))",
    #     time_elements_1    ) )


    eval(parse(text=expression_1))


  }else{
    #Only by covariate, since no time dependency.


    setDT(hazard_frame)

    # Step 1: create unique groups with row numbers
    groups <- unique(hazard_frame[, .(covariate)])
    groups[, group_i := .I]  # .I gives row number

    # Step 2: left join hazard_frame with groups
    hazard_group <- merge(
      hazard_frame,
      groups,
      by = "covariate",
      all.x = TRUE
    )

    # groups <- unique(data.frame(covariate = hazard_frame$covariate)) %>%
    #   dplyr::mutate(group_i = row_number())
    #
    # hazard_group <- hazard_frame %>%  dplyr::left_join(groups, by=c("covariate"))
  }



  #If we have to group for later output, add the relevant groups as well
  groups$group_o <- groups$group_i
  # The only time the groups will be different, is when we are including accident period as a covariate
  # As the dimension of the time periods are changing, we add the group_o output to keep track of which output group each input period corresponds to.
  if(conversion_factor != 1 & sum(c("AP_i","RP_i") %in% continuous_features)>0 ){

    time_elements_0 <- paste(sapply(time_features, function(x){
      paste0(substr(x,1,2),"_o  =ceiling(hazard_group[['",x,"']]*conversion_factor)")}),
      collapse=", ")

    time_elements_1 <- paste(sapply(time_features, function(x){paste0("",substr(x,1,2),"_o=ceiling(",x,"*conversion_factor)")}
    ), collapse=", ")

    time_elements_2 <- paste(sapply(time_features, function(x){paste0("'",substr(x,1,2),"_o'")}
    ), collapse=", ")


    # expression_0 <- paste0(sprintf(
    #   "      groups_o <- data.frame(%s, covariate = hazard_group$covariate)",
    #   time_elements_0    ),
    #   "%>% distinct() %>% dplyr::mutate(group_o = row_number())")

    expression_0 <- paste0(sprintf(
      "      groups_o <- data.table(%s, covariate = hazard_group$covariate)",
      time_elements_0    ))

    eval(parse(text=expression_0))

    groups_o <- unique(groups_o)

    # Add row numbers
    groups_o[, group_o := .I]

    groups[,group_o:=NULL]

    groups[,(paste0(substr(time_features,1,2),"_o")):=lapply(.SD, function(x) ceiling(x*conversion_factor)),.SDcols = time_features]


    groups <- merge(
      groups,
      groups_o,
      by = c(paste0(substr(time_features,1,2),"_o"), "covariate"),
      all.x = TRUE
    )

    # expression_1 <- paste0(
    #   "groups <- groups %>% dplyr::select(-group_o) %>%",
    #   sprintf(
    #     " dplyr::mutate(%s)",
    #     time_elements_1    ),
    #   " %>% ",
    #   sprintf(
    #     " dplyr::left_join(groups_o, by=c(%s, 'covariate'))",
    #     time_elements_2    ) )



    # eval(parse(text=expression_1))

  }

  return(list(hazard_group=hazard_group, groups = groups))



}


latest_observed_values_i <- function(data_reserve,
                                             groups,
                                             categorical_features,
                                             continuous_features,
                                             calendar_period_extrapolation){
  "
  Retrieve total amount of observed claims as of today.

  "
  #Max possible development time per accident period
  max_DP_i <- data_reserve %>% dplyr::group_by(AP_i) %>%
    dplyr::summarise(DP_max_rev =min(max(DP_rev_i)-DP_i)+1 ) %>%
    distinct()

  data_reserve2 <- data_reserve %>%
    dplyr::select(AP_i, AP_o, DP_rev_i, DP_i, dplyr::all_of(categorical_features), dplyr::all_of(continuous_features), I) %>%
    dplyr::mutate(AP_i = as.numeric(AP_i)) %>%
    dplyr::left_join(max_DP_i, by="AP_i")

  #The reason for the if statement is due to the !!sym logic, because !!sym(NULL) is not valid
  if(is.null(continuous_features)){ #length(continuous_features) == 1 & "AP_i" %in% continuous_features

    #latest observed pr. covariates
    observed_so_far <- data_reserve2 %>%  dplyr::group_by(pick(dplyr::all_of(categorical_features), AP_i, AP_o, DP_max_rev )) %>%
      dplyr::summarise(latest_I=sum(I), .groups = "drop")

    # observed pr. development period
    observed_dp_rev_i <- data_reserve2 %>%  dplyr::group_by(pick(AP_i, AP_o, dplyr::all_of(categorical_features),
                                                          DP_rev_i, DP_i)) %>%
      dplyr::summarise(I=sum(I), .groups = "drop")

    #Combine covariate values into single variable and add group dimension
    observed_so_far$covariate <- name_covariates(
      observed_so_far,
      categorical_features,
      continuous_features
    )
    observed_dp_rev_i$covariate <- name_covariates(
      observed_dp_rev_i,
      categorical_features,
      continuous_features
    )

    # Latest cumulative
    observed_so_far_out <- observed_so_far %>%  dplyr::left_join(groups, by=c("covariate")) %>%
      dplyr::select(AP_i, group_i, DP_max_rev, latest_I)

    observed_dp_rev_i_tmp <- observed_dp_rev_i %>%  dplyr::left_join(groups, by=c("covariate")) %>%
      dplyr::select(AP_i, group_i, DP_rev_i, DP_i, I)



  } else{
    #Very similiar to above, expect we now also take the continuous features into account

    if(( (length(continuous_features)==1 & "AP_i" %in% continuous_features) |
         (length(continuous_features)==1 & "RP_i" %in% continuous_features) |
         (length(continuous_features)==2 & sum(c("AP_i","RP_i") %in% continuous_features))==2 )){
      continuous_features_group<-NULL
    }
    else{
      continuous_features_group <- continuous_features[!(continuous_features %in% c("AP_i","RP_i") )]
    }

    #If-statment due to grouping by continous variable
    #handle <- "RP_i" %in% continuous_features
    #For now we let the handle be false, this is part of an on-going calendar-period implementation
    handle = 2
    if(is.null(continuous_features_group)){
      handle = 2
      observed_so_far <-
        switch(handle,
               data_reserve2 %>%  dplyr::group_by(pick(AP_i, AP_o,  RP_i, dplyr::all_of(categorical_features),
                                                DP_max_rev)) %>%
                 dplyr::summarise(latest_I=sum(I), .groups = "drop"),
               data_reserve2 %>%  dplyr::group_by(pick(AP_i, AP_o, dplyr::all_of(categorical_features),
                                                DP_max_rev)) %>%
                 dplyr::summarise(latest_I=sum(I), .groups = "drop")
        )

      observed_dp_rev_i <-
        switch(handle,
               data_reserve2 %>%  dplyr::group_by(pick(AP_i, AP_o, RP_i, dplyr::all_of(categorical_features),
                                                DP_rev_i, DP_i)) %>%
                 dplyr::summarise(I=sum(I), .groups = "drop"),
               data_reserve2 %>%  dplyr::group_by(pick(AP_i, AP_o, dplyr::all_of(categorical_features),
                                                DP_rev_i, DP_i)) %>%
                 dplyr::summarise(I=sum(I), .groups = "drop")
        )
    }
    else{


      observed_so_far <- switch(handle,
                                data_reserve2 %>%  dplyr::group_by(pick(AP_i, AP_o, RP_i, dplyr::all_of(categorical_features),
                                                                 dplyr::all_of(continuous_features_group)),
                                                            DP_max_rev) %>%
                                  dplyr::summarise(latest_I=sum(I), .groups = "drop"),
                                data_reserve2 %>%  dplyr::group_by(pick(AP_i, AP_o, dplyr::all_of(categorical_features),
                                                                 dplyr::all_of(continuous_features_group),
                                                                 DP_max_rev)) %>%
                                  dplyr::summarise(latest_I=sum(I), .groups = "drop")
      )

      observed_dp_rev_i <- switch(handle,
                                  data_reserve2 %>%  dplyr::group_by(pick(AP_i, AP_o, RP_i, dplyr::all_of(categorical_features),
                                                                   dplyr::all_of(continuous_features_group),
                                                                   DP_rev_i, DP_i)) %>%
                                    dplyr::summarise(I=sum(I), .groups = "drop"),
                                  data_reserve2 %>%  dplyr::group_by(pick(AP_i, AP_o, dplyr::all_of(categorical_features),
                                                                   dplyr::all_of(continuous_features_group),
                                                                   DP_rev_i, DP_i)) %>%
                                    dplyr::summarise(I=sum(I), .groups = "drop")
      )
    }

    observed_so_far$covariate <- name_covariates(
      observed_so_far,
      categorical_features,
      continuous_features_group
    )

    observed_dp_rev_i$covariate <- name_covariates(
      observed_dp_rev_i,
      categorical_features,
      continuous_features_group
    )

    time_features <- continuous_features[continuous_features %in% c("AP_i","RP_i")]


    observed_so_far_out <- observed_so_far %>%  dplyr::left_join(groups, by=c(time_features, "covariate")) %>%
      dplyr::select(AP_i, dplyr::all_of(time_features), group_i, DP_max_rev,latest_I )

    observed_dp_rev_i_tmp <- observed_dp_rev_i %>%  dplyr::left_join(groups, by=c(time_features, "covariate")) %>%
      dplyr::select(AP_i, dplyr::all_of(time_features), group_i, DP_rev_i, DP_i, I) %>%
      dplyr::inner_join(groups[,c(time_features, "group_i"), drop = FALSE], by =c(time_features, "group_i")) #filter only relevant combinations


  }

  return(list(latest_cumulative = observed_so_far_out, observed_pr_dp = observed_dp_rev_i_tmp))

}

name_covariates <-function(data, categorical_features, continuous_features){

  #
  feats <- c(categorical_features,continuous_features)

  if(is.null(feats)){return("0")}

  model_features <- data %>%
    dplyr::select(dplyr::all_of(feats)) %>%
    as.data.frame()

  mylist<- mapply(paste,
                  colnames(model_features),
                  lapply(model_features,c),
                  MoreArgs= list(sep="_"))

  mydf <- as.data.table(mylist)

  # mydf[,features.id:=paste(.SD,collapse=",")]
  #
  mydf[, features.id:=do.call(paste0,.SD)]

  # features.id <- apply(mydf, MARGIN=1, paste, collapse=",")
  return(mydf$features.id)
}




predict_i <- function(hazard_data_frame,
                              latest_cumulative,
                              grouping_method,
                              min_DP_rev_i
){
  "
  Calculate expected incremental claim number on input scale.
  Grouping is used when doing granularity increased development factors in i_to_o_development_factor
   "

  # #select relevant hazard values
  grouped_hazard_0 <- hazard_data_frame %>% #for the last development, if we included group '0', we would be extrapolating for half a parallelogram - doesn't make sense
    dplyr::left_join(latest_cumulative, by=c("group_i", "AP_i"))

  # Predict expected numbers, this is also used grouping methodology
  # For probabilty assumed ultimate = 1, otherwise calculate ultiamte.
  expected <-  grouped_hazard_0 %>%
    dplyr::select(DP_rev_i, AP_i, group_i, S_i, S_i_lag, DP_max_rev, latest_I ) %>%
    dplyr::mutate(gm = grouping_method) %>%
    dplyr::left_join(hazard_data_frame %>%
                dplyr::mutate(DP_rev_i = DP_rev_i +1) %>%
                dplyr::select(DP_rev_i, AP_i, group_i, S_i) %>%
                dplyr::rename(S_ultimate_i = S_i), by=c("DP_max_rev"="DP_rev_i",
                                                 "AP_i" = "AP_i",
                                                 "group_i" = "group_i")) %>%
    dplyr::mutate(U=dplyr::case_when(
      gm == "probability" ~ 1,
      S_i_lag == 1 ~ latest_I,
      DP_max_rev == min_DP_rev_i ~ latest_I,
      S_ultimate_i ==0 ~ 0,
      AP_i != 1 ~ 1/S_ultimate_i * latest_I,
      TRUE ~ latest_I)) %>%
    dplyr::mutate(I_expected = U*(S_i_lag-S_i)) %>%
    dplyr::mutate(IBNR = ifelse(DP_rev_i < DP_max_rev, I_expected, as.numeric(NA)) ) %>%
    dplyr::select(AP_i, group_i, DP_rev_i, I_expected, IBNR) %>%
    as.data.frame()

  return(expected)


}

retrieve_df_i <- function(hazard_data_frame,
                                  groups,
                                  adjusted=FALSE,
                                  is_baseline_model=FALSE
){
  "
  Return data frame only containing input development factors.
  "
  if(!adjusted){
    df_i <- hazard_data_frame %>%
      dplyr::select(group_i, DP_rev_i, dev_f_i) %>%
      distinct() %>%
      as.data.table() %>% dcast(DP_rev_i ~group_i, value.var="dev_f_i") %>%
      dplyr::select(-DP_rev_i)
  }else{
    df_i <- hazard_data_frame %>%
      dplyr::select(group_i, DP_rev_i, df_i_adjusted) %>%
      distinct() %>%
      as.data.table() %>% dcast(DP_rev_i ~group_i, value.var="df_i_adjusted") %>%
      dplyr::select(-DP_rev_i)
  }

  #We only have 5 columns in the case of AP being included as covariate
  if(ncol(groups) == 5){
    colnames(df_i) <- c(paste0("AP_i_",groups$AP_i,",", groups$covariate ))
  }else{
    colnames(df_i) <- c(groups$covariate )
  }
  #

  if(is_baseline_model){

    df_i <- df_i %>%
      map_df(rev) %>%
      dplyr::mutate(DP_i=row_number())

    return(df_i)

  }else{
    df_i <- as.data.frame(df_i[1:(nrow(df_i)-1),]) %>%
    map_df(rev) %>%
    dplyr::mutate(DP_i=row_number())

    }

  return(df_i)


}


input_hazard_frame <- function(
    hazard_frame,
    expected_i,
    categorical_features,
    continuous_features,
    df_i,
    groups,
    adjusted=FALSE,
    is_baseline_model=FALSE)
{
  "
  Create a hazard frame with relevant input time granularity specific values for later output.
  "



  if("AP_i" %in% continuous_features & length(continuous_features) == 1){
    continuous_features <- NULL
  }
  else{
    continuous_features <- continuous_features[!("AP_i" %in% continuous_features)]
  }
  hazard_frame_input_relevant <- hazard_frame %>%
    dplyr::select(- c(cum_dev_f_i, S_i, S_i_lead, S_i_lag, covariate))

  #If AP is included as a grouping variable
  if(ncol(groups)==5){
    df_i_long <- df_i %>%
      as.data.table() %>% melt(id.vars="DP_i") %>%
      dplyr::left_join(groups %>%
                  dplyr::mutate(covariate = paste0("AP_i_", AP_i, ",", covariate) ) %>%
                  dplyr::select(c(covariate, group_i)), by=c("variable" = "covariate")) %>%
      dplyr::mutate(DP_i = DP_i + 1)

    colnames(df_i_long) <- c("DP_i", "covariate", "df_i", "group_i")

  }
  else{

    if(is_baseline_model){

      df_i_long <- df_i %>%
        as.data.table() %>% melt(id.vars="DP_i") %>%
        dplyr::mutate(variable=0) %>%
        dplyr::left_join(groups[,c("covariate", "group_i")], by=c("variable" = "covariate")) %>%
        dplyr::mutate(DP_i = DP_i +1) #to get correct DP_i

       colnames(df_i_long) <- c("DP_i", "covariate", "df_i", "group_i")
    }else{

    df_i_long <- df_i %>%
      as.data.table() %>% melt(id.vars="DP_i") %>%
      dplyr::left_join(groups[,c("covariate", "group_i")], by=c("variable" = "covariate")) %>%
      dplyr::mutate(DP_i = DP_i +1) #to get correct DP_i

    colnames(df_i_long) <- c("DP_i", "covariate", "df_i", "group_i")

    }
  }


  max_DP_rev_i = max(expected_i$DP_rev_i)



  hazard_frame_input <- expected_i %>%
    dplyr::mutate(DP_i = max_DP_rev_i-DP_rev_i +1) %>%
    dplyr::left_join(hazard_frame_input_relevant, by =c("group_i", "AP_i", "DP_rev_i")) %>%
    dplyr::left_join(df_i_long[, c("DP_i", "group_i", "df_i")], by = c("DP_i", "group_i")) %>%
    tidyr::replace_na(list(df_i = 1))

  #Ordering
  if(adjusted == FALSE){
    hazard_frame_input <- hazard_frame_input[,c(categorical_features,
                                                continuous_features,
                                                "AP_i",
                                                "DP_rev_i",
                                                "expg",
                                                "baseline",
                                                "hazard",
                                                "df_i",
                                                "group_i",
                                                "I_expected",
                                                "IBNR")]
  }
  if(adjusted==TRUE){
    hazard_frame_input <- hazard_frame_input[,c(categorical_features,
                                                continuous_features,
                                                "AP_i",
                                                "DP_rev_i",
                                                "expg",
                                                "baseline",
                                                "hazard",
                                                "df_i",
                                                "df_i_adjusted",
                                                "group_i",
                                                "I_expected",
                                                "IBNR")]
  }
  return(hazard_frame_input)

}


predict_o <- function(
    expected_i,
    groups,
    conversion_factor,
    years,
    input_time_granularity
){
  "
  Calculate expected incremential claim number on output scale

   "
  max_dp_i <-maximum.time(years,input_time_granularity)
  # Predict expected numbers, this is also used grouping methodology
  expected <-  expected_i %>%
    dplyr::left_join(groups[,c("group_i", "group_o")], by =c("group_i")) %>%
    dplyr::mutate(DP_i =  max_dp_i-DP_rev_i + 1) %>%
    dplyr::mutate(AP_o = ceiling(AP_i*conversion_factor),
           DP_rev_o = ceiling(max_dp_i*conversion_factor)- ceiling((DP_i+(AP_i-1)%%(1/conversion_factor))*conversion_factor)+1) %>%
    #we can consider re-adding it in the future: dplyr::filter(DP_rev_o >0) %>% #since for DP_rev_o = 0, we are working with half a parallelogram in the end of the development time
    dplyr::group_by(AP_o, DP_rev_o, group_o) %>%
    summarize(I_expected = sum(I_expected,na.rm=TRUE),
              IBNR = sum(IBNR, na.rm=TRUE), .groups="drop") %>%
    dplyr::select(AP_o, group_o, DP_rev_o, I_expected, IBNR)

  return(expected)


}

i_to_o_development_factor <- function(hazard_data_frame,
                                              expected_i,
                                              dp_ranges,
                                              groups,
                                              observed_pr_dp,
                                              latest_cumulative,
                                              conversion_factor,
                                              grouping_method,
                                              min_DP_rev_i,
                                              years,
                                              input_time_granularity){
  "
  Group input development factor to output.

  "

  max_dp_i <-maximum.time(years,input_time_granularity)
  # Add output groupings to relevant frames
  hazard_data_frame <- as.data.table(hazard_data_frame) %>%
    dplyr::left_join(groups[,c("group_i", "group_o")], by =c("group_i"))

  observed_pr_dp_o  <- as.data.table(observed_pr_dp) %>%
    dplyr::left_join(groups[,c("group_i", "group_o")], by =c("group_i"))

  latest_cumulative_o <- as.data.table(latest_cumulative) %>%
    dplyr::left_join(groups[,c("group_i", "group_o")], by =c("group_i")) %>%
    dplyr::group_by(AP_i, group_o, DP_max_rev) %>%
    summarize(latest_I = sum(latest_I, na.rm=TRUE), .groups = "drop")


  #For probability approach to grouping method we assume equal exposure for each accident period
  if(grouping_method == "probability"){
    expected_i <-  predict_i(
      hazard_data_frame = hazard_data_frame,
      latest_cumulative = latest_cumulative,
      grouping_method = "probability",
      min_DP_rev_i = min_DP_rev_i
    ) %>%
      dplyr::left_join(groups[,c("group_i", "group_o")], by =c("group_i"))
  } else{
    expected_i <-  expected_i %>%
      dplyr::left_join(groups[,c("group_i", "group_o")], by =c("group_i"))

  }

  #

  # #select relevant hazard value group and add output variables, and other variables to help with grouping
  grouped_hazard_0 <- hazard_data_frame %>%
    dplyr::mutate(DP_i =  max_dp_i-DP_rev_i + 1) %>%
    dplyr::mutate( DP_rev_o = floor(max_dp_i*conversion_factor)-ceiling(DP_i*conversion_factor+((AP_i-1)%%(1/conversion_factor))*conversion_factor) +1) %>%
    dplyr::filter(DP_rev_o > 0) %>%  #for the last development, if we included group '0', we would be extrapolating for half a parallelogram - doesn't make sense
    dplyr::left_join(dp_ranges, by=c("AP_i", "DP_rev_o")) %>%
    dplyr::left_join(latest_cumulative_o, by=c("group_o", "AP_i")) %>%
    dplyr::left_join(observed_pr_dp_o, by=c("group_o", "AP_i", "DP_rev_i"))

  # Create cumulative observed to find exposure for each period
  cumulative_observed <- observed_pr_dp_o %>%
    dplyr::group_by(AP_i, group_o) %>%
    dplyr::arrange(DP_i) %>%
    dplyr::mutate(exposure = cumsum(ifelse(is.na(I),0,I) )) %>%
    dplyr::mutate(DP_rev_i = DP_rev_i -1) %>%  #as we want this as exposure we join by the previous development period
    dplyr::select(AP_i, group_o, DP_rev_i, exposure)

  exposures <- grouped_hazard_0 %>%
    dplyr::group_by(AP_i, DP_rev_o, group_o) %>%
    dplyr::filter(DP_rev_i == max(DP_rev_i)) %>%
    dplyr::left_join(cumulative_observed, by=c("AP_i", "group_o",
                                        "max_dp"="DP_rev_i"))

  #Where we do not have any observed correct exposure we extrapolate based on fitted hazard
  no_exposure <-  exposures %>%
    dplyr::select(DP_rev_i,  DP_rev_o, AP_i, group_o, S_i, DP_max_rev, latest_I ) %>%
    dplyr::mutate(gm = grouping_method) %>%
    dplyr::left_join(hazard_data_frame %>%
                dplyr::mutate(DP_rev_i = DP_rev_i +1) %>%
                dplyr::select(DP_rev_i, AP_i, group_o, S_i) %>%
                dplyr::rename(S_ultimate_i = S_i), by=c("DP_max_rev"="DP_rev_i",
                                                 "AP_i" = "AP_i",
                                                 "group_o" = "group_o")) %>%
    dplyr::mutate(U=ifelse(
      S_ultimate_i ==0, 0,
      1/S_ultimate_i * latest_I) ) %>% #handle special ultimate cases
    dplyr::mutate(U = ifelse(DP_max_rev ==min_DP_rev_i , latest_I, U))  %>%
    dplyr::mutate(U = ifelse(gm=="probability", 1 ,U)) %>%
    dplyr::mutate(U = ifelse(latest_I==0,0,U))%>%
    dplyr::mutate(exposure_expected = U*(S_i)) %>%  #in theory one could say U*S_i- ifelse(DP_max_rev==DP_rev_i-1, latest_I, U*S_i_lead ), but this might lead to negative expected as we are not sure latest equal the same as distribution estimate
    dplyr::select(AP_i, group_o, DP_rev_o, DP_rev_i, exposure_expected)

  #Take seen exposure if possible, otherwise extrapolated exposure
  exposures_combined <- exposures  %>%
    dplyr::mutate(gm = grouping_method) %>%
    dplyr::left_join(no_exposure, by  = c(   "AP_i",
                                      "DP_rev_o",
                                      "DP_rev_i",
                                      "group_o")) %>%
    dplyr::mutate(exposure_combined = ifelse(gm == "probability",
                                      dplyr::coalesce(exposure_expected,0),
                                      dplyr::coalesce(exposure, exposure_expected))
    )

  #Take seen observed if possible otherwise extrapolated observed
  grouped_hazard_1 <- grouped_hazard_0 %>%
    dplyr::mutate(gm = grouping_method) %>%
    dplyr::left_join(expected_i, by  = c("AP_i",
                                  "group_o",
                                  "DP_rev_i")) %>%
    dplyr::mutate(I_combined = ifelse(gm == "probability",
                               dplyr::coalesce(I_expected,0),
                               dplyr::coalesce(I, I_expected,0))
    )

  #group to output scale
  grouped_hazard_2 <- grouped_hazard_1 %>%
    dplyr::group_by(AP_i, DP_rev_o, group_o) %>%
    summarize(observed = sum(I_combined), .groups="drop") %>%
    dplyr::left_join(exposures_combined, by=c("AP_i", "group_o", "DP_rev_o"))%>%
    dplyr::mutate(observed=ifelse(latest_I==0,0,observed))


  output_dev_factor <- grouped_hazard_2 %>%
    dplyr::group_by(DP_rev_o, group_o) %>%
    dplyr::summarise(dev_f_o = ifelse(sum(exposure_combined)==0,
                               1,
                               (sum(observed)+  sum(exposure_combined))/sum(exposure_combined)),.groups="drop" ) %>%
    as.data.table() %>%
    dcast(DP_rev_o ~group_o, value.var="dev_f_o")




  return(output_dev_factor[,-c("DP_rev_o")])

}

output_hazard_frame <- function(
    hazard_frame_input,
    expected_o,
    categorical_features,
    continuous_features,
    df_o,
    groups,
    is_baseline_model=FALSE
    )
{
  "
  Create output hazard frame

  "
  if("AP_i" %in% continuous_features & length(continuous_features) == 1){
    continuous_features <- NULL
  }
  else{
    continuous_features <- continuous_features[!("AP_i" %in% continuous_features)]
  }

  #Relevant variables, we do not include hazard, baseline, expg as they currently only live on input-level
  hazard_frame_input_relevant <- hazard_frame_input %>%
    dplyr::select(dplyr::all_of(categorical_features), dplyr::all_of(continuous_features), group_i) %>%
    dplyr::left_join(groups[,c("group_i", "group_o")], by =c("group_i")) %>%
    dplyr::select(-c(group_i)) %>%
    distinct()

  #If AP is included as a grouping variable
  if(ncol(groups)==5){
    df_o_long <- df_o %>%
      as.data.table() %>% melt(id.vars="DP_o") %>%
      dplyr::left_join(groups[,c("AP_o","covariate", "group_o")] %>%
                  dplyr::mutate(covariate = paste0("AP_o_", AP_o, ",", covariate) ) %>%
                  distinct(), by=c("variable" = "covariate"))  %>%
      dplyr::mutate(DP_o = DP_o +1) %>%  #to get correct
      dplyr::select(DP_o, variable, value, group_o)

    colnames(df_o_long) <- c("DP_o", "covariate", "df_o", "group_o")

  }
  else{

    if(is_baseline_model){

      df_o_long <- df_o %>%
          as.data.table() %>% melt(id.vars="DP_o") %>%
          dplyr::mutate(variable=0) %>%
          dplyr::left_join(groups[,c("covariate", "group_i")], by=c("variable" = "covariate")) %>%
        dplyr::mutate(DP_o = DP_o +1) #to get correct DP_i

        colnames(df_o_long) <- c("DP_o", "covariate", "df_o", "group_o")


    }else{

    df_o_long <- df_o %>%
      as.data.table() %>% melt(id.vars="DP_o") %>%
      dplyr::left_join(groups[,c("covariate", "group_o")], by=c("variable" = "covariate")) %>%
      dplyr::mutate(DP_o = DP_o +1) #to get correct

    colnames(df_o_long) <- c("DP_o", "covariate", "df_o", "group_o")}
  }


  max_DP_rev_o = max(expected_o$DP_rev_o)

  hazard_frame_output <- expected_o %>%
    dplyr::mutate(DP_o = max_DP_rev_o-DP_rev_o +1) %>%
    dplyr::left_join(hazard_frame_input_relevant, by =c("group_o")) %>%
    dplyr::left_join(df_o_long[, c("DP_o", "group_o", "df_o")], by = c("DP_o", "group_o")) %>%
    tidyr::replace_na(list(df_o = 1))

  hazard_frame_output <- hazard_frame_output[,c(categorical_features,
                                                continuous_features,
                                                "AP_o",
                                                "DP_rev_o",
                                                "df_o",
                                                "group_o",
                                                "I_expected",
                                                "IBNR")]
  return(hazard_frame_output)

}


update_hazard_frame <- function(
    hazard_frame_input,
    hazard_frame_grouped,
    df_o,
    latest_observed_i,
    groups,
    conversion_factor,
    categorical_features,
    continuous_features,
    check_value,
    years,
    input_time_granularity
    ){
  max_dp_i <-maximum.time(years,input_time_granularity)
  #Periods where we exceed the check_value
  relevant <- hazard_frame_input %>%
    dplyr::filter(hazard > check_value & DP_rev_i < max(DP_rev_i)) %>%
    dplyr::mutate(AP_o = ceiling(AP_i*conversion_factor),
           DP_rev_o =   floor(max_dp_i*conversion_factor)-ceiling(DP_i*conversion_factor+((AP_i-1)%%(1/conversion_factor))*conversion_factor) +1)


  max_DP_rev_o = max(relevant$DP_rev_o)

  relevant <- relevant %>%
    dplyr::mutate(DP_o = max_DP_rev_o-DP_rev_o +1)  %>%
    dplyr::left_join(groups[,c("group_i", "group_o")], by =c("group_i"))


  #If AP is included as a grouping variable
  if(ncol(groups)==5){
    df_o_long <- df_o %>%
      as.data.table() %>% melt(id.vars="DP_o") %>%
      dplyr::left_join(groups[,c("AP_o","covariate", "group_o")] %>%
                  dplyr::mutate(covariate = paste0("AP_o_", AP_o, ",", covariate)) %>%
                  dplyr::select("group_o","covariate") %>%
                  distinct() , by=c("variable" = "covariate"))

    colnames(df_o_long) <- c("DP_o", "covariate", "df_o", "group_o")

  }
  else{

    df_o_long <- df_o %>%
      as.data.table() %>% melt(id.vars="DP_o") %>%
      dplyr::left_join(groups[,c("covariate", "group_o")], by=c("variable" = "covariate"))#to get correct

    colnames(df_o_long) <- c("DP_o", "covariate", "df_o", "group_o")
  }

  #Gets latest observed on output scale to predict new development
  observed_o <-  latest_observed_i %>%
    dplyr::left_join(groups[,c("group_i", "group_o")], by =c("group_i")) %>%
    dplyr::mutate(AP_o = ceiling(AP_i*conversion_factor),
           DP_rev_o =   floor(max_dp_i*conversion_factor)-ceiling(DP_i*conversion_factor+((AP_i-1)%%(1/conversion_factor))*conversion_factor) +1) %>%
    dplyr::filter(DP_rev_o >0) %>% #since for DP_rev_o = 0, we are working with half a parallelogram in the end of the development time
    dplyr::mutate(DP_o = max_DP_rev_o-DP_rev_o +1) %>%
    dplyr::group_by(AP_o,group_o) %>%
    summarize(latest_I = sum(I), DP_o_max = max(DP_o), .groups="drop") %>%
    dplyr::select(AP_o, group_o, latest_I, DP_o_max) %>%
    dplyr::mutate(DP_o_join = DP_o_max+1)

  #handle that we set these cases to 1, hence cant find exposure
  no_exposure <- latest_observed_i %>%  dplyr::group_by(group_i, DP_rev_i) %>%
    summarize(I_help=sum(I), .groups="drop") %>%
    dplyr::inner_join(relevant[relevant$hazard>check_value,c("DP_rev_i", "group_i")], by =c("DP_rev_i", "group_i"))

  df_o_long_relevant <- df_o_long %>%  dplyr::inner_join(distinct(relevant[,c("DP_o", "group_o")])
                                                  , by=c("DP_o", "group_o"))

  #Predict new level on input scale.
  predict_new <- observed_o[observed_o$group_o %in% df_o_long_relevant$group_o,] %>%
    dplyr::left_join(df_o_long, by=c("DP_o_max" = "DP_o", "group_o")) %>%
    dplyr::mutate(I_new = latest_I*df_o-latest_I) %>%
    dplyr::mutate(I_new = I_new / (1/conversion_factor)^2) #assuming equal distribution in lower granularity

  #if I_expected is zero it is because hazard>check_value, hence we draw from no_exposure help
  #Calculate new development factors, by saying (new_predict + exposure)/exposure
  if("AP_i" %in% continuous_features){
    new_df <- relevant %>%
      dplyr::left_join(predict_new[,c("group_o", "DP_o_join", "I_new")], by =c("group_o", "DP_o" =  "DP_o_join")) %>%
      dplyr::left_join(no_exposure, by=c("group_i", "DP_rev_i")) %>%
      dplyr::mutate(I_expected = ifelse(I_expected==0,I_help, I_expected)) %>%
      dplyr::mutate(df_i_adjusted = dplyr::case_when(df_i == 1 ~ (I_new + I_expected)/(I_expected),
                                       TRUE ~  (I_expected/(df_i-1) + I_new)/(I_expected/(df_i-1)) )
      ) %>%
      dplyr::mutate(IBNR = I_new,
             I_expected = I_new) %>%
      dplyr::select(AP_i, group_i, DP_rev_i, df_i_adjusted, IBNR, I_expected) %>%
      tidyr::replace_na(list(df_i_adjusted=1))
  }
  else{
    new_df <- relevant[!is.na(relevant$IBNR),] %>%
      dplyr::left_join(predict_new[,c("group_o", "DP_o_join", "I_new")], by =c("group_o", "DP_o" =  "DP_o_join")) %>%
      dplyr::left_join(no_exposure, by=c("group_i", "DP_rev_i")) %>%
      dplyr::mutate(I_expected = ifelse(I_expected==0,I_help, I_expected)) %>%
      dplyr::mutate(df_i_adjusted = (I_expected/(df_i-1) + I_new)/(I_expected/(df_i-1)) ) %>%
      dplyr::mutate(IBNR = I_new,
             I_expected = I_new) %>%
      dplyr::select(AP_i, group_i, DP_rev_i, df_i_adjusted, IBNR, I_expected) %>%
      tidyr::replace_na(list(df_i_adjusted=1))
  }

  #Update the previous development factors where relevant.
  hazard_frame_grouped_2 <- hazard_frame_grouped %>%
    dplyr::mutate(df_i_adjusted = dev_f_i) %>%
    rows_update(new_df[,c("group_i","DP_rev_i", "df_i_adjusted")], by =c("group_i", "DP_rev_i")) %>%
    dplyr::group_by(pick(dplyr::all_of(categorical_features), AP_i)) %>%
    dplyr::arrange(DP_rev_i) %>%
    dplyr::mutate(cum_dev_f_i = cumprod(df_i_adjusted)) %>%
    dplyr::mutate(S_i = ifelse(cum_dev_f_i==0,0,1/cum_dev_f_i), # to handle the ifelse statement from above
           S_i_lead = dplyr::lead(S_i, default = 0),
           S_i_lag = dplyr::lag(S_i, default = 1)) %>%
    ungroup()

  return(hazard_frame_grouped_2)

}

