# library(ReSurv)
#
#
# ## Data simulation ----
# source('~/GitHub/individual_clmplus/R/package_requirements.R')
# require(SynthETIC)
# #library(tidyverse)
# set.seed(1964)
#
#
# # Letting Claimtype 1 have decreasing frequency ----------------------------------------------
# set_parameters(ref_claim = 200000, time_unit = 1/4)
# ref_claim = 200000
# time_unit <- 1/12
# years <- 4
# I <- years / time_unit
# E <- c(rep(6000, I))
# lambda <- c(rep(0.2, I))
# #
# #Decreasing the exposure, and hence lowering the claims occurred
# E_1 <- c(rep(6000, I)) + seq(from = 0, by = -100, length = I)
# #Frequency simulation
# n_vector_0 <- claim_frequency(I = I, E = E, freq = lambda)
# n_vector_1 <- claim_frequency(I = I, E = E_1*2, freq = lambda)
# occurrence_times_0 <- claim_occurrence(frequency_vector = n_vector_0)
# occurrence_times_1 <- claim_occurrence(frequency_vector = n_vector_1)
#
# #Claim sizes
# #Not needed for reporting analysis
# S_df <- function(s) {
#   # truncate and rescale
#   if (s < 30) {
#     return(0)
#   } else {
#     p_trun <- pnorm(s^0.2, 9.5, 3) - pnorm(30^0.2, 9.5, 3)
#     p_rescaled <- p_trun/(1 - pnorm(30^0.2, 9.5, 3))
#     return(p_rescaled)
#   }
# }
# claim_sizes_0 <- claim_size(frequency_vector = n_vector_0,
#                             simfun = S_df, type = "p", range = c(0, 1e24))
# claim_sizes_1 <- claim_size(frequency_vector = n_vector_1,
#                             simfun = S_df, type = "p", range = c(0, 1e24))
#
# RTFWD_inverse <- function(n, alpha, beta, lambda, k,b){
#   U<-runif(n)
#   (-log(U)/(beta^alpha*lambda^(alpha*k))+b^(-alpha*k))^(1/(-alpha*k))
# }
#
# notidel_param_0 <- function(claim_size, occurrence_period) {
#
#
#   #c(scale = 5/(exp(-25.7845)^(1/12.15343)),
#   #  shape = 12.15343) #cv 0.1
#   c(alpha=0.5,
#     beta=2,
#     lambda=0.1*exp(1.15129)^(1/0.5),
#     k=1,
#     b=years / time_unit)
#
# }
#
# notidel_param_1 <- function(claim_size, occurrence_period) {
#
#   # target_mean <- o 3/4 / time_unit
#   #
#   # # specify the target Weibull coefficient of variation
#   # target_cv <- 0.70
#   #
#   # c(shape = get_Weibull_parameters(target_mean, target_cv)[1, ],
#   #   scale = get_Weibull_parameters(target_mean, target_cv)[2, ])
#
#   #c(scale = 5/(exp(-22.6866)^(1/12.15343)),
#   #  shape = 12.15343)
#   c(alpha=0.5,
#     beta=2,
#     lambda=0.1*exp(1.95601)^(1/0.5),
#     k=1,
#     b=years / time_unit)
#
# }
#
#
# #sum(unlist(notidel_claim_type_0)>48)/sum(unlist(notidel_claim_type_0)>0)
# ## output
# # simulate notification delays from the transformed gamma
# notidel_claim_type_0 <- claim_notification(n_vector_0, claim_sizes_0,
#                                            rfun = RTFWD_inverse,
#                                            paramfun = notidel_param_0)
#
# notidel_claim_type_1 <- claim_notification(n_vector_1, claim_sizes_1,
#                                            rfun = RTFWD_inverse,
#                                            paramfun = notidel_param_1)
#
#
# # graphically compare the result with the default Weibull distribution
# plot(ecdf(unlist(notidel_claim_type_0)), xlim = c(0, 48),
#      main = "Empirical distribution of simulated notification delays",
#      xlab = "Notification delay (in months)")
# plot(ecdf(unlist(notidel_claim_type_1)), add = TRUE, col = 2)
# legend.text <- c("ClaimType0", "ClaimType1")
# legend("bottomright", legend.text, col = 1:2, lty = 1, bty = "n")
#
# ct0 <- tibble(AT = unlist(occurrence_times_0),
#               RT = unlist(occurrence_times_0) + unlist(notidel_claim_type_0),
#               claim_type = 0)
#
# ct1 <- tibble(AT = unlist(occurrence_times_1),
#               RT = unlist(occurrence_times_1) + unlist(notidel_claim_type_1),
#               claim_type = 1)
#
# simulated_dataframe_RM_CT <- ct0 %>%  bind_rows(ct1) %>%
#   mutate(
#     claim_number = row_number(),
#   )  %>%  mutate(
#     AM = ceiling(AT),
#     RM = ceiling(RT),
#     DT = RT-AT,
#     DM = RM-AM+1,
#     DM_rev = years/time_unit - DM+1,
#     DT_rev = years/time_unit - DT,
#     TR = AM-1, #just setting truncation to max year simulated. and accounting for
#     I=1
#   ) %>%
#   select(claim_number, AT, RT, claim_type, AM, RM, DT, DM, DM_rev, DT_rev, TR, I)
#
# # fwrite(simulated_dataframe_RM_CT,  file= "C:/Users/emil_/OneDrive/Dokumenter/Uni/Speciale/R code/Simulation/simulated_dataframe_1.csv")
#
# #

# #

#
#


# # #
#
# ## Locally, run deep surv ----
#
# data = individual_data$starting.data
# categorical_variables = individual_data$categorical_features
# numeric_variables = "AP_i"
# scaler = "minmax"
# traintestsplit = 0.8
# early_stopping = FALSE
# patience = 20
# network_structure = NULL
# num_nodes = c(10,10)
# activation = "LeakyReLU"
# lr = 0.005
# xi=1
# epsilon = 0
# tie="Efron"
#
#
# require(reticulate)
# require(data.table)
# require(dplyr)
#
#
# #Import python modules
# torchtuples <- reticulate::import("torchtuples")
# torch <- reticulate::import("torch")
#
# #Source python code for left truncated deepsurv
# reticulate::source_python("C:\\Users\\gpitt\\Documents\\GitHub\\individual_clmplus\\python\\coxnetwork_custom.py")
# source("C:\\Users\\gpitt\\Documents\\GitHub\\individual_clmplus\\R\\LTdeepsurv_helperfunctions.R")
#
# if(is.null(categorical_variables) & is.null(numeric_variables)){
#   warning("Neither categorical or numeric variables specified for neural network.")
#   break()
# }
#
# #subset covariates into x-dataframe
# names_x <- names(individual_data$training.data) %in% c(categorical_variables, numeric_variables)
#
# data_x <- individual_data$training.data[,names_x ]
#
# #determine the target variables. Should be duration, event, truncation.
# data_y <- individual_data$training.data[,c("DP_rev_i", "I", "TR_i")]
#
# #If numeric variables are specified, transform by selected scaler function
# if(!is.null(numeric_variables)){
#
#   #For numeric variables scale, if nothing specifed, we use min max.
#   if(toupper(scaler) == "MINMAX" ){
#     MinMaxScaler <- function(x, na.rm = TRUE) {
#       return((x- min(x)) /(max(x)-min(x)))
#     }
#     scaler <- MinMaxScaler
#   }
#
#   #Transform by selected scaler function
#   data_scaled <- mutate(data_x, across(all_of(numeric_variables), scaler )) %>%
#     select(all_of(numeric_variables))
# }else{
#   data_scaled <- data.frame(matrix(NA, nrow = nrow(data_x), ncol = 0))
# }
#
# #If we have categorical variables, we one hot encode them.
# if(!is.null(categorical_variables)){
#   #Create temporary ID and value to perform dcast one hot encoding for categorical variables
#
#   data_x$id <- 1:nrow(data_x)
#   data_x$value <- 1
#
#   #Loop thorugh all categorical variable to encode.
#   data_onehot <- data.frame(matrix(NA, nrow = nrow(data_x), ncol = 0))
#   for ( i in categorical_variables){
#     tmp <- reshape2::dcast(data = data_x, id ~ paste0(i,"_", data_x[[i]] ), fill = 0, value.var="value")
#     data_onehot <- cbind(data_onehot, tmp[2:ncol(tmp)]) #cbind from two and forward, as we don't want the id-column
#   }
#
# }else{
#   #If no categorical variables, create empty frame
#   data_onehot <- data.frame(matrix(NA, nrow = nrow(data_x), ncol = 0))
# }
#
# #Combine transformed variables to single frame
# data_transformed <- cbind(data_scaled, data_onehot)
#
# #Split into train and validation based upon given input
#
# #asset traintestsplit has value between 0 and 1
# if(traintestsplit>1 | traintestsplit<0 ){
#   warning(paste0("Traintestsplit has been put to ", traintestsplit,". The value needs to be between 0 and 1, defaulting to 0.8."))
#   traintestsplit <- 0.8
# }
#
# id_train <- sample(c(TRUE,FALSE), nrow(data_transformed), replace=T, prob= c(traintestsplit,1-traintestsplit) )
#
# #convert to array for later numpy transforamtion
# data_train <- as.array(as.matrix(data_transformed[id_train,]))
# y_train <- as.array(as.matrix(data_y[id_train,]))
# data_val <- as.array(as.matrix(data_transformed[!id_train,]))
# y_val <- as.array(as.matrix(data_y[!id_train,]))
#
#
# #create tuples holding target and validation values. Convert to same dtype to ensure safe pytorch handling.
# y_train <- tuple(np_array(y_train[,1], dtype = "float32"), #duration
#                  np_array(y_train[,2], dtype = "float32"), #event
#                  np_array(y_train[,3], dtype = "float32")) #truncation
#
# validation_data = tuple(np_array(data_val, dtype = "float32"),
#                         tuple(np_array(y_val[,1], dtype = "float32"), #duration
#                               np_array(y_val[,2], dtype = "float32"), #event
#                               np_array(y_val[,3], dtype = "float32"))) #truncation
#
# x_train = np_array(data_train, dtype = "float32")
#
#
#
# #create network structure, if correct class specified in network_strucutre use this, otherwise use other variables.
# if(!("torch.nn.modules.container.Sequential" %in% class(network_structure$net))){
#   net <- torch$nn$Sequential()
#   input_shape =  x_train$shape[[1]]
#   for( i in 1:length(num_nodes+1)){
#     if( i > length(num_nodes)){
#       # net$append(torch$nn$Linear(input_shape, as.integer(1)))
#       net$add_module('2',torch$nn$Linear(input_shape, as.integer(1)))
#     }
#     else{
#       net$add_module('0',torch$nn$Linear(input_shape,as.integer(num_nodes[i])))
#       net$add_module('1',torch$nn[[activation]]())
#       input_shape = as.integer(num_nodes[i])
#     }
#   }
# }
#
#
# #Setup batchsize, epochs and verbose settings
# batch_size=1000
# batch_size = as.integer(batch_size)
# epochs=100
# epochs = as.integer(epochs)
#
# verbose = T
#
#
# # Setup CoxPH model, as imported from python script.
# model <- CoxPH(
#   net = net,
#   optimizer = torchtuples$optim$Adam(lr = 0.005),
#   xi=xi,
#   eps=epsilon,
#   tie = tie
# )
#
# #If early stopping specified add to callbacks.
# if(early_stopping==TRUE){
#   callbacks = list(torchtuples$callbacks$EarlyStopping(patience=patience))
# }else{
#   callbacks = NULL
# }
#
#
# #fit model
# model$fit(input = x_train,
#   target = y_train,
#   batch_size = batch_size,
#   epochs = epochs,
#   callbacks = r_to_py(callbacks),
#   verbose = verbose,
#   val_data=validation_data,
#   val_batch_size=batch_size
# )













