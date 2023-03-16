# Vector of variables

vector_of_variables <- c("DP_i",
                         "AP_i",
                         "DP_rev_i",
                         "TR_i",
                         "AP_o",
                         "DP_rev_o",
                         "TR_o",
                         "years",
                         "time_unit",
                         "p_month",
                         "time_w",
                         "weight")

globalVariables(vector_of_variables)

Sys.setenv('_R_CHECK_SYSTEM_CLOCK_' = 0)


