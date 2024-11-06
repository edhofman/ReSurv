# Vector of variables
# This script contains the global variables we define to ease the package computations.
vector_of_variables <- c("DP_i",
                         "AP_i",
                         "DP_rev_i",
                         "TR_i",
                         "AP_o",
                         "DP_rev_o",
                         "TR_o",
                         'RP_i',
                         'degree_cp',
                         "years",
                         "time_unit",
                         "p_month",
                         "time_w",
                         "weight",
                         "is_lkh",
                         "k",
                         "hazard",
                         "risks_s",
                         "events_s",
                         "ties",
                         "feature",
                         "value",
                         "IBNR",
                         "ix_group",
                         "S_i",
                         ".",
                         "id",
                         "expg",
                         "baseline",
                         "group_i",
                         "f_i",
                         "df_i",
                         "dev_f_i",
                         "group_o",
                         "cum_dev_f_i",
                         "cum_f_i",
                         "DP_o",
                         "f_o",
                         "df_o",
                         "I_expected")

globalVariables(vector_of_variables)

Sys.setenv('_R_CHECK_SYSTEM_CLOCK_' = 0)


