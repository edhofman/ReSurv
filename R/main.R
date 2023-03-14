individual_data <- ReSurv::IndividualData(simulated_dataframe_RM_CT,
                                  id="claim_number",
                                  continuous_features=NULL,
                                  categorical_features="claim_type",
                                  accident_period="AM",
                                  calendar_period="RM",
                                  input_time_unit=1/12,
                                  output_time_unit=1/4,
                                  continuous_features_spline=TRUE,
                                  degrees_of_freedom=4)

out <- ReSurv::ReSurv(individual_data)
