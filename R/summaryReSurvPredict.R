#' Summary of IBNR predictions
#'
#' Gives overview of IBNR predictions
#'
#' @param object "ReSurvPredict" object specifying hazard and development factors.
#' @param granularity \code{character}, specify if which granularity the summary should be on.
#' \itemize{
#' \item{\code{"input"}}
#' \item{\code{"output"}}
#' }
#' Default is \code{"input"}.
#' @param ... Other arguments to be passed to summary.
#'
#' @return Summary of predictions
#'
#' @import data.table
#'
#' @export
#' @method summary ReSurvPredict

summary.ReSurvPredict <- function(object, granularity = "input", ...)
{
  handle <- match(granularity, c("input","output"))

  IBNR_AP <- switch(handle,
                    data.table(object$long_triangle_format_out$input_tg)[, .(IBNR=sum(IBNR, na.rm=T)), by = AP_i],
                    data.table(object$long_triangle_format_out$output_tg)[, .(IBNR=sum(IBNR, na.rm=T)), by = AP_o]
  )

  # development_factor = switch(handle,
  #                             object$df_input,
  #                             object$df_output)

  keep = c(
    "model.out",
    "IndividualDataPP",
    "hazard_model"
  )

  summary <- list(
    IBNR_AP = IBNR_AP,
    total_IBNR = sum(IBNR_AP$IBNR),
    # development_factor=development_factor,
    grouping_method = object$grouping_method,
    granularity = granularity,
    ReSurvFit = object$ReSurvFit[keep]
  )

  class(summary) <- "summaryReSurvPredict"
  return(summary)
}

#' Print summary of IBNR predictions
#'
#' Gives overview of IBNr predictions
#'
#' @param x "ReSurvPredict" object specifying hazard and development factors.
#' @param digits \code{numeric}, number of digits to print for IBNR and Likelihood.
#' @param ... Other arguments to be passed to print.
#'
#' @return print of summary of predictions
#'
#' @export
#' @method print summaryReSurvPredict

print.summaryReSurvPredict <-
  function (x, digits = max(3L, getOption("digits") - 3L), ...)
  {
    cat("\n Hazard model:\n",
        paste(deparse(x$ReSurvFit$hazard_model), sep = "\n", collapse = "\n"), "\n\n", sep = "")

    #cat("Likelihood: \n")
    #  xx <- x$ReSurvFit$model.out$likelihood
    #print.default(xx, digits = digits, na.print = "", print.gap = 2L)

    if(is.null(x$ReSurvFit$IndividualDataPP$categorical_features) & is.null(x$ReSurvFit$IndividualDataPP$continuous_features)) {
      cat("\nNo Features \n")
    } else {
      categorical_features<-NULL
      continuous_features <- NULL
      if(!is.null(x$ReSurvFit$IndividualDataPP$categorical_features)){
        categorical_features <- sprintf("\nCategorical Features:\n%s",
                                      paste(x$ReSurvFit$IndividualDataPP$categorical_features ,  collapse="\n"))
      }
      if(!is.null(x$ReSurvFit$IndividualDataPP$continuous_features)){
        continuous_features <- sprintf("\nContinuous Features:\n%s",
                                        paste(x$ReSurvFit$IndividualDataPP$continuous_features ,  collapse="\n"))
      }
      cat(categorical_features, continuous_features)

    }
    ##
    cat("\nTotal IBNR level: \n")
    print.default(x$total_IBNR, digits = digits, na.print = "", print.gap = 2L)

    # handle <- match(x$granularity, c("input","output"))
    #
    # df_string <- paste0("\n Development factors for ", x$granularity, " granularity.",
    #                     switch(
    #                       handle,
    #                       "",
    #                       paste0(" Estimated by ", x$ReSurvPredict$grouping_method, " priciple.")
    #                     ))
    #
    # cat(df_string
    # )
    # print(x$development_factor)


    invisible(x)
  }


#' Plot of the development factors
#'
#' Plots the development factors by group code.
#'
#' @param x "ReSurvPredict" object specifying hazard and development factors.
#' @param granularity \code{character}, either \code{"input"} for \code{input_time_granularity} or \code{"output"} for \code{output_time_granularity}.
#' @param group_code \code{numeric}: Identifier for the group that will be plotted. Default is 1. The code identifiers can be find in the \code{ReSurvPredict$long_triangle_format_out} list. Depending on the granularity of interest, it will be either in \code{ReSurvPredict$long_triangle_format_out$input_tg} for \code{input_time_granularity} or \code{ReSurvPredict$long_triangle_format_out$output_tg} for \code{output_time_granularity}.
#' @param color_par \code{character}: \code{ggplot2} Colour of the line plot. Default is \code{'royalblue'}. Optional.
#' @param linewidth_par \code{numeric}: Line plot width. Optional.
#' @param ylim_par \code{numeric}: Highest intercept on the y-axis (development factors). The default is the highest predicted development factor. Optional.
#' @param ticks_by_par \code{numeric}: gap between each x-axis label (development period). Default is 2. Optional.
#' @param base_size_par \code{numeric}: base size of the plot. Default is 5. See \code{base_size} in the \code{?theme_bw} documentation. Optional.
#' @param title_par \code{character}: Title of the plot. Optional.
#' @param x_text_par \code{character}: Text on the x-axis. Optional.
#' @param plot.title.size_par \code{numeric}: size of the plot title. Default is 20. See \code{size} in the \code{?element_text} documentation. Optional.
#' @param ... Other arguments to be passed to Plot. Optional.
#'
#' @return \code{ggplot2} of the development factors
#'
#' @export
#' @method plot ReSurvPredict

plot.ReSurvPredict <-function (x,
                               granularity = "input",
                               group_code=1,
                               color_par= "royalblue",
                               linewidth_par=2.5,
                               ylim_par=NULL,
                               ticks_by_par=NULL,
                               base_size_par=NULL,
                               title_par=NULL,
                               x_text_par=NULL,
                               plot.title.size_par=NULL,
                               ...){

  if(granularity=="input"){


    dtb_2_plot <- x$long_triangle_format_out$input_tg %>%
      filter(group_i==group_code,
             DP_i>1)

    if(is.null(ticks_by_par)){
      ticks.at <- seq(1,max(dtb_2_plot$DP_i),by=2)
    }else{

      ticks.at <- seq(1,max(dtb_2_plot$DP_i),by=ticks_by_par)

    }

    labels.as <- as.character(ticks.at)


    if(is.null(ylim_par)){
      ylim_setting <- max(dtb_2_plot$df_i)
    }else{

      ylim_setting <- ylim_par

    }

    if(is.null(x_text_par)){x_char = "DP_i"}else{x_char = x_text_par}


    ggplot_definition <- dtb_2_plot  %>%
      ggplot(aes(x=DP_i,
                 y=df_i),
             ...)



  }else{

    if(granularity=="output"){

      dtb_2_plot <- x$long_triangle_format_out$output_tg %>%
        filter(group_o==group_code,
               DP_o>1)

      if(is.null(ticks_by_par)){
        ticks.at <- seq(1,max(dtb_2_plot$DP_o),by=2)
      }else{

        ticks.at <- seq(1,max(dtb_2_plot$DP_o),by=ticks_by_par)

      }

      labels.as <- as.character(ticks.at)


      if(is.null(ylim_par)){
        ylim_setting <- max(dtb_2_plot$df_o)
      }else{

        ylim_setting <- ylim_par

      }

      if(is.null(x_text_par)){
        x_char = "DP_o"
      }else{
          x_char = x_text_par
          }

      ggplot_definition <- dtb_2_plot  %>%
        ggplot(aes(x=DP_o,
                   y=df_o),
               ...)




    }else{

      stop("granularity must be either 'input' or 'output'")

    }
  }


  if(is.null(base_size_par)){

    base_size_setting=rel(5)

  }else{

    base_size_setting=base_size_par
  }

  if(is.null(plot.title.size_par)){

    plot.title.size_setting=20

  }else{

    plot.title.size_setting=plot.title.size_par
  }

  if(is.null(title_par)){

    title_setting=paste("Development factors",granularity, "granularity", "group",group_code)

  }else{

    title_setting=title_par

  }


  p <- ggplot_definition +
    geom_line(linewidth=linewidth_par,
              color=color_par) +
    ylim(1,ylim_setting+.01)+
    labs(title=title_setting,
         x = x_char,
         y = "Development factor") +
    scale_x_continuous(breaks = ticks.at,
                       labels = labels.as) +
    theme_bw(base_size=base_size_setting)+
    theme(plot.title = element_text(size=plot.title.size_setting))

  return(p)










  }






