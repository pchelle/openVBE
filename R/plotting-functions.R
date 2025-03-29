#' @export
getReferenceGoodnessOfFitPlot <- function(refAndTestSimulationsInVirtualPopulation,
                                          pkDataFilePath,
                                          referenceDrugName = "Reference",
                                          yDimension = NULL,
                                          yUnit = NULL,
                                          timeUnit = ospUnits$Time$min){

  plasmaDf <- refAndTestSimulationsInVirtualPopulation


  if (!is.null(yDimension) & !is.null(yUnit)){
    plasmaDf$conc <- toUnit(quantityOrDimension = yDimension,
                            values = plasmaDf$conc,
                            targetUnit = yUnit)
  }


  plasmaDf$time <- toUnit(quantityOrDimension = ospDimensions$Time,
                          values = plasmaDf$time,
                          targetUnit = timeUnit)

  mmdf <- NULL

  q05 <- function(x){quantile(x = x,0.05)}
  q95 <- function(x){quantile(x = x,0.95)}
  fnVec <- c("median" = median,"min" = min,"max" = max)
  for (fn in names(fnVec)){
    df <- aggregate(x = plasmaDf$conc,
                    by = list(plasmaDf$time,
                              plasmaDf$drug),
                    fnVec[[fn]])
    names(df) <- c("time","drug",fn)
    mmdf[["time"]] <- df[["time"]]
    mmdf[["drug"]] <- df[["drug"]]
    mmdf[[fn]] <- df[[fn]]
  }
  mmdf <- as.data.frame(mmdf)
  mmdf$drug <- as.character(mmdf$drug)
  mmdf$drug[mmdf$drug == "R"] <- referenceDrugName

  plt <- ggplot() + geom_ribbon(data = mmdf[mmdf$drug == referenceDrugName,],
                                mapping = aes(x = time,
                                              y = (median),
                                              ymin = (min),
                                              ymax = (max),
                                              fill = drug),alpha = 0.5) + theme_minimal()

  pkData <- read.csv(pkDataFilePath)
  pkData$drug <- referenceDrugName
  pkData$id <- as.factor(pkData$id)
  plt <- plt + geom_point(data = pkData, mapping = aes(x = time,
                                                       y = outputValues,
                                                       color = id,
                                                       shape = "Individual"))

  plt <- plt + xlab(paste0("Time (",timeUnit,")")) + ylab(paste0("Reference (",yUnit,")"))
  plt <- plt + guides(color = "none")
  plt <- plt + theme(legend.position = "bottom", legend.direction = "horizontal",
                     axis.text = element_text(size = 12),axis.title = element_text(size = 14),
                     legend.text = element_text(size = 12),legend.title = element_text(size = 14),
                     strip.text = element_text(size = 14))

  plt <- plt + scale_fill_manual(name = referenceDrugName, labels = "Simulated" , values = "darkgrey")
  plt <- plt + scale_shape_manual(name = "", labels = "Observed" , values = 16)

  methods::show(plt)

  return(plt)
}


#' @export
getReferenceVBEPlot <- function(refAndTestSimulationsInVirtualPopulation,
                                referenceDrugName = "Reference",
                                testDrugName = "Test",
                                yDimension = NULL,
                                yUnit = NULL,
                                timeUnit = ospUnits$Time$min){

  plasmaDf <- refAndTestSimulationsInVirtualPopulation


  if (!is.null(yDimension) & !is.null(yUnit)){
    plasmaDf$conc <- toUnit(quantityOrDimension = yDimension,
                            values = plasmaDf$conc,
                            targetUnit = yUnit)
  }


  plasmaDf$time <- toUnit(quantityOrDimension = ospDimensions$Time,
                          values = plasmaDf$time,
                          targetUnit = timeUnit)

  vbemmdf <- NULL

  q05 <- function(x){quantile(x = x,0.05)}
  q95 <- function(x){quantile(x = x,0.95)}
  fnVec <- c("median" = median,"min" = min,"max" = max)
  for (fn in names(fnVec)){
    df <- aggregate(x = plasmaDf$conc,
                    by = list(plasmaDf$time,
                              plasmaDf$drug),
                    fnVec[[fn]])
    names(df) <- c("time","drug",fn)
    vbemmdf[["time"]] <- df[["time"]]
    vbemmdf[["drug"]] <- df[["drug"]]
    vbemmdf[[fn]] <- df[[fn]]
  }
  vbemmdf <- as.data.frame(vbemmdf)
  vbemmdf$drug <- as.character(vbemmdf$drug)
  vbemmdf$drug[vbemmdf$drug == "R"] <- referenceDrugName
  vbemmdf$drug[vbemmdf$drug == "T1"] <- testDrugName

  vbemmdf$drug <- factor(vbemmdf$drug, levels = c(referenceDrugName, testDrugName))
  vbeplt <- ggplot() + geom_ribbon(data = vbemmdf,
                                   mapping = aes(x = time,
                                                 y = (median),
                                                 ymin = (min),
                                                 ymax = (max),
                                                 fill = drug),alpha = 0.33)

  vbeplt <- vbeplt + geom_line(data = vbemmdf,
                               mapping = aes(x = time,
                                             y = (median),
                                             color = drug),
                               size = 2) + theme_minimal()

  vbeplt <- vbeplt + theme(legend.position = "bottom", legend.direction = "horizontal",
                           axis.text = element_text(size = 12),axis.title = element_text(size = 14),
                           legend.text = element_text(size = 12),legend.title = element_text(size = 14),
                           strip.text = element_text(size = 14))

  vbeplt <- vbeplt + xlab(paste0("Time (",timeUnit,")")) + ylab(paste0(yDimension," (",yUnit,")"))
  vbeplt <- vbeplt + scale_fill_manual(name = "Drug" , labels = c(referenceDrugName, testDrugName) , values = c("#ff0000","#1e90ff"))
  vbeplt <- vbeplt + scale_color_manual(name = "Drug" , labels = c(referenceDrugName, testDrugName) , values = c("#ff0000","#1e90ff"))

  methods::show(vbeplt)
  return(vbeplt)
}


#' @export
plotVBEProbability <- function(ctsResults){
  ctsResultsSummary <- ctsResults$summary

  trialNames <- c("auc_par","cmax_par","auc_cross","cmax_cross")

  probSuccessDf <- NULL
  for (trial in trialNames){
    df <- as.data.frame(ctsResultsSummary[[trial]])
    if(!is.data.frame(df)){
      next
    }
    df$trial <- trial
    probSuccessDf <- rbind.data.frame(probSuccessDf,df)
  }

  probSuccessDf$pk <- sapply(probSuccessDf$trial,function(x){if(x %in% c("auc_par","auc_cross")){"AUC"}else{"Cmax"}})
  probSuccessDf$design <- sapply(probSuccessDf$trial,function(x){if(x %in% c("auc_par","cmax_par")){"Parallel"}else{"Crossover"}})

  proPlt <- ggplot(data = probSuccessDf,mapping = aes(x = n_subj,y = pBE_rep)) + geom_line() + geom_point() +  facet_grid( pk ~ design ) + scale_x_continuous(name = "Number of subjects") + scale_y_continuous(name = "Probability of BE")
  proPlt <- proPlt + geom_hline(yintercept = 0.9,linetype = "dashed")
  proPlt <- proPlt + annotate("text", x = 10, y = 0.82, label = "0.9")
  proPlt <- proPlt + theme(axis.title = element_text(size = 14),axis.text = element_text(size = 12), strip.text = element_text(size = 16))
  return(proPlt)

}
getMeanMinMaxDataframe <- function(dataframe, valuesColumnName, byColumnNames,
                                   returnGeometricMean = FALSE,
                                   yminFunction = min, ymaxFunction = max) {
  # Define mean function
  meanFn <- if (returnGeometricMean) {
    function(x) exp(mean(log(x), na.rm = TRUE))
  } else {
    function(x) mean(x, na.rm = TRUE)
  }

  # Create list of aggregation functions
  aggFunctions <- list(ymean = meanFn, ymin = yminFunction, ymax = ymaxFunction)

  # Aggregate data
  aggregated <- lapply(aggFunctions, function(fn) {
    aggregate(dataframe[[valuesColumnName]],
              by = dataframe[byColumnNames],
              FUN = fn)
  })

  # Combine results into a single dataframe
  result <- aggregated[[1]][, byColumnNames, drop = FALSE] # Base structure
  for (name in names(aggFunctions)) {
    result[[name]] <- aggregated[[name]]$x
  }

  return(result)
}


#' @export
plotVBEBands <- function(ctsResults){

  loEnd <- 0.05

  yminFunction <- function(x){quantile(x = x, loEnd )}
  ymaxFunction <- function(x){quantile(x = x, 1 - loEnd)}

  dfAll <- ctsResults$trials
  if("trials" %in% names(dfAll)){
    dfAll <- dfAll[["trials"]]
  }

  trialLabels <- c("auc_par_all" = paste0("AUC (Parallel)"),
                   "auc_cross_all" = paste0("AUC (Crossover)"),
                   "cmax_par_all" = paste0("Cmax (Parallel)"),
                   "cmax_cross_all" = paste0("Cmax (Crossover)"))

  trialNames <- names(trialLabels)

  pltList <- list()

  for (trial in trialNames){

    df <- dfAll[[trial]]

    df <- reshape2::melt(data = df,measure.vars = c("de_ratio","ci_lo","ci_up"))

    df <- getMeanMinMaxDataframe(dataframe = df,
                                 valuesColumnName = "value",
                                 byColumnNames = c("n_subj","variable"),
                                 returnGeometricMean = TRUE,
                                 yminFunction = yminFunction,
                                 ymaxFunction = ymaxFunction)

    plt <- ggplot() + geom_ribbon(data = df,mapping = aes(x = n_subj,y= ymean,ymin = ymin, ymax = ymax, fill = variable),alpha = 0.33)

    plt <- plt + geom_line(data = df,mapping = aes(x = n_subj,y= ymean, color = variable),size = 0.5)
    plt <- plt + scale_fill_manual(name = "X", values =  c(ci_up = "#1e90ff",ci_lo = "#1e90ff",de_ratio  = "#ff0000"))
    plt <- plt + scale_color_manual(name = "X", values =  c(ci_up = "#1e90ff",ci_lo = "#1e90ff",de_ratio  = "#ff0000"))
    plt <- plt + geom_hline(yintercept = c(0.8, 1/0.8),
                            color = "#000000", linetype = "dashed", size = 0.5)
    plt <- plt + xlab("Sample size") + ylab("Ratio")
    plt <- plt + theme_minimal()
    plt <- plt + theme(axis.text = element_text(size = 12),legend.position = "none",
                       panel.background = element_rect(fill = "white", colour = "white"),
                       plot.background = element_rect(fill = "white", colour = "white"),
                       axis.title = element_text(size = 14), strip.text = element_text(size = 14))

    plt <- plt + ggtitle(label = trialLabels[[trial]])

    pltList[[trial]] <- plt

  }

  return(pltList)
}

#' @export
plotClinicalTrialSimulationResults <- function(ctsResults,addPlotHeaders = TRUE){
  pkDictionary <- list(auc_par = "AUC (Parallel)",
                       cmax_par = "Cmax (Parallel)",
                       auc_cross = "AUC (Crossover)",
                       cmax_cross = "Cmax (Crossover)",
                       auc_rep = "AUC (Replicated)",
                       cmax_rep = "Cmax (Replicated)")
  dataDf <- NULL
  pltList <- list()
  for (trial in names(pkDictionary)){
    df <- ctsResults$summary[[trial]]
    if(!is.data.frame(df)){
      next
    }
    df$trial <- pkDictionary[[trial]]
    plt <- ggplot(data = df)
    plt <- plt + geom_ribbon(mapping = aes(x = n_subj , ymin = ci_lo , ymax = ci_up ), fill = "#1e90ff", color = "#1e90ff", alpha = 0.33)
    plt <- plt + geom_line( mapping = aes(x = n_subj , y = de_ratio) )
    plt <- plt + geom_point( mapping = aes(x = n_subj , y = de_ratio) , size = 3 )
    plt <- plt + geom_hline(yintercept = c(0.8,1/0.8) , color = "#ff0000" , linetype = "dashed" , size = 1 )
    plt <- plt + xlab("Sample size") + ylab("Ratio")
    plt <- plt + theme_minimal()
    plt <- plt +  theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      plot.background = element_rect(fill = "white", colour = "white"),
      axis.text = element_text(size = 12),axis.title = element_text(size = 14),
      strip.text = element_text(size = 14))


    if (addPlotHeaders){
      plt <- plt + facet_grid( . ~ trial)
    }
    pltList[[trial]] <- plt
  }
  return(pltList)
}
