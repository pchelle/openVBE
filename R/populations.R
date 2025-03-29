unifyUnits <- function(values,units,dimension){

  if( is.null(values)){
    return(list(values = NULL, units = NULL))
  }

  allUnits <- unique(units)
  if( length(allUnits) == 1 ){
    return(list(values = values, units = allUnits))
  }

  values <- sapply(seq_along(values),function(n){
    ospsuite::toBaseUnit(quantityOrDimension = dimension,
                         values = values[n],
                         unit = units[n])})
  units <- getBaseUnit(dimension)
  return(list(values = values, units = units ))
}

unifyPopulationDataUnits <- function(demographicData){
  cofactorUnits <- NULL

  weightsDimension <- cofactorDimensions[[standardColumnNames$weightColumn]]
  weights <- unifyUnits(values = demographicData[[standardColumnNames$weightColumn]],
                        units = demographicData[[standardColumnNames$weightUnitsColumn]],
                        dimension = weightsDimension)
  demographicData[[standardColumnNames$weightColumn]] <- weights$values
  demographicData[[standardColumnNames$weightUnitsColumn]] <- weights$units
  cofactorUnits[[standardColumnNames$weightColumn]] <- weights$units

  heightsDimension <- cofactorDimensions[[standardColumnNames$heightColumn]]
  heights <- unifyUnits(values = demographicData[[standardColumnNames$heightColumn]],
                        units = demographicData[[standardColumnNames$heightUnitsColumn]],
                        dimension = heightsDimension)
  demographicData[[standardColumnNames$heightColumn]] <- heights$values
  demographicData[[standardColumnNames$heightUnitsColumn]] <- heights$units
  cofactorUnits[[standardColumnNames$heightColumn]] <- heights$units

  agesDimension <- cofactorDimensions[[standardColumnNames$ageColumn]]
  ages <- unifyUnits(values = demographicData[[standardColumnNames$ageColumn]],
                     units = demographicData[[standardColumnNames$ageUnitsColumn]],
                     dimension = agesDimension)
  demographicData[[standardColumnNames$ageColumn]] <- ages$values
  demographicData[[standardColumnNames$ageUnitsColumn]] <- ages$units
  cofactorUnits[[standardColumnNames$ageColumn]] <- ages$units

  gestationalAgesDimension <- cofactorDimensions[[standardColumnNames$gestationalAgeColumn]]
  gestationalAges <- unifyUnits(values = demographicData[[standardColumnNames$gestationalAgeColumn]],
                                units = demographicData[[standardColumnNames$gestationalAgeUnitsColumn]],
                                dimension = gestationalAgesDimension)
  demographicData[[standardColumnNames$gestationalAgeColumn]] <- gestationalAges$values
  demographicData[[standardColumnNames$gestationalAgeUnitsColumn]] <- gestationalAges$units
  cofactorUnits[[standardColumnNames$gestationalAgeColumn]] <- gestationalAges$units

  return(list(demographicData = demographicData,cofactorUnits = cofactorUnits))
}

getClusteringNPOD <- function(inferredDistribution,parameterNames,cofactorNames,numberOfClusters){

  print("Building cluster model based on NPOD results.")

  weighted_points <- getWeightedPoints(res = inferredDistribution,
                                       params = cofactorNames)

  sampledPoints <- weighted_points$points[sample(1:length(weighted_points$weights),
                                                 size = 100,
                                                 replace = TRUE,
                                                 prob = weighted_points$weights),]
  mclustBIC <- mclust::mclustBIC
  clusters <- mclust::Mclust(data = sampledPoints,G = numberOfClusters)
  plot(clusters, what = "classification")


  return(clusters)
}

getClustersFunctions <- list(NPOD = getClusteringNPOD)

getParameterListNPOD <- function(inferredDistribution){
  return(inferredDistribution$parameters)
}

getParameterListFunctions <- list(NPOD = getParameterListNPOD)


#' @export
getClusters <- function(inferredDistribution,numberOfClusters){
  parameterList <- getParameterListFunctions[[inferredDistribution$method]](inferredDistribution)
  parameterNames <- sapply(parameterList,function(par){par$displayName})
  cofactorNames <- intersect(allCofactorNames,names(inferredDistribution$demographicData))
  clusters <- getClustersFunctions[[inferredDistribution$method]](inferredDistribution,parameterNames,cofactorNames,numberOfClusters)
  return(clusters)
}


#' @export
createVirtualPopulation <- function(inferredDistribution,
                                    proportionOfFemales = 50,
                                    numberOfVirtualIndividuals = 100,
                                    clusters,
                                    demographyRanges){
  print("Creating virtual population characteristics.")
  # create virtual population

  if( "cofactorPaths" %in% names(inferredDistribution) ){

    populationDataframe <- data.frame(IndividualId = seq(0,(numberOfVirtualIndividuals-1))) # data.frame(matrix(ncol = 0, nrow = numberOfVirtualIndividuals))

    for (cofactor in names(inferredDistribution$cofactorPaths)){
      geomean <-  mean(c( log( min(demographyRanges[[cofactor]]$range) ) , log(max(demographyRanges[[cofactor]]$range)) ))
      geosd <- (log(max(demographyRanges[[cofactor]]$range)) - log(min(demographyRanges[[cofactor]]$range)))/(2*1.96)
      values <- rlnorm(n = numberOfVirtualIndividuals,meanlog = geomean,sdlog = geosd)
      populationDataframe[[ cofactorPaths[[cofactor]] ]] <- toBaseUnit(quantityOrDimension = cofactorDimensions[[cofactor]],
                                                                       values = values,
                                                                       unit = demographyRanges[[cofactor]]$units)
    }

    cofactorPathDictionary <- cofactorPaths

  } else {

    weightRange = demographyRanges$weight$range
    weightUnits = demographyRanges$weight$units
    heightRange = demographyRanges$height$range
    heightUnits = demographyRanges$height$units
    ageRange = demographyRanges$age$range
    ageUnits = demographyRanges$age$units
    gestationalAgeRange = demographyRanges$gestationalAge$range
    gestationalAgeUnits = demographyRanges$gestationalAge$units
    BMIRange = demographyRanges$BMI$range
    BMIUnits = demographyRanges$BMI$units

    virtualPopnChars <- createPopulationCharacteristics(species = unique(inferredDistribution$demographicData$species),
                                                        population = unique(inferredDistribution$demographicData$population),
                                                        numberOfIndividuals = numberOfVirtualIndividuals,
                                                        proportionOfFemales = proportionOfFemales,
                                                        weightMin = weightRange[1] %||% minNull(inferredDistribution$demographicData[[standardColumnNames$weightColumn]]),
                                                        weightMax = weightRange[2] %||% maxNull(inferredDistribution$demographicData[[standardColumnNames$weightColumn]]),
                                                        weightUnit = weightUnits %||% inferredDistribution$cofactorUnits[[standardColumnNames$weightColumn]] %||%  ospUnits$Mass$kg,
                                                        heightMin = heightRange[1] %||% minNull(inferredDistribution$demographicData[[standardColumnNames$heightColumn]]),
                                                        heightMax = heightRange[2] %||% maxNull(inferredDistribution$demographicData[[standardColumnNames$heightColumn]]),
                                                        heightUnit = heightUnits %||% inferredDistribution$cofactorUnits[[standardColumnNames$heightColumn]] %||%  ospUnits$Length$cm,
                                                        ageMin = ageRange[1] %||% minNull(inferredDistribution$demographicData[[standardColumnNames$ageColumn]]),
                                                        ageMax = ageRange[2] %||% maxNull(inferredDistribution$demographicData[[standardColumnNames$ageColumn]]),
                                                        ageUnit = ageUnits %||% inferredDistribution$cofactorUnits[[standardColumnNames$ageColumn]] %||%  ospUnits$`Age in years`$`year(s)`,
                                                        gestationalAgeMin = gestationalAgeRange[1] %||% minNull(inferredDistribution$demographicData[[standardColumnNames$gestationalAgeColumn]]),
                                                        gestationalAgeMax = gestationalAgeRange[2] %||% maxNull(inferredDistribution$demographicData[[standardColumnNames$gestationalAgeColumn]]),
                                                        gestationalAgeUnit =  gestationalAgeUnits %||% inferredDistribution$cofactorUnits[[standardColumnNames$gestationalAgeColumn]] %||%  ospUnits$`Age in weeks`$`week(s)`,
                                                        BMIMin = BMIRange[1],
                                                        BMIMax = BMIRange[2],
                                                        BMIUnit = BMIUnits %||% ospUnits$BMI$`kg/m²` )

    print("Creating virtual population.")
    virtualPopn <- createPopulation(populationCharacteristics = virtualPopnChars)
    populationDataframe <- populationToDataFrame(population = virtualPopn$population)
    cofactorPathDictionary <- individualParameterPaths
  }
  referencePopulationDataframe <- populationDataframe
  testPopulationDataframe <- populationDataframe

  # update the parameter paths to match new simulation
  parameterList <- getParameterListFunctions[[inferredDistribution$method]](inferredDistribution)
  referenceSimulationParameterPaths <- getParameterPathsInReferenceSimulation(parameterList)
  testSimulationParameterPaths <- getParameterPathsInTestSimulation(parameterList)

  parameterPaths <- names(testPopulationDataframe)
  for (j in seq_along(parameterList)){
    parameterPaths[parameterPaths == referenceSimulationParameterPaths[j]] <- testSimulationParameterPaths[j]
  }
  colnames(testPopulationDataframe) <- parameterPaths

  cofactorNames <- intersect(allCofactorNames,names(inferredDistribution$demographicData))

  # find theta values of closest individual in reference population
  for (i in 1:numberOfVirtualIndividuals) {
    # read in ith person's weight and height

    cofactorValues <- sapply(cofactorNames,function(cof){
      dimension <- cofactorDimensions[[cof]]
      toUnit(quantityOrDimension = dimension,
             values = populationDataframe[[cofactorPathDictionary[[cof]]]][i],
             sourceUnit = getBaseUnit(dimension),
             targetUnit = inferredDistribution$cofactorUnits[[cof]])})

    givenPoints <- c(rep(NA,length(parameterList)),cofactorValues)

    #Sample from conditional distribution and ensure that sample is within parameter bounds
    print(paste0("Sampling from conditional distribution for individual ",i," of ",numberOfVirtualIndividuals,"."))
    thetaSamples <- samplesFromMixture(mclstResults = clusters,
                                       givenPoints = givenPoints,
                                       numberOfSamples = 1,
                                       lowerBounds = sapply(parameterList, function(x){x$lowerBound}),
                                       upperBounds = sapply(parameterList, function(x){x$upperBound}))


    for (j in seq_along(parameterList)){
      referencePopulationDataframe[[ referenceSimulationParameterPaths[j] ]][i] <- thetaSamples[,j]
      testPopulationDataframe[[ testSimulationParameterPaths[j] ]][i] <- thetaSamples[,j]
    }

  }

  return(list(referencePopulationDataframe = referencePopulationDataframe,testPopulationDataframe = testPopulationDataframe))
}


#' @export
simulateVirtualPopulation <- function(referenceSimulationFilePath,
                                      testSimulationFilePath,
                                      referencePopulationDataframe,
                                      testPopulationDataframe,
                                      outputPath,
                                      startTime,
                                      endTime,
                                      resolutionPtsMin){

  # load updated population
  write.csv(referencePopulationDataframe, file = "referenceVirtualPopulation.csv", row.names = FALSE)
  write.csv(testPopulationDataframe, file = "testVirtualPopulation.csv", row.names = FALSE)
  referenceVirtualPopulation <- loadPopulation("referenceVirtualPopulation.csv")
  testVirtualPopulation <- loadPopulation("testVirtualPopulation.csv")
  file.remove("referenceVirtualPopulation.csv")
  file.remove("testVirtualPopulation.csv")

  print("Simulating virtual population.")
  referenceSimulation <- loadSimulation(referenceSimulationFilePath)
  referenceSimulation$outputSelections$clear()
  addOutputs(quantitiesOrPaths = outputPath,simulation = referenceSimulation)

  testSimulation <- loadSimulation(testSimulationFilePath)
  testSimulation$outputSelections$clear()
  addOutputs(quantitiesOrPaths = outputPath,simulation = testSimulation)


  # set output interval
  setOutputInterval(simulation = referenceSimulation, startTime = startTime, endTime = endTime, resolution = resolutionPtsMin)
  setOutputInterval(simulation = testSimulation, startTime = startTime, endTime = endTime, resolution = resolutionPtsMin)

  # generate plasma concentration time profiles
  referenceSimulationResults <- runSimulation(simulation = referenceSimulation, population = referenceVirtualPopulation)
  testSimulationResults <- runSimulation(simulation = testSimulation, population = testVirtualPopulation)

  referenceResultsData <- getOutputValues(referenceSimulationResults, quantitiesOrPaths = outputPath)$data
  referenceResultsData$drug <- "R"
  referenceResultsData$period <- 1

  testResultsData <- getOutputValues(testSimulationResults, quantitiesOrPaths = outputPath)$data
  testResultsData$drug <- "T1"
  testResultsData$period <- 1

  combinedResultsData <- rbind(referenceResultsData,testResultsData)

  resultsData <- data.frame(id = 1 + combinedResultsData[["IndividualId"]] ,
                            time = combinedResultsData[["Time"]],
                            conc = combinedResultsData[[outputPath]],
                            per = combinedResultsData[["period"]],
                            drug = combinedResultsData[["drug"]])

  resultsData$id <- as.factor(resultsData$id)
  resultsData$per <- as.factor(resultsData$per)
  resultsData$drug <- as.factor(resultsData$drug)

  return(resultsData)

}

#' @export
runClinicalTrialSimulation <- function(virtualPopulationSimulationResults,
                                       n_trials = 50,
                                       subj_min = 5,
                                       subj_max = 50,
                                       subj_step = 5){
  print("Running clinical trial simulation.")
  res <- calc_be(virtualPopulationSimulationResults, n_trials = n_trials,subj_min = subj_min,subj_max = subj_max,subj_step = subj_step)
  extRes <- extract_be(res)
  # make_report(sim=virtualPopulationSimulationResults,
  #             be=extRes)
  return(extRes)
}



