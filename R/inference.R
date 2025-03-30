#' @export
runDistributionInference <- function(method = "NPOD",
                                     referenceSimulationFilePath,
                                     outputPath,
                                     pkDataFilePath,
                                     studyPopulationDataFilePath = NULL,
                                     cofactorPaths = NULL,
                                     inferenceParameters = NULL,
                                     numberOfIterations = 10,
                                     useLogNormalLikelihood = NULL,
                                     initialGridSize = NULL,
                                     saveResultsPath = NULL,
                                     cacheFolder = NULL
                                     ){


  return(runDistributionInferenceFunctionsList[[method]](referenceSimulationFilePath = referenceSimulationFilePath,
                                                         outputPath = outputPath,
                                                         pkDataFilePath = pkDataFilePath,
                                                         studyPopulationDataFilePath = studyPopulationDataFilePath,
                                                         cofactorPaths = cofactorPaths,
                                                         inferenceParameters = inferenceParameters,
                                                         numberOfIterations = numberOfIterations,
                                                         useLogNormalLikelihood = useLogNormalLikelihood,
                                                         initialGridSize = initialGridSize,
                                                         saveResultsPath = saveResultsPath,
                                                         cacheFolder = cacheFolder
                                                         ))

}

runDistributionInferenceNPOD <- function(referenceSimulationFilePath,
                                         outputPath,
                                         pkDataFilePath,
                                         studyPopulationDataFilePath,
                                         cofactorPaths,
                                         inferenceParameters,
                                         numberOfIterations,
                                         useLogNormalLikelihood,
                                         initialGridSize,
                                         saveResultsPath,
                                         cacheFolder,
                                         ...){

  print("Starting NPOD")

  individualCharacteristics <- NULL
  cofactorUnits <- NULL
  if(!is.null(studyPopulationDataFilePath)){
    demographicData <- read.csv(studyPopulationDataFilePath)
    unifiedCofactorsUnits <- unifyPopulationDataUnits(demographicData)
    cofactorUnits <- unifiedCofactorsUnits$cofactorUnits
  }

  npodSettings <- NPODRunSettings$new(MaxIter = numberOfIterations,
                                      cache_folder_name = cacheFolder)

  npod <- NPODObject$new(simulationFilePath = referenceSimulationFilePath,
                         optimizationParameterList = inferenceParameters,
                         outputPath = outputPath,
                         pkDataFilePath = pkDataFilePath,
                         studyPopulationDataFilePath = studyPopulationDataFilePath,
                         cofactorPaths = cofactorPaths,
                         initialGridSize = initialGridSize,
                         useLogNormalLikelihood = useLogNormalLikelihood,
                         npodRunSettings = npodSettings)

  npodResults <- npod$runNPOD()
  npodResults$parameters <- inferenceParameters
  npodResults$cofactorUnits <- cofactorUnits
  npodResults$cofactorPaths <- cofactorPaths
  npodResults$method <- "NPOD"

  if(!is.null(saveResultsPath)){
    saveRDS(object = npodResults,file = saveResultsPath)
  }
  return(npodResults)
}


runDistributionInferenceFunctionsList <- list(NPOD = runDistributionInferenceNPOD)


