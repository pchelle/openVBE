rm(list=ls())
graphics.off()
library(openVBE)
dateTime <- paste0(format(Sys.Date(), "%Y%m%d"), "_", format(Sys.time(), "%H%M%S"))
exampleName <- "bupropion-oral"
subfolder <- file.path("examples",exampleName)

#Data set D1
pkDataFilePath <- file.path(subfolder, "bupropion_pk_data_reference.csv")

#Data set D2
studyPopulationDataFilePath <- file.path(subfolder,"bupropion_study_population_data.csv")

########>   STEP S1   <########
referenceSimulationFilePath <- file.path(subfolder,"bupropion-oral-reference-formulation-SR-PO-150 mg.pkml")
testSimulationFilePath <- file.path(subfolder,"bupropion-oral-test-formulation-ER-PO-150 mg.pkml")



########>   STEP S2   <########

outputPath <- "Organism|PeripheralVenousBlood|Bupropion human 1st order|Plasma (Peripheral Venous Blood)"
#Identified individual-specific parameters that are uncertain and influential with respect to plasma concentration profiles
ScaleFactorX <- VBEParameter$new(pathInReferenceSimulation = "Applications|PO 150 mg - human|SR PO 150 mg - FDA table|ScaleFactorX",
                                 pathInTestSimulation = "Applications|PO 150 mg - human|ER PO 150 mg - FDA table|ScaleFactorX",
                                 displayName = "ScaleFactorX",
                                 dimension = ospDimensions$Dimensionless,
                                 unit = ospUnits$Dimensionless$Unitless,
                                 lowerBound = 0.5,
                                 upperBound = 1.5,
                                 logIncrement = FALSE)

ScaleFactorY <- VBEParameter$new(pathInReferenceSimulation = "Applications|PO 150 mg - human|SR PO 150 mg - FDA table|ScaleFactorY",
                                 pathInTestSimulation =  "Applications|PO 150 mg - human|ER PO 150 mg - FDA table|ScaleFactorY",
                                 displayName = "ScaleFactorY",
                                 dimension = ospDimensions$Dimensionless,
                                 unit = ospUnits$Dimensionless$Unitless,
                                 lowerBound = 0.5,
                                 upperBound = 1.5,
                                 logIncrement = FALSE)

LiverIntestinalCL <- VBEParameter$new(pathInReferenceSimulation = "Liver and Intestinal CL|Reference concentration",
                                      pathInTestSimulation = "Liver and Intestinal CL|Reference concentration",
                                      displayName = "LiverIntestinalCL",
                                      dimension = ospDimensions$`Concentration (molar)`,
                                      unit = ospUnits$`Concentration [molar]`$`µmol/l`,
                                      lowerBound = 0.5,
                                      upperBound = 1.5,
                                      logIncrement = FALSE)
inferenceParameters <- list(ScaleFactorX,ScaleFactorY,LiverIntestinalCL)
inferredDistribution <- runDistributionInference(method = "NPOD",
                                                 referenceSimulationFilePath = referenceSimulationFilePath,
                                                 outputPath = outputPath,
                                                 pkDataFilePath = pkDataFilePath,
                                                 studyPopulationDataFilePath = studyPopulationDataFilePath,
                                                 cofactorPaths = NULL,
                                                 inferenceParameters = inferenceParameters,
                                                 initialGridSize = NULL,
                                                 numberOfIterations = 10,
                                                 saveResultsPath = file.path(subfolder,paste0(exampleName,"-npod-results-",dateTime,".rds")),
                                                 cacheFolder = file.path("examples/cacheFolder"))


########>   STEP S3   <########
demographyRanges <- list(weight = list(range = c(50, 110),
                                       units = "kg"),
                         height = list(range = c(16,20),
                                       units = "dm"))

clusters <- getClusters(inferredDistribution = inferredDistribution,
                        numberOfClusters = 2)
virtualPopulation <- createVirtualPopulation(inferredDistribution = inferredDistribution,
                                             proportionOfFemales = 50,
                                             numberOfVirtualIndividuals = 100,
                                             clusters = clusters,
                                             demographyRanges = demographyRanges)

refAndTestSimulationsInVirtualPopulation <- simulateVirtualPopulation(referenceSimulationFilePath = referenceSimulationFilePath,
                                                                      testSimulationFilePath = testSimulationFilePath,
                                                                      referencePopulationDataframe = virtualPopulation$referencePopulationDataframe,
                                                                      testPopulationDataframe = virtualPopulation$testPopulationDataframe,
                                                                      outputPath = outputPath,
                                                                      startTime = 0,
                                                                      endTime = 5760,
                                                                      resolutionPtsMin =1/15)  # pts/min


########>   PLOTTING   <########
source(file.path(subfolder,"plotting-script-vbe-bupropion.R"))


########>   STEP S4   <########
ctsResults <- runClinicalTrialSimulation(virtualPopulationSimulationResults = refAndTestSimulationsInVirtualPopulation,n_trials = 500)

pltList <- plotClinicalTrialSimulationResults(ctsResults)
show(pltList)

bandsPlot <- plotVBEBands(ctsResults)
show(bandsPlot)

probabilityPlot <- plotVBEProbability(ctsResults)
show(probabilityPlot)




