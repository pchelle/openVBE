#' @export
NPODObject <- R6::R6Class(
  classname = "NPODObject",
  public = list(
    simulation = NULL,
    optimizationParameterList = NULL,
    outputPath = NULL,
    pkData = NULL,
    demographicData = NULL,
    populationData = NULL,
    npodRunSettings = NULL,
    cached_mu = NULL,
    n_err = 0,
    theta_0 = NULL,
    pyl = NULL,
    err_log = rep(list(), 5),
    initialize = function(simulationFilePath,
                          optimizationParameterList,
                          outputPath,
                          pkDataFilePath,
                          studyPopulationDataFilePath = NULL,
                          cofactorPaths = NULL,
                          initialGridSize = NULL,
                          npodRunSettings = NULL) {
      # if there is a populationData check that the individuals in it are the same as in the pkData
      self$npodRunSettings <- npodRunSettings
      if (is.null(npodRunSettings)) {
        self$npodRunSettings <- NPODRunSettings$new()
      }

      self$outputPath <- outputPath

      self$simulation <- ospsuite::loadSimulation(simulationFilePath)
      # self$simulation$outputSelections$clear()
      addOutputs(
        quantitiesOrPaths = self$outputPath,
        simulation = self$simulation
      )


      checkOptimizationParameterList <- all(sapply(optimizationParameterList, inherits, "VBEParameter"))
      if(!checkOptimizationParameterList){
        stop("Not all optimization parameters are of class VBEParameter.")
      }
      self$optimizationParameterList <- optimizationParameterList # check that this is a list of VBEParameter class:

      self$parsePKData(pkDataFilePath)

      self$simulation$outputSchema$addTimePoints(sort(as.numeric(unique(unlist( self$pkData$timeVectorList ))))) # add times

      self$parsePopulationData(studyPopulationDataFilePath, cofactorPaths)

      self$cached_mu <- self$multi_mu
      if (!is.null(self$npodRunSettings$npod_cache)) {
        self$cached_mu <- memoise::memoise(self$multi_mu, cache = self$npodRunSettings$npod_cache)
      }

      self$setGrid(initialGridSize)
    },
    parsePKData = function(pkDataFilePath) {

      pkDataDf <- read.csv(pkDataFilePath)

      pkDataDf$time <- sapply(1:nrow(pkDataDf), function(nr) {
        toBaseUnit(
          quantityOrDimension = "Time",
          values = pkDataDf$time[nr],
          unit = pkDataDf$timeUnits[nr]
        )
      })
      pkDataDf$timeUnits <- getBaseUnit("Time")

      outputQuantity <- getQuantity(path = self$outputPath, container = self$simulation)

      pkDataDf$outputValues <- sapply(1:nrow(pkDataDf), function(nr) {
        value <- toBaseUnit(
          quantityOrDimension = outputQuantity$dimension,
          values = pkDataDf$outputValues[nr],
          unit = pkDataDf$outputUnits[nr]
        )
        return(value)
      })

      pkDataDf$outputUnits <- sapply(1:nrow(pkDataDf), function(nr) {
        unit <- getBaseUnit(outputQuantity$dimension)
        return(unit)
      })

      uniqueIndividualIds <- unique(pkDataDf$id)

      numberOfIndividuals <- length(uniqueIndividualIds)

      timeVectorList <- vector(mode = "list", length = numberOfIndividuals)
      plasmaConcentrationVectorList <- vector(mode = "list", length = numberOfIndividuals)
      idList <- vector(mode = "list", length = numberOfIndividuals)

      min_y <- min(na.omit(pkDataDf$outputValues))

      for (idNo in seq_along(uniqueIndividualIds)) {
        individualId <- uniqueIndividualIds[idNo]
        individualTimeVector <- pkDataDf$time[pkDataDf$id == individualId & pkDataDf$outputValues != 999 & !is.na(pkDataDf$outputValues)]
        individualPlasmaConcentrationProfile <- pkDataDf$outputValues[pkDataDf$id == individualId & pkDataDf$outputValues != 999 & !is.na(pkDataDf$outputValues)]

        if (length(individualPlasmaConcentrationProfile) == 0) {
          next
        }

        timeVectorList[[idNo]] <- individualTimeVector
        plasmaConcentrationVectorList[[idNo]] <- individualPlasmaConcentrationProfile
        idList[[idNo]] <- individualId
      }

      self$pkData <- list(
        timeVectorList = timeVectorList,
        plasmaConcentrationVectorList = plasmaConcentrationVectorList,
        idList = idList
      )
    },
    parsePopulationData = function(studyPopulationDataFilePath, cofactorPaths = NULL) {
      if (is.null(studyPopulationDataFilePath)) {
        self$populationData <- NULL # data.frame(IndividualId = sapply(self$pkData$idList,function(x){x}))
        return()
      }

      self$demographicData <- read.csv(studyPopulationDataFilePath)

      if (!is.null(cofactorPaths)) {
        populationDataFrame <- data.frame(IndividualId = self$demographicData$id)
        for (cofactor in names(cofactorPaths)) {
          populationDataFrame[[cofactorPaths[[cofactor]]]] <- self$demographicData[[cofactor]]
        }
        self$populationData <- populationDataFrame
        return()
      }

      individualCharacteristics <- getIndividualCharacteristics(studyPopulationData = self$demographicData)

      individualList <- list()
      populationDataFrame <- NULL
      for (indNo in seq_along(individualCharacteristics)) {
        id <- individualCharacteristics[[indNo]]$id
        individual <- createIndividual(individualCharacteristics[[indNo]]$characteristics)
        individualDf <- cbind(data.frame(t(individual$derivedParameters$values)), data.frame(t(individual$distributedParameters$values)))
        names(individualDf) <- c(individual$derivedParameters$paths, individual$distributedParameters$paths)
        individualDf <- cbind(data.frame(IndividualId = id), individualDf)
        populationDataFrame <- rbind.data.frame(populationDataFrame, individualDf)
      }

      self$populationData <- populationDataFrame
    },
    setGrid = function(initialGridSize) {

      numberOfParameters <- length(self$optimizationParameterList)

      size_theta0 <- initialGridSize
      if (is.null(size_theta0)) {
        size_theta0 <- length(self$pkData$idList)
      }

      if (numberOfParameters == 1) {
        uniformMatrix <- t(DiceDesign::runif.faure(size_theta0, 2)$design)
        uniformMatrix <- matrix(uniformMatrix[1, ], ncol = length(uniformMatrix[1, ]))
      } else {
        uniformMatrix <- t(DiceDesign::runif.faure(size_theta0, numberOfParameters)$design)
      }

      parameterMatrix <- uniformMatrix

      for (parNo in seq_along(self$optimizationParameterList)){

        a <- self$optimizationParameterList[[parNo]]$lowerBound
        b <- self$optimizationParameterList[[parNo]]$upperBound
        logInc <- self$optimizationParameterList[[parNo]]$logIncrement

        if(logInc){
          log_a <- log(a)
          log_b <- log(b)
          parameterMatrix[parNo,] <- exp(log_a + uniformMatrix[parNo,]*(log_b - log_a))
        } else{
          parameterMatrix[parNo,] <- a + uniformMatrix[parNo,]*(b-a)
        }

      }

      self$theta_0 <- parameterMatrix

    },
    runNPOD = function() {
      npodResults <- self$Dopt()
      return(npodResults)
    },
    Dopt = function() {
      print("Entering Dopt")

      old_theta <- self$theta_0
      counter <- 1
      F0 <- -10^(30)
      F1 <- 2 * F0

      options <- neldermead::optimset(MaxFunEvals = 2000000000, TolX = 1e-14, MaxIter = self$npodRunSettings$MaxIter, TolFun = 1e-14)
      objective_function_values <- c() #Initialize a vector to store the objective function values at each iteration
      objective_function_values[counter] <- F0
      objective_function_values[counter + 1] <- F1
      while (abs(objective_function_values[counter + 1] - objective_function_values[counter]) > self$npodRunSettings$theta_F) {

        print("Starting optimization run:")
        print("Old theta:")
        print(old_theta)
        #Calculate the matrix of log likelihoods P1 for the set of support points old_theta
        P1 <- self$PSI_2(theta = old_theta)

        burke_results_old_theta <- burke(P1)
        lam1 <- burke_results_old_theta$lambda  #Dual variable obtained from burke optimization with theta = old_theta

        #CONDENSE - remove low probability support points
        L1 <- 0.00000001
        L2 <- max(lam1)/1000
        ind1 <- (lam1 > L1) & (lam1 > L2)  #Find elements of lambda that are greater than the lower bounds L1 and L2
        inb_theta <- old_theta[, ind1, drop = FALSE]  #Select the columns of theta corresponding to values of lambda that fall within the constraints L_lower and L_upper
        print("After Condense 2")
        print(inb_theta)

        #Calculate the matrix of log likelihoods P2 for the set of support points inb_theta
        P2 <- self$PSI_2(theta = inb_theta)

        burke_results_inb_theta <- burke(P2)
        lam2 <- burke_results_inb_theta$lambda  #Dual variable obtain from burke optimization with theta = inb_theta
        objective_function_values[counter + 2] <- burke_results_inb_theta$fobj # Add new objective value to objective_function_values

        #CONDENSE - remove low probability support points
        L3 <- max(lam2)/1000
        ind2 <- lam2 > L3 #Find elements of lambda that are greater than the lower bound L3
        new_weights <- lam2[ind2] / sum(lam2[ind2])
        new_theta <- inb_theta[, ind2, drop = FALSE]
        print("After Condense 2")
        print(new_theta)

        if (abs(objective_function_values[counter + 2] - objective_function_values[counter + 1]) <= self$npodRunSettings$theta_F) {
          print("Stopping: Absolute change in objective function is smaller than theta_F.")
          print(new_theta)
          break
        }

        if (counter >= self$npodRunSettings$ncycles) {
          print("Stopping: Maximum number of cycles reached.")
          print(new_theta)
          break
        }

        K <- length(new_theta[1, ])

        self$pyl <- P2[, ind2, drop = FALSE] %*% new_weights # P2 %*% new_weights

        print("Starting fminsearch and prune")
        for (l in 1:K) {
          print("Starting fminsearch...")

          print(l)
          print(new_theta)

          cand_theta <- tryCatch(
            {
              neldermead::fminsearch(self$multi_D, new_theta[, l], options)
            },
            error = function(e) {
              NULL
            }
          )

          print("fminserarch complete")
          print("prune:")

          print("prune...")
          new_theta <- self$prune(theta = new_theta,theta_plus = cand_theta$optbase$xopt)
          print("...complete.")
          print(new_theta)
        }
        print("Completed fminsearch and prune")

        old_theta <- new_theta

        counter <- counter + 1
        print("Counter: ")
        print(counter)
      }
      print("Exiting Dopt")
      return(list("count" = counter, "theta" = new_theta, "w" = new_weights, "LogLikelihood" = objective_function_values[length(objective_function_values)], "PSI" = P2, demographicData = self$demographicData))
    },
    prune = function(theta, theta_plus){
      if(is.null(theta_plus)){
        return(theta)
      }
      print("Entering prune")
      # The `prune` function checks if a candidate support point (`theta_plus`)
      # is sufficiently different from the existing set of support points (`theta`) and within
      # specified bounds (`a`, `b`). If these conditions are met, `theta_plus` is
      # added as a new column to `theta`.

      # Ensure theta_plus is treated as a column vector
      theta_plus <- as.numeric(theta_plus)

      # Determine lower and upper bounds with log scaling where applicable
      a <- sapply(self$optimizationParameterList, function(par){ifelse(test = isTRUE(par$logIncrement),log10Safe(par$lowerBound),par$lowerBound)})
      b <- sapply(self$optimizationParameterList, function(par){ifelse(test = isTRUE(par$logIncrement),log10Safe(par$upperBound),par$upperBound)})

      # Scale theta_plus based on logIncrement flag
      theta_plus_scaled <- sapply(seq_along(theta_plus), function(n){ifelse(test = isTRUE(self$optimizationParameterList[[n]]$logIncrement),log10Safe(theta_plus[n]),theta_plus[n])})

      # Scale theta matrix based on logIncrement flag
      theta_scaled <- theta
      for (n in seq_len(nrow(theta_scaled))){
        if(self$optimizationParameterList[[n]]$logIncrement){
          theta_scaled[n,] <- log10Safe(theta_scaled[n,])
        }
      }

      if (ncol(theta_scaled) > 0) {
        # Broadcast theta_plus_scaled to the same dimensions as theta_scaled
        distances <- colSums(abs(theta_scaled - matrix(theta_plus_scaled, nrow = nrow(theta_scaled), ncol = ncol(theta_scaled), byrow = FALSE)) / (b - a))
        min_dist <- min(distances)
      } else {
        min_dist <- Inf
      }

      # Boundary checks
      within_bounds <- all(theta_plus_scaled > a & theta_plus_scaled < b)
      print("Exiting prune")
      # Add theta_plus_scaled if conditions are met
      if (min_dist > self$npodRunSettings$theta_d && within_bounds) {
        print(paste0("Adding new support point: (", paste(theta_plus, collapse=", "), ")"))
        theta <- cbind(theta, theta_plus)
        return(theta)
      }

      #Otherwise reject candidate support point and provide reason.
      print(paste0("Rejecting candidate support point: (", paste(theta_plus, collapse=", "), ")."))
      if(!within_bounds){
        print("Reason: Outside of bounds.")
      } else {
        print("Reason: too close to existing support point.")
      }
      return(theta)

    },
    multi_D = function(theta_parameter) {
      print("Entering multi_D")
      # multi_D computes a likelihood-based metric (D_comp) by summing the
      # probabilities from PSI_2, normalized by self$pyl. It starts
      # with a penalty term (-N), calls PSI_2, and updates D_comp.
      N <- length(self$pkData$plasmaConcentrationVectorList) # nsub
      D_comp <- -N
      D_comp <- D_comp + sum(self$PSI_2(theta = theta_parameter)/self$pyl)
      print("Exiting multi_D")
      return(D_comp)
    },
    PSI_2 = function(theta, ySimList = NULL, useLog = FALSE) {
      print("Entering PSI_2")
      # mprob is an N x K matrix where each entry [i, l] represents the probability
      # of observing y[[i]] given the model prediction ySimList[[i, l]].
      # A small lower bound (1e-100) is enforced for numerical stability.

      y <- self$pkData$plasmaConcentrationVectorList
      t <- self$pkData$timeVectorList

      # Determine K efficiently
      K <- if (is.matrix(theta)) ncol(theta) else length(theta)

      # Compute ySimList if not provided
      if (is.null(ySimList)) {
        t1 <- system.time({
          ySimList <- self$cached_mu(theta)
        })
      }

      #Optimize sigma
      sigma <- list()
      #sigma1 <- list()
      for (ind in seq_len(nrow(ySimList))){
        sig <- 0
        for (sp in seq_len(ncol(ySimList))){
          if(useLog){
            sig  <- sig + (((log(ySimList[[ind,sp]]) - log(y[[ind]]))^2)/ncol(ySimList))
          } else {
            sig  <- sig + (((ySimList[[ind,sp]] - y[[ind]])^2)/ncol(ySimList))
          }
        }
        #sigma1[[ind]] <- sqrt(sig)
        sigma[[ind]] <- sqrt(sig)
        sigma[[ind]] <- 0.01 + (0.4*y[[ind]]) + 0.5*min(unlist(y))
      }

      # Vectorized probability computation
      t1 <- system.time({
        mprob <- sapply(1:K, function(l) {
          sapply(1:length(y), function(i) {
            max(1e-100, self$getLikelihood(y = y[[i]] , sigma = sigma[[i]] , ySim = ySimList[[i, l]] , useLog = useLog))
          })
        })
      })
      print("Exiting PSI_2")
      return(mprob)
    },
    getLikelihood = function(y, t, theta, sigma, ySim, useLog = FALSE) {
      # The getLikelihood function calculates the likelihood of the observed data (y)
      # given the model predictions (ySim) and measurement error (sigma).
      # The default case uses arithmetic error, but if use_log = TRUE,
      # it calculates the likelihood using the logarithmic (ratio) error.
      j <- length(y)
      z <- matrix(0, 1, j)

      # Ensure sizes match
      if (length(y) != length(ySim)) {
        warning("Observed and predicted sizes don't match")
        stop()
      }

      if (useLog) {
        # Log-transformed error (ratio error)
        z <- (log(y) - log(ySim))^2
        logProb <- -log((sqrt(2 * pi) * sigma * y)) - (z / (2 * (sigma^2)))
        # Calculate the final probability by exponentiating the log likelihood and summing over all data points
        prob <- exp(sum(logProb[sigma != 0 & y != 0]))
      } else {
        # Arithmetic error (default case)
        z <- (y - ySim)^2
        logProb <- -log((sqrt(2 * pi) * sigma)) - (z / (2 * (sigma^2)))
        # Calculate the final probability by exponentiating the log likelihood and summing over all data points
        prob <- exp(sum(logProb[sigma != 0]))
      }


      return(prob)
    },
    multi_mu = function(theta) {

      timePointsList <- self$pkData$timeVectorList

      print("Entering Multi Mu")
      n_ind <- length(timePointsList)

      inferenceParameterPaths <- sapply(self$optimizationParameterList, function(x) {
        x$pathInReferenceSimulation
      })

      popDf <- data.frame(matrix(ncol = 0, nrow = n_ind))

      if (!is.null(self$populationData)) {
        populationParameterPaths <- setdiff(colnames(self$populationData), c("IndividualId", inferenceParameterPaths))
        for (path in populationParameterPaths) {
          popDf[[path]] <- self$populationData[[path]]
        }
      }

      for (path in inferenceParameterPaths) {
        popDf[[path]] <- 0
      }

      ySimList <- matrix(rep(list(), n_ind * length(theta[1, ])), nrow = n_ind, ncol = length(theta[1, ]))

      bigPopDf <- NULL
      for (sp in seq_len(ncol(theta))){
        df <- popDf
        for (pathNo in seq_along(inferenceParameterPaths)) {
          df[[  inferenceParameterPaths[pathNo]  ]] <- theta[pathNo,sp]
        }
        bigPopDf <- rbind.data.frame(bigPopDf,df)
      }
      bigPopDf <- cbind(data.frame(IndividualId = seq(nrow(bigPopDf))), bigPopDf)

      write.csv(x = bigPopDf, file = "population.csv", row.names = FALSE)

      population <- loadPopulation("population.csv")
      file.remove("population.csv")

      print("Simulating...")
      res <- runSimulation(simulation = self$simulation, population = population)
      print("...complete.")
      resData <- getOutputValues(res, quantitiesOrPaths = self$outputPath)$data

      for (ind in unique(resData$IndividualId)) {
        df <- resData[resData$IndividualId == ind, ]
        individualNumber <- ((ind - 1) %% n_ind) + 1
        ptNumber <- ((ind - 1) %/% n_ind) + 1
        ySimList[[individualNumber, ptNumber]] <- df[[self$outputPath]][df$Time %in% timePointsList[[individualNumber]]]
      }
      # })
      print("Exiting Multi Mu")
      return(ySimList)
    }
  )
)


#' @export
NPODRunSettings <- R6::R6Class(
  classname = "NPODRunSettings",
  public = list(
    noise = NULL,
    model = NULL,
    ncycles = NULL,
    MaxIter = NULL,
    theta_F = NULL,
    theta_d = NULL,
    size_theta0 = NULL,
    cache_folder_name = NULL,
    npod_cache = NULL,
    initialize = function(noise = FALSE,
                          model = NULL,
                          ncycles = Inf,
                          MaxIter = 10,
                          theta_F = 10e-2,
                          theta_d = 10e-4,
                          size_theta0 = NULL,
                          cache_folder_name = NULL) {
      self$noise <- noise
      self$model <- model
      self$ncycles <- ncycles
      self$MaxIter <- MaxIter
      self$theta_F <- theta_F
      self$theta_d <- theta_d
      self$size_theta0 <- size_theta0

      self$npod_cache <- cachem::cache_mem()
      if (!is.null(cache_folder_name)) {
        self$cache_folder_name <- cache_folder_name
        self$npod_cache <- cachem::cache_disk(cache_folder_name)
      }
    }
  )
)


getIndividualCharacteristics <- function(studyPopulationData) {
  individualCharacteristics <- list()
  for (indNo in seq_along(studyPopulationData[[standardColumnNames$idColumn]])) {
    individualCharacteristics[[indNo]] <- list()
    id <- studyPopulationData[[standardColumnNames$idColumn]][indNo]
    individualCharacteristics[[indNo]]$id <- id
    individualCharacteristics[[indNo]]$characteristics <- ospsuite::createIndividualCharacteristics(
      species = studyPopulationData[[standardColumnNames$speciesColumn]][studyPopulationData$id == id],
      population = studyPopulationData[[standardColumnNames$populationColumn]][studyPopulationData$id == id],
      gender = studyPopulationData[[standardColumnNames$genderColumn]][studyPopulationData$id == id],
      weight = studyPopulationData[[standardColumnNames$weightColumn]][studyPopulationData$id == id],
      weightUnit = studyPopulationData[[standardColumnNames$weightUnitsColumn]][studyPopulationData$id == id] %||% ospUnits$Mass$kg,
      height = studyPopulationData[[standardColumnNames$heightColumn]][studyPopulationData$id == id],
      heightUnit = studyPopulationData[[standardColumnNames$heightUnitsColumn]][studyPopulationData$id == id] %||% ospUnits$Length$cm,
      age = studyPopulationData[[standardColumnNames$ageColumn]][studyPopulationData$id == id],
      ageUnit = studyPopulationData[[standardColumnNames$ageUnitsColumn]][studyPopulationData$id == id] %||% ospUnits$`Age in years`$`year(s)`,
      gestationalAge = studyPopulationData[[standardColumnNames$gestationalAgeColumn]][studyPopulationData$id == id] %||% 40,
      gestationalAgeUnit = studyPopulationData[[standardColumnNames$gestationalAgeUnitsColumn]][studyPopulationData$id == id] %||% ospUnits$`Age in weeks`$`week(s)`
    )
  }
  return(individualCharacteristics)
}
