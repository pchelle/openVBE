#' @export
VBEParameter <- R6::R6Class(classname = "VBEParameter",
                         public = list(
                           initialize = function(pathInReferenceSimulation,
                                                 pathInTestSimulation,
                                                 displayName = NULL,
                                                 dimension,
                                                 unit,
                                                 lowerBound = 1e-1,
                                                 upperBound = 1e1,
                                                 logIncrement = TRUE){
                             self$pathInReferenceSimulation <- pathInReferenceSimulation
                             self$pathInTestSimulation <- pathInTestSimulation
                             self$displayName <- displayName %||% pathInReferenceSimulation
                             self$dimension <- dimension
                             self$unit <- unit
                             self$lowerBound <- lowerBound
                             self$upperBound <- upperBound
                             self$logIncrement <- logIncrement
                           }
                         ),
                         private = list(
                           .pathInReferenceSimulation = NULL,
                           .pathInTestSimulation = NULL,
                           .displayName = NULL,
                           .dimension = NULL,
                           .unit = NULL,
                           .lowerBound = NULL,
                           .upperBound = NULL,
                           .logIncrement = NULL
                         ),
                         active = list(

                           pathInReferenceSimulation = function(value){
                             if(missing(value)){
                               return(private$.pathInReferenceSimulation)
                             }
                             error(condition = !is.character(value),
                                   errorMessage = "Parameter 'pathInReferenceSimulation' must be of type 'character'.")
                             private$.pathInReferenceSimulation <- value
                           },

                           pathInTestSimulation = function(value){
                             if(missing(value)){
                               return(private$.pathInTestSimulation)
                             }
                             error(condition = !is.character(value),
                                   errorMessage = "Parameter 'pathInTestSimulation' must be of type 'character'.")
                             private$.pathInTestSimulation <- value
                           },

                           displayName = function(value){
                             if(missing(value)){
                               return(private$.displayName)
                             }
                             error(condition = !is.character(value),
                                   errorMessage = "Parameter 'displayName' must be of type 'character'.")
                             private$.displayName <- value
                           },

                           dimension = function(value){
                             if(missing(value)){
                               return(private$.dimension)
                             }
                             error(condition = !is.character(value),
                                   errorMessage = "Parameter 'dimension' must be of type 'character'.")
                             private$.dimension <- value
                           },

                           unit = function(value){
                             if(missing(value)){
                               return(private$.unit)
                             }
                             error(condition = !is.character(value),
                                   errorMessage = "Parameter 'unit' must be of type 'character'.")
                             private$.unit <- value
                           },

                           lowerBound = function(value){
                             if(missing(value)){
                               return(private$.lowerBound)
                             }
                             error(condition = !is.numeric(value),
                                   errorMessage = "Parameter 'lowerBound' must be of type 'numeric'.")
                             private$.lowerBound <- value
                           },

                           upperBound = function(value){
                             if(missing(value)){
                               return(private$.upperBound)
                             }
                             error(condition = !is.numeric(value),
                                   errorMessage = "Parameter 'lowerBound' must be of type 'numeric'.")
                             private$.upperBound <- value
                           },

                           logIncrement = function(value){
                             if(missing(value)){
                               return(private$.logIncrement)
                             }
                             error(condition = !is.logical(value),
                                   errorMessage = "Parameter 'logIncrement' must be of type 'logical'.")
                             private$.logIncrement <- value
                           }

                         )
)


getParameterBounds <- function(parameterList){
  parameterBounds <- list()
  parameterBounds[[1]] <- sapply(parameterList, function(par){par$lowerBound})
  parameterBounds[[2]] <- sapply(parameterList, function(par){par$upperBound})
  return(parameterBounds)
}

getParameterPathsInReferenceSimulation <- function(parameterList){
  parameterPaths <- sapply(parameterList, function(par){par$pathInReferenceSimulation})
  return(parameterPaths)
}

getParameterPathsInTestSimulation <- function(parameterList){
  parameterPaths <- sapply(parameterList, function(par){par$pathInTestSimulation})
  return(parameterPaths)
}

