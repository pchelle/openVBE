#' @export
SimOutput <- R6::R6Class(
  "SimOutput",
  public = list(

    initialize = function(name,
                          path,
                          dimension,
                          displayUnit) {
      private$.name <- name
      private$.path <- path
      private$.dimension <- dimension
      private$.displayUnit <- displayUnit
    }


  ),

  active = list(

    name = function(value){
      if(missing(value)){
        return(private$.name)
      }
    },

    path = function(value){
      if(missing(value)){
        return(private$.path)
      }
    },

    dimension = function(value){
      if(missing(value)){
        return(private$.dimension)
      }
    },

    displayUnit = function(value){
      if(missing(value)){
        return(private$.displayUnit)
      }
    }


  ),

  private = list(
    .name = NULL,
    .path = NULL,
    .dimension = NULL,
    .displayUnit = NULL
  )
)
