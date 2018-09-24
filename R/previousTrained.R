setGeneric("previousTrained", function(classifyResult, ...)
{standardGeneric("previousTrained")})

setMethod("previousTrained", "ClassifyResult", 
          function(classifyResult, .iteration, verbose = 3)
{
  if(verbose == 3)
    message("Using existing classifier for classification.")
  
  if(length(.iteration) == 1)
    previousModel <- models(classifyResult)[[.iteration]]
  else # Resample index and fold index.
    previousModel <- models(classifyResult)[[.iteration[[1]]]][[.iteration[[2]]]]

  previousModel
})