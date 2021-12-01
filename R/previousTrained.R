setGeneric("previousTrained", function(classifyResult, ...)
standardGeneric("previousTrained"))

setMethod("previousTrained", "ClassifyResult", 
          function(classifyResult, .iteration, verbose = 3)
{
  if(verbose == 3)
    message("Using existing classifier for classification.")
  
  models(classifyResult)[[.iteration]]
})