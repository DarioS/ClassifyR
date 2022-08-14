# Automated Usage of Previously Created Classifiers
previousTrained <- function(classifyResult, .iteration, verbose = 3)
{
  if(verbose == 3)
    message("Using existing classifier for classification.")
  
  models(classifyResult)[[.iteration]]
}
attr(previousTrained, "name") <- "previousTrained"