classifyInterface <- function(..., verbose = 3)
{
  classification  <- Classify(...)
  
  if(verbose == 3)
    message("Poisson LDA classification completed.")
  classification
}