classifyInterface <- function(..., verbose = 3)
{
  if(!requireNamespace("PoiClaClu", quietly = TRUE))
    stop("The package 'PoiClaClu' could not be found. Please install it.")
  
  classification  <- PoiClaClu::Classify(...)
  
  if(verbose == 3)
    message("Poisson LDA classification completed.")
  classification
}