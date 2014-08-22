classify <- function(...)
{
  params <- list(...)
  params <- params[names(params) != "verbose"]
  do.call(Classify, params)
}