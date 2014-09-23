setGeneric("getLocationsAndScales", function(expression, ...)
           {standardGeneric("getLocationsAndScales")})

setMethod("getLocationsAndScales", "matrix", function(expression, ...)
{
  colnames(expression) <- NULL
  rownames(expression) <- NULL
  getLocationsAndScales(ExpressionSet(expression))
})

setMethod("getLocationsAndScales", "ExpressionSet",
          function(expression, location = c("mean", "median"),
                               scale = c("SD", "MAD", "Qn"))
{
  if(scale == "Qn" && !requireNamespace("robustbase", quietly = TRUE))
    stop("The package 'robustbase' could not be found. Please install it.")
  location <- match.arg(location) 
  scale <- match.arg(scale)
  expression <- exprs(expression)
 
  list(switch(location,
              mean = rowMeans(expression),
              median = apply(expression, 1, median)),
       switch(scale,
              SD = apply(expression, 1, sd),
              MAD = apply(expression, 1, mad),
              Qn = apply(expression, 1, robustbase::Qn)))
})