setGeneric("limmaSelection", function(measurements, ...)
           {standardGeneric("limmaSelection")})

# Matrix of numeric measurements.
setMethod("limmaSelection", "matrix", function(measurements, classes, ...)
{
  if(is.null(names(classes))) names(classes) <- colnames(measurements)
  limmaSelection(MultiAssayExperiment(list(matrix = measurements), DataFrame(class = classes)),
                 targets = "matrix", ...)
})

# One or more omics datasets, possibly with clinical data.
setMethod("limmaSelection", "MultiAssayExperiment", 
          function(measurements, datasetName, targets = NULL,
                   trainParams, predictParams, resubstituteParams, ...,
                   selectionName = "Moderated t-test", verbose = 3)
{
  if(!requireNamespace("limma", quietly = TRUE))
    stop("The package 'limma' could not be found. Please install it.")
  if(is.null(targets))
    stop("'targets' must be specified but was not.")
  if(length(setdiff(targets, names(measurements))))
    stop("Some values of 'targets' are not names of 'measurements' but all must be.")                            
            
  if(verbose == 3)
    message("Doing feature selection.")
            
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  measurementsTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]

  fitParams <- list(t(as.matrix(measurementsTable)), model.matrix(~ classes))
  if(!missing(...))
    fitParams <- append(fitParams, ...)
  linearModel <- do.call(limma::lmFit, fitParams)
  linearModel <- limma::eBayes(linearModel)
  orderedFeatures <- match(rownames(limma::topTable(linearModel, 2, number = Inf, sort.by = "p")),
                           colnames(measurementsTable))

  .pickFeatures(measurementsTable, classes,
                datasetName, trainParams, predictParams, resubstituteParams,
                orderedFeatures, selectionName, verbose)
})