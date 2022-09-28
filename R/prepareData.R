#' Convert Different Data Classes into DataFrame and Filter Features
#' 
#' Input data could be of matrix, MultiAssayExperiment, or DataFrame format and this
#' function will prepare a DataFrame of features and a vector of outcomes and help
#' to exclude nuisance features such as dates or unique sample identifiers from
#' subsequent modelling.
#' 
#' @aliases prepareData prepareData,matrix-method prepareData,DataFrame-method
#' prepareData,MultiAssayExperiment-method
#' @param measurements Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing all of the data. For a
#' \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns
#' are features.
#' @param outcome Either a factor vector of classes, a \code{\link{Surv}} object, or
#' a character string, or vector of such strings, containing column name(s) of column(s)
#' containing either classes or time and event information about survival.
#' @param outcomeColumns If \code{measurements} is a \code{MultiAssayExperiment}, the
#' names of the column (class) or columns (survival) in the table extracted by \code{colData(data)}
#' that contain(s) the each individual's outcome to use for prediction.
#' @param useFeatures If \code{measurements} is a \code{MultiAssayExperiment},
#' a two-column table of features to use. The first column must have assay names
#' and the second column must have feature names found for that assay. \code{"clinical"} is
#' also a valid assay name and refers to the clinical data table. \code{"all"} is a special
#' keyword that means all features (passing any other filters) of that assay will be used 
#' for modelling. Otherwise, a character vector of feature names to use suffices.
#' @param maxMissingProp Default: 0.0. A proportion less than 1 which is the maximum
#' tolerated proportion of missingness for a feature to be retained for modelling.
#' @param topNvariance Default: NULL. An integer number of most variable features to subset to.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method.
#' @return A list of length two. The first element is a \code{\link{DataFrame}} of features
#' and the second element is the outcomes to use for modelling.
#' @author Dario Strbenac

#' @usage NULL
#' @export
setGeneric("prepareData", function(measurements, ...)
           standardGeneric("prepareData"))

#' @rdname prepareData
#' @export
setMethod("prepareData", "matrix",
  function(measurements, outcome, ...)
{
  prepareData(S4Vectors::DataFrame(measurements, check.names = FALSE), outcome, ...)
})

#' @rdname prepareData
#' @export
setMethod("prepareData", "DataFrame",
  function(measurements, outcome, useFeatures = "all", maxMissingProp = 0.0, topNvariance = NULL)
{
  if(useFeatures != "all") # Subset to only the desired ones.
    measurements <- measurements[, useFeatures]

  # Won't ever be true if input data was MultiAssayExperiment because wideFormat already produces valid names.  
  if(all.equal(colnames(measurements), make.names(colnames(measurements))) != TRUE)
  {
    warning("Unsafe feature names in input data. Converted into safe names.")
    mcols(measurements)$feature <- colnames(measurements) # Save the originals.
    colnames(measurements) <- make.names(colnames(measurements)) # Ensure column names are safe names.
  }
      
 # DataFrame's outcome variable can be character or factor, so it is a bit involved.
  if(is.character(outcome) && length(outcome) > 3 && length(outcome) != nrow(measurements))
    stop("'outcome' is a character variable but has more than one element. Either provide a\n",
         "       one to three column names or a factor of the same length as the number of samples.")

  ## String specifies the name of a single outcome column, typically a class.
  if(is.character(outcome) && length(outcome) == 1)
  {
    outcomeColumn <- match(outcome, colnames(measurements))
    if(is.na(outcomeColumn))
      stop("Specified column name of outcome is not present in the data table.")
    outcome <- measurements[, outcomeColumn]
    measurements <- measurements[, -outcomeColumn, drop = FALSE]
    # R version 4 and greater no longer automatically casts character columns to factors because stringsAsFactors
    # is FALSE by default, so it is more likely to be character format these days. Handle it.
    if(class(outcome) != "factor") # Assume there will be no ordinary regression prediction tasks ... for now.
      outcome <- factor(outcome)
  }
  
  # survival's Surv constructor has two inputs for the popular right-censored data and
  # three inputs for less-common interval data.
  if(is.character(outcome) && length(outcome) %in% 2:3)
  {
    outcomeColumns <- match(outcome, colnames(measurements))
    if(any(is.na(outcomeColumns)))
      stop("Specified column names of outcome is not present in the data table.")
    outcome <- measurements[, outcomeColumns]
    measurements <- measurements[, -outcomeColumns, drop = FALSE]
  }
  
  if(is(outcome, "factor") && length(outcome) > 3 & length(outcome) < nrow(measurements))
    stop("The length of outcome is not equal to the number of samples.")
  
  ## A vector of characters was input by the user. Ensure that it is a factor.
  if(is.character(outcome) & length(outcome) == nrow(measurements))
    outcome <- factor(outcome)
  
  # Outcome has columns, so it is tabular. It is inferred to represent survival data.
  if(!is.null(ncol(outcome)) && ncol(outcome) %in% 2:3)
  {
    # Assume that event status is in the last column (second for two columns, third for three columns)
    numberEventTypes <- length(unique(outcome[, ncol(outcome)]))
    # Could be one or two kinds of events. All events might be uncensored or censored
    # in a rare but not impossible scenario.
    if(numberEventTypes > 2)
      stop("Number of distinct event types in the last column exceeds 2 but must be 1 or 2.")
      
    if(ncol(outcome) == 2) # Typical, right-censored survival data.
      outcome <- survival::Surv(outcome[, 1], outcome[, 2])
    else # Three columns. Therefore, counting process data.
      outcome <- survival::Surv(outcome[, 1], outcome[, 2], outcome[, 3])
  }
  
  # Remove samples with indeterminate outcome.
  dropSamples <- which(is.na(outcome) | is.null(outcome))
  if(length(dropSamples) > 0)
  {
    measurements <- measurements[-dropSamples, ]
    outcome <- outcome[-dropSamples]  
  }
  
  # Remove features with more missingness than allowed.
  nSamples <- nrow(measurements)
  dropFeatures <- which(apply(measurements, 2, function(featureMeasurements) sum(is.na(featureMeasurements)))
                        / nrow(measurements) > maxMissingProp)
  if(length(dropFeatures) > 0)
    measurements <- measurements[, -dropFeatures]
  
  # Use only the most variable features.
  if(!is.null(topNvariance))
  {
    mostVariance <- order(apply(measurements, 2, var, na.rm = TRUE), decreasing = TRUE)[1:topNvariance]
    measurements <- measurements[, mostVariance]
  }
  
  list(measurements = measurements, outcome = outcome)
})

#' @rdname prepareData
#' @export
setMethod("prepareData", "MultiAssayExperiment",
  function(measurements, outcomeColumns = NULL, useFeatures = "all", ...)
{
  if(is.character(useFeatures)) useFeatures <- data.frame(assay = names(measurements), feature = "all")
  omicsTargets <- setdiff(useFeatures[, "assay"], "clinical")
  if(length(omicsTargets) > 0)
  {
    if(any(anyReplicated(measurements[, , omicsTargets])))
      stop("Data set contains replicates. Please remove or average replicate observations and try again.")
  }
  
  if(!is.null(outcomeColumns) && !all(outcomeColumns %in% colnames(MultiAssayExperiment::colData(measurements))))
    stop("Not all column names specified by 'outcomeColumns' found in clinical table.")  
  if(!all(useFeatures[, "assay"] %in% c(names(measurements), "clinical")))
    stop("Some assay names in first column of 'useFeatures' are not assay names in 'measurements' or \"clinical\".")

  clinicalColumns <- colnames(MultiAssayExperiment::colData(measurements))
  if("clinical" %in% useFeatures[, "assay"])
  {
    clinicalRows <- useFeatures[, "assay"] == "clinical"      
    clinicalColumns <- useFeatures[clinicalRows, "feature"]
    useFeatures <- useFeatures[!clinicalRows, ]
  } else {
    clinicalColumns <- NULL
  }
  
  if(nrow(useFeatures) > 0)
  {
    measurements <- measurements[, , unique(useFeatures[, "assay"])]
  
    # Get all desired measurements tables and clinical columns (other than the columns representing outcome).
    # These form the independent variables to be used for making predictions with.
    # Variable names will have names like RNA_BRAF for traceability.
    dataTable <- MultiAssayExperiment::wideFormat(measurements, colDataCols = union(clinicalColumns, outcomeColumns))
    rownames(dataTable) <- dataTable[, "primary"]
    S4Vectors::mcols(dataTable)[, "sourceName"] <- gsub("colDataCols", "clinical", S4Vectors::mcols(dataTable)[, "sourceName"])
    colnames(S4Vectors::mcols(dataTable))[1] <- "assay"
            
    # Sample information variable names not included in column metadata of wide table but only as row names of it.
    # Create a combined column named "feature" which has feature names of the assays as well as the clinical.
    S4Vectors::mcols(dataTable)[, "feature"] <- as.character(S4Vectors::mcols(dataTable)[, "rowname"])
    missingIndices <- is.na(S4Vectors::mcols(dataTable)[, "feature"])
    S4Vectors::mcols(dataTable)[missingIndices, "feature"] <- colnames(dataTable)[missingIndices]
    
    # Finally, a column annotation recording variable name and which table it originated from for all of the source tables.
    S4Vectors::mcols(dataTable) <- S4Vectors::mcols(dataTable)[, c("assay", "feature")]
    
    # Subset to only the desired features.
    useFeaturesSubset <- useFeatures[useFeatures[, "feature"] != "all", ]
    if(nrow(useFeaturesSubset) > 0)
    {
      uniqueAssays <- unique(useFeatures[, "assay"])
      for(filterAssay in uniqueAssays)
      {
        dropFeatures <- S4Vectors::mcols(dataTable)[, "assay"] == filterAssay &
                        !S4Vectors::mcols(dataTable)[, "feature"] %in% useFeatures[useFeatures[, 1] == filterAssay, 2]
        dataTable <- dataTable[, !dropFeatures]
      }
    }
    dataTable <- dataTable[, -match("primary", colnames(dataTable))]
  } else { # Must have only been clinical data.
    dataTable <- MultiAssayExperiment::colData(measurements)
    S4Vectors::mcols(dataTable) <- DataFrame(assay = "clinical", feature = colnames(dataTable))
  }
    
  # Do other filtering and preparation in DataFrame function.
  prepareData(dataTable, outcomeColumns, useFeatures = "all", ...)
})