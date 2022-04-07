

#' Union of BiocParallelParam and missing
#' 
#' Allows a method parameter to be either a BiocParallelParam class object or
#' missing.
#' 
#' 
#' @name BiocParallelParamOrMissing
#' @aliases BiocParallelParamOrMissing BiocParallelParamOrMissing-class
#' @docType class
#' @author Dario Strbenac
NULL






#' Union of character and missing
#' 
#' Allows a method parameter to be either a character or missing.
#' 
#' 
#' @name characterOrMissing
#' @aliases characterOrMissing characterOrMissing-class
#' @docType class
#' @author Dario Strbenac
#' @examples
#' 
#'   ModellingParams("upsample")
#'   ModellingParams() # Default is "downsample" if it is missing.
#' 
NULL



#' Union of A DataFrame Object and NULL
#' 
#' Allows a slot to be either a DataFrame class object or empty. No
#' constructor.
#' 
#' 
#' @name DataFrameOrNULL
#' @aliases DataFrameOrNULL DataFrameOrNULL-class
#' @docType class
#' @author Dario Strbenac
NULL










#' Union of numeric and missing
#' 
#' Allows a method parameter to be either a numeric value or missing.
#' 
#' 
#' @name numericOrMissing
#' @aliases numericOrMissing numericOrMissing-class
#' @docType class
#' @author Dario Strbenac
#' @examples
#' 
#'   CrossValParams(permutations = 10) # numeric
#'   CrossValParams() # missing
#' 
NULL


