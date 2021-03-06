\name{plotFeatureClasses}
\alias{plotFeatureClasses}
\alias{plotFeatureClasses,matrix-method}
\alias{plotFeatureClasses,DataFrame-method}
\alias{plotFeatureClasses,MultiAssayExperiment-method}
\title{Plot Density, Scatterplot, Parallel Plot or Bar Chart for Features By Class}
\description{Allows the visualisation of measurements in the data set. If \code{targets} is of type
             \code{\link{Pairs}}, then a parallel plot is automatically drawn. If it's a single categorical
             variable, then a bar chart is automatically drawn.}
\usage{
  \S4method{plotFeatureClasses}{matrix}(measurements, classes, targets, ...)
  \S4method{plotFeatureClasses}{DataFrame}(measurements, classes, targets, groupBy = NULL,
                        groupingName = NULL, whichNumericFeaturePlots = c("both", "density", "stripchart"),
                    measurementLimits = NULL, lineWidth = 1, dotBinWidth = 1,
                    xAxisLabel = NULL, yAxisLabels = c("Density", "Classes"),
                   showXtickLabels = TRUE, showYtickLabels = TRUE,
                   xLabelPositions = "auto", yLabelPositions = "auto",
                   fontSizes = c(24, 16, 12, 12, 12),
                   colours = c("#3F48CC", "#880015"), showDatasetName = TRUE, plot = TRUE)
  \S4method{plotFeatureClasses}{MultiAssayExperiment}(measurements, targets, groupBy = NULL, groupingName = NULL,
                   showDatasetName = TRUE, ...)
}
\arguments{
  \item{measurements}{A \code{\link{matrix}}, \code{\link{DataFrame}} or a
                     \code{\link{MultiAssayExperiment}} object containing the data.
                     For a matrix, the rows are for features and the columns are for samples.
                     A column with name \code{"class"} must be present in the \code{DataFrame}
                     stored in the \code{colData} slot.}
  \item{classes}{Either a vector of class labels of class \code{\link{factor}} or if
                 the measurements are of class \code{DataFrame} a character vector of
                 length 1 containing the column name in \code{measurement} is also permitted.
                 Not used if \code{measurements} is a \code{MultiAssayExperiment} object.}
  \item{targets}{If \code{measurements} is a \code{matrix} or \code{DataFrame}, then
                a vector of numeric or character indicies or the feature identifiers corresponding to the
                feature(s) to be plotted. If \code{measurements} is a \code{MultiAssayExperiment}, then
                a \code{DataFrame} of 2 columns must be specified. The first column contains the
                names of the tables and the second contains the names of the variables, thus
                each row unambiguously specifies a variable to be plotted.}
  \item{groupBy}{If \code{measurements} is a \code{DataFrame}, then a character vector of
                length 1, which contains the name of a categorical feature, may be specified.
                If \code{measurements} is a \code{MultiAssayExperiment}, then a
                character vector of length 2, which contains the name of a data table as the first
                element and the name of a categorical feature as the second element, may be specified.
                Additionally, the value \code{"clinical"} may be used to refer to the column 
                annotation stored in the \code{colData} slot of the of the \code{MultiAssayExperiment}
                object. A density plot will have additional lines of different line types for
                each category. A strip chart plot will have a separate strip chart created for
                each category and the charts will be drawn in a single column on the
                graphics device. A parallel plot and bar chart plot will similarly be laid out.}
  \item{groupingName}{A label for the grouping variable to be used in plots.}                
  \item{...}{Unused variables by the three top-level methods passed to the internal method
             which generates the plot(s).}
  \item{whichNumericFeaturePlots}{If the feature is a single feature and has numeric measurements,
                           this option specifies which types of plot(s) to draw. The default value is
                           \code{"both"}, which draws a density plot and also a stip chart below
                           the density plot. Other options are \code{"density"} for drawing only
                           a density plot and \code{"stripchart"} for drawing only a strip chart.}
  \item{measurementLimits}{The minimum and maximum expression values to plot. Default: \code{NULL}.
                           By default, the limits are automatically computed from the data values.}
  \item{lineWidth}{Numeric value that alters the line thickness for density plots. Default: 1.}
  \item{dotBinWidth}{Numeric value that alters the diameter of dots in the strip chart. Default: 1.}                  
  \item{xAxisLabel}{The axis label for the plot's horizontal axis. Default: \code{NULL}.}
  \item{yAxisLabels}{A character vector of length 1 or 2. If the feature's measurements are numeric an
                     \code{whichNumericFeaturePlots} has the value \code{"both"}, the first value is the
                     y-axis label for the density plot and the second value is the y-axis label
                     for the strip chart. Otherwise, if the feature's measurements are numeric and
                     only one plot is drawn, then a character vector of length 1 specifies the
                     y-axis label for that particular plot. Ignored if the feature's measurements
                     are categorical.}
  \item{showXtickLabels}{Logical. Default: \code{TRUE}. If set to \code{FALSE}, the x-axis labels
                         are hidden.}
  \item{showYtickLabels}{Logical. Default: \code{TRUE}. If set to \code{FALSE}, the y-axis labels
                         are hidden.}
  \item{xLabelPositions}{Either "auto" or a vector of values. The positions of labels on the x-axis.
                         If "auto", the placement of labels is automatically calculated.}
  \item{yLabelPositions}{Either "auto" or a vector of values. The positions of labels on the y-axis.
                         If "auto", the placement of labels is automatically calculated.}
  \item{fontSizes}{A vector of length 5. The first number is the size of the title.
                   The second number is the size of the axes titles. The third number is
                   the size of the axes values. The fourth number is the size of the
                   legends' titles. The fifth number is the font size of the legend labels.}
  \item{colours}{The colours to plot data of each class in. The length of this vector must be as
                 long as the distinct number of classes in the data set.}
  \item{showDatasetName}{Logical. Default: \code{TRUE}. If \code{TRUE} and the data is in a
                         \code{MultiAssayExperiment} object, the the name of the table in
                         which the feature is stored in is added to the plot title.}                 
  \item{plot}{Logical. Default: \code{TRUE}. If \code{TRUE}, a plot is produced on the
              current graphics device.}  
}
\value{
  Plots are created on the current graphics device and a list of plot objects is
  invisibly returned. The classes of the plot object are determined based on the
  type of data plotted and the number of plots per feature generated. If the
  plotted variable is discrete or if the variable is numeric and one plot type was specified,
  the list element is an object of class \code{ggplot}. Otherwise, if the variable
  is numeric and both the density and stripchart plot types were made, the list element
  is an object of class \code{TableGrob}.
  
  Settling \code{lineWidth} and \code{dotBinWidth} to the same value doesn't result in
  the density plot and the strip chart having elements of the same size. Some manual 
  experimentation is required to get similarly sized plot elements.
}
\author{Dario Strbenac}

\examples{
  # First 25 samples and first 5 genes are mixtures of two normals. Last 25 samples are
  # one normal.
  genesMatrix <- sapply(1:15, function(geneColumn) c(rnorm(5, 5, 1)))
  genesMatrix <- cbind(genesMatrix, sapply(1:10, function(geneColumn) c(rnorm(5, 15, 1))))
  genesMatrix <- cbind(genesMatrix, sapply(1:25, function(geneColumn) c(rnorm(5, 9, 2))))
  genesMatrix <- rbind(genesMatrix, sapply(1:50, function(geneColumn) rnorm(95, 9, 3)))
  rownames(genesMatrix) <- paste("Gene", 1:100)
  colnames(genesMatrix) <- paste("Sample", 1:50)
  classes <- factor(rep(c("Poor", "Good"), each = 25), levels = c("Good", "Poor"))
  plotFeatureClasses(genesMatrix, classes, targets = "Gene 4",
                     xAxisLabel = bquote(log[2]*'(expression)'), dotBinWidth = 0.5)
                     
                     
  
  infectionResults <- c(rep(c("No", "Yes"), c(20, 5)), rep(c("No", "Yes"), c(5, 20)))
  genders <- factor(rep(c("Male", "Female"), each = 10, length.out = 50))
  clinicalData <- DataFrame(Gender = genders, Sugar = runif(50, 4, 10),
                              Infection = factor(infectionResults, levels = c("No", "Yes")),
                            row.names = colnames(genesMatrix))
  plotFeatureClasses(clinicalData, classes, targets = "Infection")
  plotFeatureClasses(clinicalData, classes, targets = "Infection", groupBy = "Gender")
  
  dataContainer <- MultiAssayExperiment(list(RNA = genesMatrix),
                                        colData = cbind(clinicalData, class = classes))
  targetFeatures <- DataFrame(table = "RNA", feature = "Gene 50")                                     
  plotFeatureClasses(dataContainer, targets = targetFeatures,
                     groupBy = c("clinical", "Gender"),
                     xAxisLabel = bquote(log[2]*'(expression)'))
}