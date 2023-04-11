\name{ClassifyRnews}
\title{ClassifyR News}
\encoding{UTF-8}

\section{Version 3.4.0}{\itemize{
  \item
  Companion website with more in-depth explanation and examples.
  \item
  Default random forest classifier based on ranger now does two-step classification; one for variable importance and one for model fitting, as recommended by ranger's developer.
  \item
  More functions use automatic parameter value selection as their defaults.
  \item
  \code{randomSelection} function to choose random sets of features in cross-validation.
  \item
  \code{crossValidate} function now permits custom parameter tuning via \code{extraParams}.
  \item
  Elastic net GLM and ordinary GLM now calculate class weights be default, so as to perform well in class-imbalanced scenarios.
  \item
  \code{precisionPathwaysTrain} and \code{precisionPathwaysPredict} functions for building tree-like models of multiple assays and their accessory functions \code{calcCostsAndPerformance}, \code{bubblePlot}, \code{flowchart}, \code{strataPlot} for model performance evaluation.
  \item
  \code{crissCrossValidate} function that takes a list of data sets with the same set of features and the same set of outcomes and does all possible pairs of training and prediction to evaluate generalisability. \code{crissCrossPlot} for visual evaluation.
  \item
  
}}

\section{Version 3.2.0}{\itemize{
  \item
  Fast Cox survival analysis.
  \item
  Simple parameter sets, as used by crossVaildate, now come with tuning parameter grid as standard.
  \item
  Wrappers are greatly simplified. Now, there is only one method for a data frame and they are not exported because they are not used directly by the end-user anyway.
  \item
  \code{prepareData} function to filter and subset input data using common ways, such as missingness and variability.
  \item
  Invalid column names of data (e.g. spaces, hyphens) are converted into safe names before modelling but converted back into original names for tracking ranked and selected features.
  \item
  \code{available} function shows the keywords corresponding to transformation, selection, classifier functions.
  \item
  More functions have automatically-selected parameters based on input data, reducing required user-specified parameters.
  \item
  New classifiers added for random survival forests and extreme gradient boosting.
  \item
  Adaptive sampling for modelling with uncertainty of class labels can be enabled with \code{adaptiveResamplingDelta}.
  \item
  Parameter tuning fixed to only use samples from the training set.
}}

\section{Version 3.0.0}{\itemize{
  \item
  Now supports survival models and their evaluation, in addition to existing classification functionality.
  \item
  Cross-validation no longer requires specific annotations like data set name and classifier name. Now, the user can specify any characteristics they want and use these as variables to group by or change line appearances. Also, characteristics like feature selection name and classifier name are automatically filled in from an internal table.
  \item
  Ease of use greatly inproved by \code{crossValidate} function which allows specification of classifiers by a single keyword. Previously, parameter objects such as \code{SelectParams} and \code{TrainParams} had to be explicitly specified, making it challenging for users not familar with S4 object-oriented programming.
  \item
  Basic multi-omics data integration functionality available via \code{crossValidate} which allows combination of different tables. Pre-validation and PCA dimensionality techniques provide a fair way to compare high-dimensional omics data with low-dimensional clinical data. Also, it is possible to simply concatenate all data tables.
  \item
  Model-agnostic variable importance calculated by training when leaving out one selected variable at a time. Turned off by default as it substantially increases run time. See \code{doImportance} parameter of \code{ModellingParams} for more details.
  \item
  Parameters specifying the cross-validation procedure and data modelling formalised as \code{CrossValParams} and \code{ModellingParams} classes.
  \item
  Feature selection can now be done either based a on resubstitution metric (i.e. train and test on the training data) or a cross-validation metric (i.e. split the training data into training and testing partitions to tune the selected features). all feature selection functions have been converted into feature ranking functions, because the selection procedure is a feature of cross-validation.
  \item
  All function and class documentation coverted from manually written Rd files to Roxygen format.
  \item
  Human Reference Interactome (binary experimental PPI) included in bundled data for pairs-based classification. See \code{?HuRI} for more details.
  \item
  Performance plots can now do either box plots or violin plots. Box plot remains the default style.  
}}

\section{Version 2.14.0}{\itemize{
  \item
  Upsampling and downsampling to equalise class sizes added.
}}

\section{Version 2.6.0}{\itemize{
  \item
  Two-stage easy-hard classifier added.
}}

\section{Version 2.2.0}{\itemize{
  \item
  getClasses is no longer a slot of PredictParams. Every predictor function needs to return either a
  factor vector of classes, a numeric vector of class scores for the second class, or a data frame with
  a column for the predicted classes and another for the second-class scores.
  \item
  Cross-validations which use folds ensure that samples belonging to each class are in approximately
  the same proportions as they are for the entire data set.
  \item
  Classification can reuse fitted model from previous classification by using previousTrained function.  
  \item
  Feature selection using gene sets and networks. Classification can use meta-features derived
  from the individual features used for feature selection.
  \item
  tTestSelection function for feature selection based on ordinary t-test statistic ranking. Now the default
  feature selection function, if none is specified.
  \item
  Tuning parameter optimisation metric is specified by providing a tuneOptimise parameter to TrainParams rather than depending on ResubstituteParams being used during feature selection.
}}

\section{Version 2.0.0}{\itemize{
  \item
  Broad support for DataFrame and MultiAssayExperiment data sets by feature selection and classification functions.
  \item
  The majority of processing is now done in the DataFrame method for functions that implement methods for
  multiple kinds of inputs.
  \item
  Elastic net GLM classifier and multinomial logistic regression classifier wrapper functions.
  \item
  Plotting functions have a new default style using a white background with black axes.
  \item
  Vignette simplified and uses a new mass cytometry data set with clearer differences between classes to
  demonstrate classification and its performance evaluation.
}}

\section{Version 1.12.0}{\itemize{
  \item
  Alterations to make plots compatible with ggplot versions 2.2 and greater.
  \item
  calcPerformance can calculate some performance metrics for classification tasks based on data sets
  with more than two classes.
  \item
  Sample-wise metrics, like sample-specific error rate and sample-specific accuracy are calculated by
  calcPerformance and added to the ClassifyResult object, rather than by samplesMetricMap and being
  inaccessible to the end-user.
}}

\section{Version 1.10.0}{\itemize{
  \item
  errorMap replaced by samplesMetricMap. The plot can now show either error rate or accuracy.
}}

\section{Version 1.8.0}{\itemize{
  \item
  Ordinary k-fold cross-validation option added.
  \item
  Absolute difference of group medians feature selection function added.
}}

\section{Version 1.4.0}{\itemize{
  \item
  Weighted voting mode that uses the distance from an observation to the nearest crossover point
  of the class densities added.
  \item
  Bartlett Test selection function included.
  \item
  New class SelectResult. rankPlot and selectionPlot can additionally work with lists
  of SelectResult objects. All feature selection functions now return a SelectResult object
  or a list of them.
  \item
  priorSelection is a new selection function for using features selected in a prior cross
  validation for a new data set classification.
  \item
  New weighted voting mode, where the weight is the distance of the x value from the
  nearest crossover point of the two densities. Useful for predictions with skewed 
  features.
}}

\section{Version 1.2.0}{\itemize{
  \item
  More classification flexibility, now with parameter tuning integrated into the process.
  \item
  New performance evaluation functions, such as a ROC curve and a performance plot.
  \item
  Some existing predictor functions are able to return class scores, not just class labels.
}}

\section{Version 1.0.0}{\itemize{
  \item
  First release of the package, which allows parallelised and customised classification, with many convenient performance evaluation functions.
}}