\name{ClassifyRnews}
\title{ClassifyR News}
\encoding{UTF-8}

\section{Version 2.2.0}{\itemize{
  \item
  Feature selection using gene sets and networks. Classification can use meta-features derived
  from the individual features used for feature selection.
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
  Vignette simplified and uses a new mass cytometry dataset with clearer differences between classes to
  demonstrate classification and its performance evaluation.
}}

\section{Version 1.12.0}{\itemize{
  \item
  Alterations to make plots compatible with ggplot versions 2.2 and greater.
  \item
  calcPerformance can calculate some performance metrics for classification tasks based on datasets
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
  validation for a new dataset classification.
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