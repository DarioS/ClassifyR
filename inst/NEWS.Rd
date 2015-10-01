\name{ClassifyRnews}
\title{ClassifyR News}
\encoding{UTF-8}

\section{Version 1.4.0}{\itemize{

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