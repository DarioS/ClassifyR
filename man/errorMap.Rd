\name{errorMap}
\alias{errorMap}
\alias{errorMap,list-method}
\title{Plot a Grid of Sample Error Rates}
\description{A grid of coloured tiles is drawn. There is one column for each sample
and one row for each classification result.
       }
\usage{
  \S4method{errorMap}{list}(results,
                     errorColours = c("#FFFFFF", "#E0E0E0", "#BABABA", "#888888", "#000000"),
                     classColours = c("blue", "red"), fontSizes = c(24, 16, 12, 12, 12),
                     mapHeight = 4)
}
\arguments{
  \item{results}{A list of \code{\link{ClassifyResult}} objects.}
  \item{errorColours}{A vector of colours for error levels.}
  \item{classColours}{A vector of colours for class levels.}
  \item{fontSizes}{A vector of length 5. The first number is the size of the title.
                   The second number is the size of the axes titles. The third number is
                   the size of the axes values. The fourth number is the size of the
                   legends' titles. The fifth number is the font size of the legend labels.}
  \item{mapHeight}{Height of the map, relative to the height of the class colour bar.}                 
}
\details{
The names of \code{results} determine the row names that will be in the plot.
The length of \code{errorColours} determines how many bins the error rates will
be discretised to.
}
\value{
  A plot is produced and a grob is returned that can be saved to a graphics device.
}
\author{Dario Strbenac}

\examples{
  predicted <- data.frame(sample = sample(10, 100, replace = TRUE),
                            predicted = rep(c("Healthy", "Cancer"), each = 50))
  actual <- factor(rep(c("Healthy", "Cancer"), each = 5))
  result1 <- ClassifyResult(LETTERS[1:10], LETTERS[10:1], list(sample(10, 10)),
                            list(predicted), actual, list("fold", 100, 5))
  predicted[, "predicted"] <- sample(predicted[, "predicted"])
  result2 <- ClassifyResult(LETTERS[1:10], LETTERS[10:1], list(sample(10, 10)),
                            list(predicted), actual, validation = list("leave", 1))
  # wholePlot <- errorMap(list(Gene = result1, Protein = result2)) # Wait for namespace problems to be fixed.
  # ggsave("wholePlot.png", wholePlot)
}