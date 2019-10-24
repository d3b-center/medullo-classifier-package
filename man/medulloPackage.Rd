\docType{package}
\name{medulloPackage}
\alias{medulloPackage}
\title{MedulloClassifier: R Package to classify medulloblastoma datasets}
\description{
The goal of this package is to develop a classifier that can accurately
predict amongst 4 molecular subtypes of Medulloblastoma, Sonic Hedgehog (SHH), WNT, Group 3,
and Group 4 from RNA-Seq or microarray data (and potentially any transcriptomic data).
}
\details{
These are the functions available in the package:

  \code{\link{classify}}: returns a character vector containing the predictions made by the
  classifier (i.e. 'WNT', 'SHH', 'Group3', 'Group4').
  \cr
  \cr
  \code{\link{calcScore}}: Function to calculate the column mean for each sample with respect to the   gene signatures associated with each of the 4 Medulloblastoma subtypes.
  \cr
  \cr
  \code{\link{calcStats}}: Function to analyze the prediction class given that the correct
  classification of the provided Medulloblastoma dataset is known. Function
  calculates an accuracy score, a confidence interval, creates a confusion
  matrix, and displays other useful statistics (i.e. Sensitivity, Specificity).
  \cr
  \cr
  \code{\link{signatureGenes}}: Function that calculates the gene ratios for each gene signature.
}
\examples{
## Following is an example of how to use the functions
## in this package:

library(medulloPackage)

## Load data containing example expression matrix
data(exprs_109401)

## Load vector of observed classifications
data(actual_109401)

## use classification function on expression matrix
pred_109401 <- medulloPackage::classify(exprs_109401)
pred_109401.classes <- pred_109401$best.fit # get the predicted classes
pred_109401.pvals <- pred_109401$p.value # get associated p-values

## View output of classifier
head(pred_109401.classes)

## Calculate statistics
stats <- medulloPackage::calcStats(actual_109401, pred_109401.classes)

## View stats
confusion.matrix <- stats[[1]]
overall.stats <- stats[[2]]
class.stats <- stats[[3]]
accuracy <- stats[[4]]

## Print confusion matrix
head(confusion.matrix)

## Prints accuracy stats
head(overall.stats)

## Prints sensitivity, specificity, etc. for each class
head(class.stats)

## Prints overall accuracy of classifier in making correct predictions
head(accuracy)

}

