#' @title Calc Stats
#' @author Pichai Raman
#' @author Sherjeel Arif
#' @name calcStats
#' @description Function to calculate statistics on prediction class and
#' actual class of a medulloblastoma dataset.
#' @details This function can be used to analyze prediction class if the correct
#' classification of the provided Medulloblastoma dataset is known. Function
#' calculates an accuracy score, a confidence interval, creates a confusion
#' matrix, and displays other useful statistics (i.e. Sensitivity, Specificity).
#'
#' The parameters are character vectors and each value in the vector corresponds to
#' one of the four molecular subtypes of Medulloblastoma (i.e. Sonic Hedgehog (SHH),
#' WNT, Group 3, and Group 4).
#'
#' @param myClassActual a vector of characters corresponding to the actual
#' classification for each sample in the dataset.
#' @param myClassPred a vector of characters corresponding to the prediction made
#' for each sample.
#'
#' @examples
#' ## mock prediction class
#' predictionClass <- c("Group3", "WNT", "WNT", "Group4", "SHH", "Group4")
#'
#' ## mock actual class
#' actualClass <- c("Group4", "WNT", "SHH", "Group4", "SHH", "Group4")
#'
#' ## call function to obtain accuracy and other statistics
#' medulloPackage::calcStats(actualClass, predictionClass)
#'
#' @export


#######################
# Required Libraries
#######################

library(pheatmap)
library(GSVA)
library(caret)
library(lattice)
library(preprocessCore)
source("~/medulloPackage/R/calcScore.R")


calcStats <- function(myClassActual = NULL, myClassPred = NULL) {
  x <- sum(lengths(regmatches(myClassActual, gregexpr("U", myClassActual)))) # counts "U"'s : the number of unknowns
  myScore <- sum(myClassPred==myClassActual)/(length(myClassActual)-sum(x)) # calculate accuracy score

  sampAnnot <- data.frame(myClassPred, myClassActual);
  sampAnnot[,"Correct"] <- myClassPred==myClassActual
  sampAnnot <- sampAnnot[sampAnnot[,2]!="U",]
  sampAnnot[,2] <- factor(sampAnnot[,2], levels=c("Group3", "Group4", "WNT", "SHH"))
  sampAnnot[,1] <- factor(sampAnnot[,1], levels=c("Group3", "Group4", "WNT", "SHH"))

  print(caret::confusionMatrix(sampAnnot[,1], sampAnnot[,2]));

  writeLines("")
  print(paste("Accuracy: ", myScore))
}
