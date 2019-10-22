#' @title Calc Stats
#' @author Pichai Raman
#' @author Sherjeel Arif
#' @author Komal Rathi
#' @name calcStats
#' @description Function to calculate statistics on prediction class and
#' actual class of a medulloblastoma dataset.
#' @details This function can be used to analyze prediction class if the correct
#' classification of the provided Medulloblastoma dataset is known. Function
#' calculates an accuracy score, a confidence interval, creates a confusion
#' matrix, and displays other useful statistics (i.e. Sensitivity, Specificity).
#'
#' The input parameters are character vectors and each value in the vector corresponds to
#' one of the four molecular subtypes of Medulloblastoma (i.e. Sonic Hedgehog (SHH),
#' WNT, Group 3, and Group 4).
#'
#' The expected output is a list of length 4. The first element
#' in the outputted list is the confusion matrix, which is a data frame of dimension 4x4.
#' The second element is a data frame containing overall statistics on the accuracy of the classifier.
#' The third element in the the list is also a data frame containing class statistics (e.g.
#' Sensitivity, Specificity, etc.). The last element in the list is a character which corresponds
#' to the the prediction accuracy of the medulloblastoma classifier.
#'
#' @param myClassActual a vector of characters corresponding to the actual
#' classification for each sample in the dataset. The vector can only include the characters
#' 'WNT', SHH', 'Group3', 'Group4', and 'U'.
#' @param myClassPred a vector of characters corresponding to the prediction made
#' for each sample. The vector can only include the characters
#' 'WNT', SHH', 'Group3', and 'Group4'.
#'
#' @examples
#' ## prediction class
#' predictionClass <- c("Group3", "WNT", "WNT", "Group4", "SHH", "Group4")
#'
#' ## observed class
#' actualClass <- c("Group4", "WNT", "SHH", "Group4", "SHH", "Group4")
#'
#' ## call function to obtain accuracy and other statistics
#' stats <- medulloPackage::calcStats(actualClass, predictionClass)
#'
#' ## View the contents of stats
#' confusion.matrix <- stats[[1]]
#' overall.stats <- stats[[2]]
#' class.stats <- stats[[3]]
#' accuracy <- stats[[4]]
#'
#' head(confusion.matrix)
#' head(overall.stats)
#' head(class.stats)
#' head(accuracy)
#'
#' @export


#######################
# Required Libraries
#######################

 calcStats <- function(myClassActual = NULL, myClassPred = NULL) {
  x <- sum(lengths(regmatches(myClassActual, gregexpr("U", myClassActual)))) # counts "U"'s : the number of unknowns
  myScore <- sum(myClassPred==myClassActual)/(length(myClassActual)-sum(x)) # calculate accuracy score
  myScore <- format(myScore*100, digits = 2)

  sampAnnot <- data.frame(myClassPred, myClassActual);
  sampAnnot[,"Correct"] <- myClassPred==myClassActual
  sampAnnot <- sampAnnot[sampAnnot[,2]!="U",]
  sampAnnot[,2] <- factor(sampAnnot[,2], levels=c("Group3", "Group4", "WNT", "SHH"))
  sampAnnot[,1] <- factor(sampAnnot[,1], levels=c("Group3", "Group4", "WNT", "SHH"))

  res <- caret::confusionMatrix(sampAnnot[,1], sampAnnot[,2])

  # format confusion matrix
  confusion.matrix <- res$table
  confusion.matrix <- as.data.frame.matrix(confusion.matrix)
  rownames(confusion.matrix) <- paste0('Pred_', rownames(confusion.matrix))
  colnames(confusion.matrix) <- paste0('Ref_', colnames(confusion.matrix))

  # overall stats
  overall.stats <- data.frame(stats = res$overall)
  overall.stats <- format(overall.stats, digits =3)
  overall.stats[1:5,] <- format(paste0(as.numeric(overall.stats[1:5,])*100,"%"), 2)

  # stats by class
  class.stats <- as.data.frame(res$byClass)
  class.stats <- format(class.stats, digits = 3)
  class.stats[] <- lapply(class.stats, function(x) as.numeric(as.character(x)))
  class.stats[] <- lapply(class.stats, function(x) paste0(x*100, "%"))

  # accuracy in %
  acc <- paste0("Accuracy: ", myScore, "%")

  res <- list(confusion.matrix, overall.stats, class.stats, acc)

  print(acc)
  return(res)
}
