################################
#Purpose: Calcuate Score
#Authors: Pichai Raman
#Date: 8/5/2019
################################

#' @title Calc score
#' @author Pichai Raman
#' @description  Function to calculate the gene ratio for each sample.
#' @details None
#' @param myMat matrix containing gene ratios for each gene signature and sample. The row names
#' contain the gene signatures and the column names contain the individual samples.
#' @param mySetsUp list of 4 containing the gene ratios associated with
#' each of the 4 molecular subtypes of Medulloblastoma. Provided in the data directory.
#' @export
#'

calcScoreMat <- function(myMat=NULL, mySetsUp=NULL) {
  getScoreSet <- function(x, myMat=myMat) {
    return(colMeans(myMat[x,]))
  }

  myMatUp <- data.frame(lapply(mySetsUp,FUN=getScoreSet, myMat));
  return(myMatUp)
}


