################################
#Purpose: Calcuate Score
#Authors: Pichai Raman
#Date: 8/5/2019
################################

#' @title Calc score
#' @author Pichai Raman
#' @author Sherjeel Arif
#' @description  Function to calculate the gene ratio for each sample.
#' @details This function takes 2 input parameters. The first parameter is a gene ratio matrix
#' in which the entries are the expression values for each gene signature. The second parameter
#' is a list of 4, where each list has the gene ratios associated with the 4 Medulloblastoma
#' subtypes.
#'
#' The function takes the gene ratio matrix and computes the column mean. However, function does
#' not use all the values in the column. Only the values of the gene ratios associated with a given
#' Medulloblastoma subtype are used to compute a column mean.
#'
#' The expected output is a matrix data frame where the rownames are the Sample Identifiers and
#' the column names are the 4 Medulloblastoma subtypes.
#'
#' @param myMat matrix containing gene ratios for each gene signature and sample. The row names
#' contain the gene signatures and the column names contain the individual samples.
#' @param mySetsUp list of 4 containing the gene ratios associated with
#' each of the 4 molecular subtypes of Medulloblastoma. Provided in the data directory as
#' 'medulloSetsUp'.
#'
#' @examples
#' ## load provided data
#' data(geneRatOut_109401)
#'
#' ## Use calcScore function
#' myMat <- calcScore(geneRatOut_109401)
#'
#' ## View contents of matrix
#' head(myMat[1:4])
#'
#' @export
#'

calcScore <- function(myMat=NULL, mySetsUp=NULL) {
  # if no value is supplied, use default data
  if(is.null(mySetsUp)) {
    mySetsUp <- get(utils::data("medulloSetsUp"))
  }

  getScoreSet <- function(x, myMat=myMat) {
    return(colMeans(myMat[x,]))
  }

  myMatUp <- data.frame(lapply(mySetsUp,FUN=getScoreSet, myMat));
  return(myMatUp)
}


