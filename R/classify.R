################################
#Purpose: Classify Medullo Samples
#Authors: Pichai Raman, Sherjeel Arif
#Date: 8/5/2019
################################

#' @title Classify
#' @author Pichai Raman
#' @author Sherjeel Arif
#' @description  Function to classify
#' @details Classifier for predicting amongst the 4 molecular subtypes of
#' Medulloblastoma, Sonic Hedgehog (SHH), WNT, Group 3, and Group 4 from
#' RNA-Seq or microarray data.
#' @param geneRatioOut Matrix containing the gene ratios for different samples.
#' @param medulloGeneSetsUp list of 4 containing the gene ratios associated with
#' each of the 4 molecular subtypes of Medulloblastoma.
#' @export
#'

classify <- function(geneRatioOut = NULL, medulloGeneSetsUp="../../data/medulloSetsUp.RDS")
{
  medulloGeneSetsUp$WNT <- intersect(medulloGeneSetsUp$WNT, rownames(geneRatioOut))
  medulloGeneSetsUp$SHH <- intersect(medulloGeneSetsUp$SHH, rownames(geneRatioOut))
  medulloGeneSetsUp$Group3 <- intersect(medulloGeneSetsUp$Group3, rownames(geneRatioOut))
  medulloGeneSetsUp$Group4 <- intersect(medulloGeneSetsUp$Group4, rownames(geneRatioOut))

  myMat <- calcScoreMat(geneRatioOut, medulloGeneSetsUp);
  myClassPred <- colnames(myMat)[max.col(myMat,ties.method="first")]
  return(myClassPred)

}


