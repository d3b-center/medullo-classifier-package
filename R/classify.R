################################
#Purpose: Classify Medullo Samples
#Authors: Pichai Raman, Sherjeel Arif
#Date: 8/5/2019
################################

#' @title Classify
#' @name classify
#' @author Pichai Raman
#' @author Sherjeel Arif
#' @description  Function to classify
#' @details Classifier for predicting amongst the 4 molecular subtypes of
#' Medulloblastoma, Sonic Hedgehog (SHH), WNT, Group 3, and Group 4 from
#' RNA-Seq or microarray data.
#' @param geneRatioOut matrix containing gene ratios for each gene signature and sample.
#' The row names contain the gene signatures and the column names contain the individual samples.
#' @param medulloGeneSetsUp list of 4 containing the gene signature associated with
#' each of the 4 molecular subtypes of Medulloblastoma.
#' @export
#'

library(reshape2)
source("~/medulloPackage/R/signatureGenes.R")

classify <- function(exprs = NULL)
{
  medulloGeneSetsUp <- readRDS("data/medulloSetsUp.RDS")

  # calculate gene ratio matrix
  geneRatioOut <- signatureGenes(exprs)


  medulloGeneSetsUp$WNT <- intersect(medulloGeneSetsUp$WNT, rownames(geneRatioOut))
  medulloGeneSetsUp$SHH <- intersect(medulloGeneSetsUp$SHH, rownames(geneRatioOut))
  medulloGeneSetsUp$Group3 <- intersect(medulloGeneSetsUp$Group3, rownames(geneRatioOut))
  medulloGeneSetsUp$Group4 <- intersect(medulloGeneSetsUp$Group4, rownames(geneRatioOut))

  myMat <- calcScoreMat(geneRatioOut, medulloGeneSetsUp);
  myClassPred <- colnames(myMat)[max.col(myMat,ties.method="first")]
  return(myClassPred)
}
