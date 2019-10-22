################################
#Purpose: Classify Medullo Samples
#Authors: Pichai Raman, Sherjeel Arif
#Date: 8/5/2019
################################

#' @title Classify
#' @name classify
#' @author Pichai Raman
#' @author Sherjeel Arif
#' @author Komal Rathi
#' @description  Function to Classify samples into Medulloblastoma subtypes using
#' transcriptomic data.
#' @details Classifier for predicting amongst the 4 molecular subtypes of
#' Medulloblastoma, Sonic Hedgehog (SHH), WNT, Group 3, and Group 4 from
#' RNA-Seq or microarray data. The input is an expression matrix. The following types of data
#' are allowed as inputs: (1) FPKM (2) TPM (3) quantile normalized data (4) microarray data. The
#' expected output of this function is a character vector containing the Medulloblastoma subtypes
#' (i.e. 'WNT', 'SHH', 'Group3', 'Group4').
#' @param exprs matrix containing gene expression values.
#' The row names contain HUGO/HGNC gene symbols and the column names contain the sample identifiers.
#' @param medulloGeneSetsUp list of 4 containing the gene signature associated with
#' each of the 4 molecular subtypes of Medulloblastoma.
#' @examples
#' # Load example data containing expression matrix
#' data(exprs_109401)
#'
#' # use classification function on expression matrix
#' prediction <- medulloPackage::classify(exprs_109401)
#'
#' # View classifier output
#' head(prediction)
#'
#' @export
#'

classify <- function(exprs = NULL, medulloGeneSetsUp = NULL) {

  # if no value is supplied, use default data
  if(is.null(medulloGeneSetsUp)) {
    medulloGeneSetsUp <- get(utils::data("medulloSetsUp"))
  }

  # calculate gene ratio matrix
  geneRatioOut <- signatureGenes(exprs)


  medulloGeneSetsUp$WNT <- intersect(medulloGeneSetsUp$WNT, rownames(geneRatioOut))
  medulloGeneSetsUp$SHH <- intersect(medulloGeneSetsUp$SHH, rownames(geneRatioOut))
  medulloGeneSetsUp$Group3 <- intersect(medulloGeneSetsUp$Group3, rownames(geneRatioOut))
  medulloGeneSetsUp$Group4 <- intersect(medulloGeneSetsUp$Group4, rownames(geneRatioOut))

  myMat <- calcScore(geneRatioOut, medulloGeneSetsUp);
  myClassPred <- colnames(myMat)[max.col(myMat,ties.method="first")]
  return(myClassPred)
}
