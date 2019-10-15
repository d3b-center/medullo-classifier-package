#' @title Signature Genes
#' @name signatureGenes
#' @author Pichai Raman
#' @author Sherjeel Arif
#' @author Komal Rathi
#' @description  Function to calculate gene ratio matrix
#' @details Function for determining the gene ratio matrix. The function expects an expression matrix
#' as the only parameter.
#' @param exprs expression matrix in which the rows contains the gene symbols and the
#' columns are the samples. The values inside the matrix are the respective gene signatures.
#' @param signatureProbes character vector of gene signature
#' @export
#'

signatureGenes <- function(exprs = NULL, signatureProbes = NULL) {
  ################################
  #Now read in signature genes
  #Filter matrix to signature genes
  #and Create Gene Ratios
  ################################

  # if no value is supplied, use default data
  if(is.null(signatureProbes)) {
    signatureProbes <- get(utils::data("bestFeaturesNew"))
  }

  getGenes <- function(x) {
    out <- strsplit(x, split="_")
    output <- c(out[[1]])
  }

  signatureGenes <- sapply(signatureProbes, FUN=getGenes)
  signatureGenes <- as.character(signatureGenes);
  signatureGenes <- unique(signatureGenes);

  #Filter Matrix
  exprs_SG <- exprs[intersect(rownames(exprs), signatureGenes), ]

  #Create Ratios
  createRatio <- function(x) {
    g1 <- x[1];
    g2 <- x[2];
    g1g2_ratio <- 2^(exprs[g1,]-exprs[g2,])
    return(g1g2_ratio)
  }

  corGenes <- stats::cor(t(exprs_SG));
  corGenes <- data.frame(reshape2::melt(corGenes));
  corGenes <- corGenes[corGenes[,"value"]<.99,];

  exprs <- as.matrix(exprs);
  geneRatioOut <- apply(corGenes, FUN=createRatio, MARGIN=1)
  geneRatioOut <- data.frame(t(geneRatioOut));

  rownames(geneRatioOut) <- paste(corGenes[,1], corGenes[,2], sep="_");
  colnames(geneRatioOut) <- colnames(exprs)

  ################################
  #Filter to signature ratios
  #Create Heatmap
  ################################

  geneRatioOut <- geneRatioOut[intersect(rownames(geneRatioOut), signatureProbes),]
  return(geneRatioOut)

}
