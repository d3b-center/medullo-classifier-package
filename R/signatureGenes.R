#' @title Signature Genes
#' @name signatureGenes
#' @author Pichai Raman
#' @author Komal S. Rathi
#' @author Sherjeel Arif
#' @description  Function to calculate gene ratio matrix
#' @details Function for determining the gene ratio matrix. The function expects an expression matrix
#' as the only necessary input parameter. The expression matrix must be one of the following formats:
#' (1) FPKM (2) TPM (3)quantile normalized data (4) microarray data.
#'
#' The expected output is a matrix data frame in which the row names are the gene ratio symbols and
#' the column names are the sample identifiers.
#'
#' @param exprs expression matrix in which the rows contains the HUGO/HGNC gene symbols and the
#' columns are the samples identifiers. The entries inside the matrix are the respective gene
#' expression values. (see details section for more info)
#' @param signatureProbes character vector of gene signatures
#' @examples
#' # Load example data containing expression matrix
#' data(exprs_109401)
#'
#' # use function
#' geneRatioOut_109401 <- signatureGenes(exprs_109401)
#'
#' # View gene ratio matrix
#' head(geneRatioOut_109401[1:5])
#'
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
