################################
#Purpose: Calcuate Score
#Authors: Pichai Raman
#Date: 8/5/2019
################################

#' @title Calc score
#' @author Pichai Raman
#' @author Komal S. Rathi
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
#' The expected output is a list of two dataframes:
#' one is a data frame where the rownames are the Sample Identifiers and the column names are the 4 Medulloblastoma subtypes.
#' and the second is a data frame of t-test pvalues across all combination of subtypes per samples
#'
#' @param myMat matrix containing gene ratios for each gene signature and sample. The row names
#' contain the gene signatures and the column names contain the individual samples.
#' @param mySetsUp list of 4 containing the gene ratios associated with
#' each of the 4 molecular subtypes of Medulloblastoma. Provided in the data directory as
#' 'medulloSetsUp'.
#'
#' @examples
#' ## load provided data
#' data(geneRatioOut_109401)
#'
#' ## Use calcScore function
#' myMat <- calcScore(geneRatioOut_109401)
#'
#' ## View contents of matrix
#'
#' # mean scores
#' head(myMat[[1]][1:4])
#'
#' # t-test pvalues
#' head(myMat[[2]][1:4])
#'
#' @export
#'

calcScore <- function(myMat = NULL, mySetsUp = NULL) {
  # if no value is supplied, use default data
  if(is.null(mySetsUp)) {
    mySetsUp <- get(utils::data("medulloSetsUp"))
  }

  getScoreSet <- function(x, myMat = myMat) {
    return(colMeans(myMat[x,]))
  }

  myMatUp <- data.frame(lapply(mySetsUp, FUN = getScoreSet, myMat))

  # for each column, do a t-test across all combinations of medullo subtypes
  calc.pvalues <- function(x, myMat, mySetsUp) {

    # create a dataframe of four subtypes per sample
    wnt = myMat[mySetsUp$WNT, x]
    ssh = myMat[mySetsUp$SHH, x]
    gr3 = myMat[mySetsUp$Group3, x]
    gr4 = myMat[mySetsUp$Group4, x]
    df <- rowr::cbind.fill(wnt, ssh, gr3, gr4, fill = NA)
    colnames(df) <- c("WNT","SHH","Group3","Group4")

    # compute t-test on all combinations of subtypes
    combinations <- utils::combn(colnames(df), 2, simplify = FALSE)
    results <- lapply(seq_along(combinations), function (n) {
      df <- df[,colnames(df) %in% unlist(combinations[n])]
      result <- stats::t.test(df[,1], df[,2])$p.value
      return(result)})

    df1 <- data.frame(sample = x,
                      one = matrix(unlist(combinations), ncol = 2, byrow = TRUE)[,1],
                      two = matrix(unlist(combinations), ncol = 2, byrow = TRUE)[,2],
                      p.value = unlist(results))

    df2 <- data.frame(sample = x,
                      one = matrix(unlist(combinations), ncol = 2, byrow = TRUE)[,2],
                      two = matrix(unlist(combinations), ncol = 2, byrow = TRUE)[,1],
                      p.value = unlist(results))

    results <- unique(rbind(df1, df2))
    results$p.value[is.na(results$p.value)] <- 0 # added
    results <- plyr::ddply(.data = results, .variables = "one", .fun = function(x) x[which(x$p.value == max(x$p.value)),])
  }

  # apply function on all columns
  pval.list <- lapply(colnames(myMat), FUN = function(x) calc.pvalues(x, myMat, mySetsUp))
  pval.list <- do.call(rbind, pval.list)
  pval.list <- unique(pval.list[,c("sample","one","p.value")]) # added

  return(list(pred = myMatUp, pval = pval.list))
}


