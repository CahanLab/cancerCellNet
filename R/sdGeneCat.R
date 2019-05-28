#' @title
#' Standard Deviation Category Gene Expressions
#' @description
#' Finding the standard deviation the gene expression level for each category
#'
#' @param expDat the normalized expression table from \code{\link{trans_prop}}
#' @param sampTab the sample table
#' @param dLevel the name of the column that contains categories
#' @param sampID the name of the column that contains sample ID
#'
#' @return a matrix of the standard deviation for gene expressions
#' @export
sdGeneCat<-function(expDat, sampTab, dLevel, sampID = NULL){

  if (is.null(sampID) == TRUE) {
    sampID = "sampID"

    sampTab[, sampID] = rownames(sampTab)
  }

  stdMatrix = matrix(nrow=nrow(expDat), ncol=0)
  rownames(stdMatrix) = rownames(expDat)


  for (cat in unique(sampTab[, dLevel])) {

    tempExpDat = expDat[, as.vector(sampTab[sampTab[, dLevel] == cat, sampID])]

    #calculates the standard deviation
    tempStd = matrix(apply(tempExpDat, 1, sd), ncol = 1, nrow = length(apply(tempExpDat, 1, sd)))
    rownames(tempStd) = names(apply(tempExpDat, 1, sd))
    colnames(tempStd) = paste0(cat, "_Std")

    if(all(rownames(stdMatrix) == rownames(tempStd))) {
      stdMatrix = cbind(stdMatrix, tempStd)
    }
  }

  #return
  stdMatrix
}
