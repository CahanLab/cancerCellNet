#' @title
#' Average Category Gene Expressions
#' @description
#' Averaging the gene expression level for each category
#'
#' @param expDat the normalized expression table from \code{\link{trans_prop}}
#' @param sampTab the sample table
#' @param dLevel the name of the column that contains categories
#' @param sampID the name of the column that contains sample ID
#'
#' @return an average matrix of the gene expressions
#' @export
avgGeneCat<-function(expDat, sampTab, dLevel, sampID = NULL){

  if (is.null(sampID) == TRUE) {
    sampID = "sampID"

    sampTab[, sampID] = rownames(sampTab)
  }

   returnMatrix = matrix(nrow=nrow(expDat), ncol=0)
   rownames(returnMatrix) = rownames(expDat)

   for (cat in unique(sampTab[, dLevel])) {

      tempExpDat = expDat[, as.vector(sampTab[sampTab[, dLevel] == cat, sampID])]

      # calculates the mean
      tempMatrix = matrix(apply(tempExpDat, 1, mean), ncol = 1, nrow = length(apply(tempExpDat, 1, mean)))
      rownames(tempMatrix) = names(apply(tempExpDat, 1, mean))
      colnames(tempMatrix) = paste0(cat, "_Avg")

      if(all(rownames(returnMatrix) == rownames(tempMatrix))) {
         returnMatrix = cbind(returnMatrix, tempMatrix)
      }
   }

   #return
   returnMatrix
}

