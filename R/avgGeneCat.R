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
#' @return a matrix containing the average of
#' @export
avgGeneCat<-function(expDat, sampTab, dLevel, sampID){

   returnMatrix = matrix(nrow=nrow(expDat), ncol=0)
   rownames(returnMatrix) = rownames(expDat)

   for (cat in unique(sampTab[, dLevel])) {

      tempExpDat = expDat[, sampTab[sampTab[, dLevel] == cat, sampID]]

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

