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
#' @return a list containing average matrix and standard matrix
#' @export
avgGeneCat<-function(expDat, sampTab, dLevel, sampID){

   returnList = list()

   returnMatrix = matrix(nrow=nrow(expDat), ncol=0)
   rownames(returnMatrix) = rownames(expDat)

   stdMatrix = matrix(nrow=nrow(expDat), ncol=0)
   rownames(stdMatrix) = rownames(expDat)


   for (cat in unique(sampTab[, dLevel])) {

      tempExpDat = expDat[, sampTab[sampTab[, dLevel] == cat, sampID]]

      # calculates the mean
      tempMatrix = matrix(apply(tempExpDat, 1, mean), ncol = 1, nrow = length(apply(tempExpDat, 1, mean)))
      rownames(tempMatrix) = names(apply(tempExpDat, 1, mean))
      colnames(tempMatrix) = paste0(cat, "_Avg")

      #calculates the standard deviation
      tempStd = matrix(apply(tempExpDat, 1, sd), ncol = 1, nrow = length(apply(tempExpDat, 1, sd)))
      rownames(tempStd) = names(apply(tempExpDat, 1, sd))
      colnames(tempStd) = paste0(cat, "_Std")

      if(all(rownames(returnMatrix) == rownames(tempMatrix))) {
         returnMatrix = cbind(returnMatrix, tempMatrix)
      }

      if(all(rownames(stdMatrix) == rownames(tempStd))) {
        stdMatrix = cbind(stdMatrix, tempStd)
      }
   }

   returnList[["Avg_Matrix"]] = returnMatrix
   returnList[["Std_Matrix"]] = stdMatrix
   #return
   returnList
}

