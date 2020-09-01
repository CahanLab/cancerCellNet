#' @title
#' Predict sub class query using broad class and subclass classifier (newer version)
#' @description
#' The function predicts the query data using the broad class and subclass classifier
#' @param cnProc_broad the cnProc of the broad classifier
#' @param cnProc_subclass the cnProc of the subclass classifier
#' @param expDat the expression data of query data
#' @param nrand the number of random profiles generate for evaluation process
#' @param weight_broadClass the weight that will be put on the broadclass classification
#' @return a classification score matrix
#'
#' @export

subClass_predict<-function(cnProc_broad, cnProc_sub, expDat, nrand = 2, weight_broadClass = 1) {
   ccnList = cnProc_sub # only concerns with subclass classifier

   rf_tsp = ccnList[['classifier']]
   cgenes = ccnList[['cgenes']]
   xpairs = ccnList[['xpairs']]

   classMatrix = broadClass_predict(cnProc_broad, expDat = expDat, nrand = nrand)
   classMatrix = classMatrix[, -grep("rand", colnames(classMatrix))]

  if (weight_broadClass > 1) {
    print("Time to add weights")
    originalRowNames = rownames(classMatrix)
    originalClassMatrix = classMatrix
    weightIter = c(1:(weight_broadClass - 1))

    for (weight in weightIter) {
      newRownames = paste0(originalRowNames, "-", weight)
      additionalClassMatrix = originalClassMatrix
      rownames(additionalClassMatrix) = newRownames
      classMatrix = rbind(classMatrix, additionalClassMatrix)
    }

    print("finished adding weights")
  }

   expValTrans = subClassQuery_transform(expDat = expDat, cgenes = cgenes, xpairs = xpairs, classMatrix = classMatrix)

   returnMatrix = rf_classPredict(rf_tsp, expValTrans, numRand=nrand)

   return(returnMatrix)
}
