#' @title
#' Predict query using broad class classifier
#' @description
#' The function predicts the query data using the broad class classifier
#' @param cnProc the cnProc of the broad classifier
#' @param expDat the expression data of query data
#' @param nrand the number of random profiles generate for evaluation process
#' @return a classification score matrix
#'
#' @export
broadClass_predict<-function(cnProc, expDat, nrand = 2) {
   ccnList = cnProc
   expVal = expDat

   rf_tsp = ccnList[['classifier']]
   cgenes = ccnList[['cgenes']]
   xpairs = ccnList[['xpairs']]

   cat("Loaded in the cnProc\n")

   expValTrans = query_transform(expVal[cgenes,], xpairs)
   classRes_val = rf_classPredict(rf_tsp, expValTrans, numRand=nrand)

   cat("All Done\n")

   #return
   classRes_val
}
