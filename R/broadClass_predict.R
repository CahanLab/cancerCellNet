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

   rf_tsp<-ccnList[['classifier']]
   cgenes<-ccnList[['cgenes']]
   xpairs<-ccnList[['xpairs']]

   print("Loaded in the cnProc")

   expValTrans<-query_transform(expVal[cgenes,], xpairs)
   classRes_val<-rf_classPredict(rf_tsp, expValTrans, numRand=nrand)

   print("All Done")

   #return 
   classRes_val
}