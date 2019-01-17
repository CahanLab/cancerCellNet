#' @title
#' Sub-class Query Transform
#'
#' @description
#' Transform the sub-class query expression data 
#'
#' @param expDat query expression data 
#' @param cgenes classification genes 
#' @param xpairs classification gene pairs 
#' @param classMatrix classification matrix of the query expression data from broad classifier 
#' @return transformed sub-class query matrix 
#'
#' @export
subClassQuery_transform<-function(expDat, cgenes, xpairs, classMatrix) {
   expValTrans = query_transform(expDat[cgenes,], xpairs)

   extraFeatures = classMatrix[, colnames(expValTrans)]
   if (all(colnames(expValTrans) == colnames(extraFeatures))) {
      print("All Good")
   }

   returnMatrix = rbind(expValTrans, extraFeatures)
   
   #return 
   returnMatrix
}