
#' @title
#' Classify the samples
#'
#' @description
#' This functions will classify the query samples using the random forest object created from \code{\link{makeClassifier}}
#'
#' @param rfObj random forest object created from \code{\link{makeClassifier}}
#' @param expQuery transformed query data created from \code{\link{query_transform}}
#' @param numRand number of random profiles wish to create
#' @return a classification matrix of the query data
#' @export
rf_classPredict<-function(rfObj, expQuery, numRand=50) {

  randDat<-randomize(expQuery, num=numRand) # generate the random profiles
  expQuery<-cbind(expQuery, randDat)

  preds<-rownames(rfObj$importance)
  xpreds<-t(predict(rfObj, t(expQuery[preds,]), type='prob'))
  colnames(xpreds)<-colnames(expQuery)
  xpreds
}
