#' @title
#' Make Classifier
#'
#' @description
#' Create a random forest classifier with the transformed training data from \code{\link{query_transform}}.
#'
#' @param expTrain transformed training data from \code{\link{query_transform}}
#' @param genes vector of gene pairs from \code{\link{ptGetTop}} used as predictors
#' @param groups named vector of cells to cancer categories
#' @param nRand number of randomized profiles to make
#' @param ntrees number of trees to build
#' @importFrom randomForest randomForest
#'
#' @return Random Forest Classifier object
#' @export
makeClassifier<-function(expTrain, genes, groups, nRand=50, ntrees=2000){
  randDat<-randomize(expTrain, num=nRand) # randomizes the
  expTrain<-cbind(expTrain, randDat)

  allgenes<-rownames(expTrain)
  missingGenes<-setdiff(unique(genes), allgenes) # I don't think this line does anyhting
  cat("Number of missing genes ", length(missingGenes),"\n")
  ggenes<-intersect(unique(genes), allgenes)

  randomForest::randomForest(t(expTrain[ggenes,]), as.factor(c(groups, rep("rand", ncol(randDat)))), ntree=ntrees)
}
