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
#' @param stratify whether to use stratified sampling or not
#' @param samplesize the samplesize for straified sampling
#' @importFrom randomForest randomForest
#'
#' @return Random Forest Classifier object
#' @export
makeClassifier<-function(expTrain, genes, groups, nRand=50, ntrees=2000, stratify=FALSE, sampsize=40){
  randDat = randomize(expTrain, num=nRand)
  #randDat = ModifiedRandomize(expTrain, num=nRand)

  expTrain = cbind(expTrain, randDat)

  allgenes = rownames(expTrain)
  missingGenes = setdiff(unique(genes), allgenes)
  cat("Number of missing genes ", length(missingGenes),"\n")
  ggenes = intersect(unique(genes), allgenes)

  # return random forest object
  if(!stratify){
    return(randomForest::randomForest(t(expTrain[ggenes,]), as.factor(c(groups, rep("rand", ncol(randDat)))), ntree=ntrees))
  }else{
    return(randomForest::randomForest(t(expTrain[ggenes,]), as.factor(c(groups, rep("rand", ncol(randDat)))), ntree=ntrees, strata = as.factor(c(groups, rep("rand", ncol(randDat)))), sampsize=rep(sampsize, length(c(unique(groups), "rand")))))
  }
}
