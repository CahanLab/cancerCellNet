#' @title
#' Make Classifier tailored for subclassifier 
#' @description
#' Create a random forest classifier with the transformed training data from \code{\link{query_transform}}.
#'
#' @param expTrain transformed training data from \code{\link{query_transform}}
#' @param genes vector of gene pairs from \code{\link{ptGetTop}} used as predictors
#' @param groups named vector of cells to cancer categories
#' @param nRand number of randomized profiles to make
#' @param ntrees number of trees to build
#' @param classMatrix_rand the broad classification for the random category 
#' @param pairTransformedMatrix the pairtransformed matrix of the samples 
#' @importFrom randomForest randomForest
#'
#' @return Random Forest Classifier object for sub-classifier
makeSubClassifier<-function(expTrain, genes, groups, nRand, ntrees=2000, classMatrix_rand, pairTransformedMatrix){
  #randDat<-randomize(expTrain, num=nRand) # OG randomize

  randDat<-ModifiedRandomize(pairTransformedMatrix, num=ncol(classMatrix_rand))
  randDat = rbind(randDat, classMatrix_rand)
  expTrain<-cbind(expTrain, randDat)

  allgenes<-rownames(expTrain)
  missingGenes<-setdiff(unique(genes), allgenes) # I don't think this line does anyhting
  cat("Number of missing genes ", length(missingGenes),"\n")
  ggenes<-intersect(unique(genes), allgenes)

  rf = randomForest::randomForest(t(expTrain[ggenes,]), as.factor(c(groups, rep("rand", ncol(randDat)))), ntree=ntrees)

  returnList = list(tspRF = rf, trainingExp = expTrain[ggenes,], namedVector = as.factor(c(groups, rep("rand", ncol(randDat)))))

  # return 
  returnList
}
