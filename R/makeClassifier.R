
#' make a classifier, with a randomized class too
#' make a classifier, with a randomized class too
#' @param expTrain training data
#' @param genes vector of genes to use as predictors
#' @param groups named vector of cells to groups or classes
#' @param nRand =50 num of randomized profiles to make
#' @param ntrees =2000 number of trees to build
#' @importFrom randomForest randomForest
#' @return Random Forest Classifier object
#' @export
makeClassifier<-function(expTrain, genes, groups, nRand=50, ntrees=2000){
  randDat<-randomize(expTrain, num=nRand) # randomizes the
  expTrain<-cbind(expTrain, randDat)
  allgenes<-rownames(expTrain)
  missingGenes<-setdiff(unique(genes), allgenes)
  cat("Number of missing genes ", length(missingGenes),"\n")
  ggenes<-intersect(unique(genes), allgenes)
  randomForest(t(expTrain[ggenes,]), as.factor(c(groups, rep("rand", ncol(randDat)))), ntree=ntrees)
}
