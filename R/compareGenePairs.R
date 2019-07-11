#' @title
#' Compile genePairs comparison matrix
#' @description
#' compile genePairs comparison matrix for users to elucidate the biological significance behind classification results
#'
#' @param query_exp the expression matrix of query samples
#' @param training_exp the expression matrix of training samples
#' @param training_st the sample table of training data
#' @param classCol the column name of the sample table that indicates the classes
#' @param sampleCol the column name of the sample table that indicates the sample names. NULL if sample names are indcated in rownames of the sample table
#' @param cnProc the cnProc from the training
#' @param numPairs the number of genes to extract for comparison from most important to least important in the classifier
#'
#' @return gene pair comparison matrix
#'
#' @export
compareGenePairs<-function(query_exp, training_exp, training_st, classCol, sampleCol = NULL, cnProc, numPairs = 20) {

  RF_classifier = cnProc$classifier

  importantPairs = RF_classifier$importance
  importantPairs = sort(importantPairs[, 1], decreasing = TRUE)
  importantPairs = importantPairs[grep("_", names(importantPairs))]

  if (numPairs < length(importantPairs)) {
    userPairs = importantPairs
  }
  else {
    userPairs = importantPairs[1:numPairs]
  }

  query_pairs = query_transform(expDat = query_exp, genePairs = names(userPairs))

  training_pairs = query_transform(expDat = training_exp, genePairs = names(userPairs))

  avg_training_pairs = avgGeneCat(expDat = training_pairs, sampTab = training_st, dLevel = classCol, sampID = sampleCol)

  PairCompareMatrix = makeGeneCompareTab(queryExpTab = query_pairs,
                                         avgGeneTab = avg_training_pairs,
                                         geneSamples = names(userPairs))

  #return
  PairCompareMatrix
}
