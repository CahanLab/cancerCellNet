#' @title
#' Broad Class Ensemble Training (not commonly used)
#' @description
#' Tranining an ensemble of broad class classifiers
#' @param stTrain a dataframe that matches the samples with category
#' @param expTrain the expression matrix
#' @param colName_cat the name of the column that contains categories
#' @param colName_samp the name of the column that contains sample names
#' @param numClassifier the number of classifiers in the ensemble
#' @param ncells the number of samples per category for training
#' @param nTopGenes the number of classification genes per category
#' @param nTopGenePairs the number of top gene pairs per category
#' @param nRand number of random profiles generate for training
#' @param nTrees number of trees for random forest classifier
#' @param weightDown_total numeric post transformation sum of read counts for weighted_down function
#' @param weightedDown_dThresh the threshold at which anything lower than that is 0 for weighted_down function
#' @param transprop_xFact scaling factor for transprop
#'
#' @return a list of broadclass_returns with classifiers
#' @export
broadClass_ensembleTrain <- function(stTrain, expTrain, colName_cat, colName_samp="row.names", numClassifier = 5, ncells = 60, nTopGenes = 20, nTopGenePairs = 50, nRand = 40, nTrees = 1000, weightedDown_total = 5e5, weightedDown_dThresh = 0.25, transprop_xFact = 1e5) {

  if (class(stTrain) != "data.frame") {
    stTrain = as.data.frame(stTrain)
  }

  if (colName_samp != "row.names") {
    rownames(stTrain) = stTrain[, colName_samp]
  }

  returnList = list()
  for(i in 1:numClassifier) {
    stList = splitCommon(stTrain, ncells = ncells, dLevel = colName_cat)

    this_stTrain = stList[[1]]
    this_expTrain = expTrain[, rownames(this_stTrain)]

    tempBroadReturn = broadClass_train(stTrain = this_stTrain,
                                       expTrain = this_expTrain,
                                       colName_cat = colName_cat,
                                       colName_samp = colName_samp,
                                       nTopGenes = nTopGenes,
                                       nTopGenePairs = nTopGenePairs,
                                       nRand = nRand,
                                       nTrees = nTrees,
                                       weightedDown_total = weightedDown_total,
                                       weightedDown_dThresh = weightedDown_dThresh,
                                       transprop_xFact = transprop_xFact)

    returnList[[paste0("broadReturn_", i)]] = tempBroadReturn

  }

  #return
  returnList
}
