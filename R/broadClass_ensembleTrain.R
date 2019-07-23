library(cancerCellNet)

stTrain = utils_loadObject("~/Dropbox (CahanLab)/Data_Repo/GDCTrainingData/Named_stGDC_20181218.rda")
expTrain = utils_loadObject("~/Dropbox (CahanLab)/Data_Repo/GDCTrainingData/Named_expGDC_20181218.rda")

iGenes = utils_loadObject("~/Dropbox (CahanLab)/Pure_Testing/Training_WO_GEMMSARC/iGenes_noSarc.rda")

expTrain = expTrain[iGenes, ]
splitCommon(stTrain, ncells = 4000, dLevel = "project_id")

myTest = broadClass_ensembleTrain(stTrain = stTrain, expTrain = expTrain, colName_cat = "project_id", numClassifier = 3, ncells = 60, nTopGenes = 5, nTopGenePairs = 20)

broadClass_ensembleTrain <- function(stTrain, expTrain, colName_cat, numClassifier = 5, ncells = 60, colName_samp="row.names", nTopGenes = 20, nTopGenePairs = 50, nRand = 40, nTrees = 1000, weightedDown_total = 5e5, weightedDown_dThresh = 0.25, transprop_xFact = 1e5) {

  if (class(stTrain) != "data.frame") {
    stTrain<-as.data.frame(stTrain)
  }

  if (colName_samp != "row.names") {
    rownames(stTrain)<-stTrain[, colName_samp]
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
