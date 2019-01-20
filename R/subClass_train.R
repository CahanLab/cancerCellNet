#' @title
#' Sub-Class Training
#' @description
#' Tranining sub-class classifier
#' @param cnProc_broad the broad cnProc from \code{\link{broadClass_train}}
#' @param stTrain a dataframe that matches the samples with broad category and sub-class
#' @param expTrain the expression matrix
#' @param colName_broadCat the name of the column in sample table that contains broad categories
#' @param colName_subClass the name of the column in sample table that contains sub class
#' @param name_broadCat the name of the broad class in which the subclasses are
#' @param colName_samp the name of the column that contains sample names
#' @param nTopGenes the number of classification genes per category
#' @param nTopGenePairs the number of top gene pairs per category
#' @param nRand number of random profiles generate for training
#' @param nTrees number of trees for random forest classifier
#' @param weightDown_total numeric post transformation sum of read counts for weighted_down function
#' @param weightedDown_dThresh the threshold at which anything lower than that is 0 for weighted_down function
#' @param transprop_xFact scaling factor for transprop
#'
#' @return a list containing normalized expression data, classification gene list, cnPRoc
#' @export
subClass_train<-function(cnProc_broad, stTrain, expTrain, colName_broadCat, colName_subClass, name_broadCat, colName_samp, nTopGenes = 20, nTopGenePairs = 50, nRand = 40, nTrees = 1000, weightedDown_total = 5e5, weightedDown_dThresh = 0.25, transprop_xFact = 1e5) {
   if (class(stTrain) != "data.frame") {
      stTrain = as.data.frame(stTrain)
   }
   rownames(stTrain) = stTrain[, colName_samp]
   cat("Sample table has been prepared\n")

   expTnorm = trans_prop(weighted_down(expTrain, weightedDown_total, dThresh = weightedDown_dThresh), transprop_xFact)
   cat("Expression data has been normalized\n")

   stTrain_sub = stTrain[stTrain[, colName_broadCat] == name_broadCat, ]
   expTnorm_sub = expTnorm[, rownames(stTrain_sub)]
   cat("The sub-class expression data has been partitioned\n")

   system.time(cgenes<-findClassyGenes(expDat = expTnorm_sub, sampTab = stTrain_sub, dLevel = colName_subClass, topX = nTopGenes))
   cat("Finished finding classification genes\n")
   cgenesA<-cgenes[['cgenes']]
   grps<-cgenes[['grps']]

   cgenes_list <- cgenes[['labelled_cgenes']]

   cat("There are ", length(cgenesA), " classification genes\n")

   system.time(xpairs<-ptGetTop(expTnorm_sub[cgenesA,], grps, topX=nTopGenePairs, sliceSize=2000))
   cat("Finished finding top gene pairs\n")

   # some of these might include selection cassettes; remove them
   xi<-setdiff(1:length(xpairs), grep("selection", xpairs))
   xpairs<-xpairs[xi]

   system.time(pdTrain<-query_transform(expTrain[cgenesA, ], xpairs))
   cat("Finished pair transforming the data\n")

   classMatrix = broadClass_predict(cnProc = cnProc_broad, expDat = expTrain)
   classMatrix = classMatrix[, -grep("rand", colnames(classMatrix))]

   cat("Start SubClass Query Transform\n")
   expValTrans = subClassQuery_transform(expDat = expTrain, cgenes = cgenesA, xpairs = xpairs, classMatrix = classMatrix)
   cat("Features have been selected\n")

   newFeatures = c(xpairs, rownames(classMatrix))

   newGrps = as.vector(stTrain[, colName_subClass])
   names(newGrps) = rownames(stTrain)

   system.time(tspRF<-makeClassifier(expValTrans[newFeatures,], genes=newFeatures, groups=newGrps, nRand=nRand, ntrees=nTrees))
   cat("Finished making the classifier \n")

   cnProc_subClass = list("cgenes"= cgenesA, "xpairs"=xpairs, "grps"= newGrps, newFeatures = "newFeatures",  "classifier" = tspRF)

   returnList = list("expTnorm" = expTnorm, "sampTab" = stTrain, "cgenes_list" = cgenes_list, "cnProc_subClass" = cnProc_subClass)

   cat("All Done \n")
   # return
   returnList
}