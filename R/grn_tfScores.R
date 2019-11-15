#' @title calculate TF scores
#' @description calculate the TF that one should target to yield target cell type
#'
#' @param expQuery the logRanked expression matrix for query data
#' @param subnetName the name of the target cell type
#' @param grnAll result of running \code{\link{ccn_makeGRN}}
#' @param trainNorm result of running \code{\link{ccn_trainNorm}}
#' @param geneImportance the importance of gene pairs from classifier. Result of running \code{\link{processImportance}}
#' @param classWeight whether to take the weight of the classifier into calculation
#' @param classWeightVal the value of the classifier weight
#' @param exprWeight whether to take the weight of the expression into calculation
#' @param exprWeightVal the value of the expression weight
#' @param correlationFactor the weight of correlation direction in the network
#'
#' @return matrix of TF scores and query samples
#' @export
ccn_tfScores <- function(expQuery, subnetName, grnAll, trainNorm, geneImportance, classWeight=TRUE, classWeightVal = 3, exprWeight=TRUE, exprWeightVal = 3, correlationFactor = 1) {
  TF_targetList = grnAll[["ctGRNs"]][["tfTargets"]][[subnetName]]

  tfs = names(TF_targetList)
  netGenes = grnAll[["ctGRNs"]][["geneLists"]][[subnetName]]

  grnTable = grnAll$overallGRN$grnTable
  tVals = trainNorm$tVals

  # get the z scores matrix
  z_scoreMat = matrix(0, nrow = length(netGenes), ncol = ncol(expQuery))
  rownames(z_scoreMat) = names(netGenes)
  colnames(z_scoreMat) = colnames(expQuery)

  # get the tf scores matrix
  tf_scoreMat = matrix(0, nrow = length(tfs), ncol = ncol(expQuery))
  rownames(tf_scoreMat) = tfs
  colnames(tf_scoreMat) = colnames(expQuery)

  print(tfs)
  for(querySample in colnames(expQuery)) {
    xvals = as.vector(expQuery[names(netGenes), querySample])
    names(xvals) = names(netGenes)
    z_scoreMat[, querySample] = ccn_zscoreVect(genes = names(netGenes), xvals = xvals, tVals = tVals, ctt = subnetName)

  }

  # assign weights to each gene
  weights = rep(1, length(netGenes))
  names(weights) = names(netGenes)
  if(exprWeight){
    meanVect = unlist(tVals[[subnetName]][['mean']][names(netGenes)]);
    weights = (exprWeightVal*meanVect)/sum(exprWeightVal*meanVect); # also arbritary value on the weight you are putting on the expression
  }

  if(classWeight){ #TODO modify this to fit the gene pairs
    classImp = weights
    for(gene in names(classList[[subnetName]])) {
      if(gene %in% names(classImp)) {
        classImp[gene] = classList[[subnetName]][gene] + classWeightVal # arbritary value
      }
    }
    weights = weights*classImp;
  }


  for(sampleID in colnames(z_scoreMat)) {

    tfScore = vector() # the vector of tf z-scores
    # to calculate TFs
    cat("Calculating TF scores for", sampleID, "\n")

    for(i in seq(1, length(tfs))) {
      tf = tfs[i]
      targs = TF_targetList[[tf]]; #
      targs = intersect(targs, rownames(expQuery))

      temp_tfScore = calc_tfScores(tf, targs, sampleID, z_scoreMat, netGenes, weights, grnTable, correlationFactor)
      tfScore = c(tfScore, temp_tfScore)
    }

    tf_scoreMat = tf_scoreMat[names(tfScore), ]
    tf_scoreMat[, sampleID] = tfScore
  }

  return(tf_scoreMat)
}

#' Calculate the TF
#' Calculate the TF score given all the parameters
#'
#' @param tf TF
#' @param targs a list of targets for TF
#' @param sampleID the sample ID
#' @param z_scoreMat the z value matrix
#' @param netGenes genes selected for building the subnetwork
#' @param weights a vector of weights for each gene
#' @param grnTable the grn table with correlation values and signs
#' @param correlationFactor the weight of correlation between
calc_tfScores <- function(tf, targs, sampleID, z_scoreMat, netGenes, weights, grnTable, correlationFactor = 1) {

  if(z_scoreMat[tf, sampleID] > 0 & netGenes[tf] == 1) { # if the gene suppose to be higher and is above z
    temp_zscore = 1
  } else if(z_scoreMat[tf, sampleID] < 0 & netGenes[tf] == -1) { # if the gene is suppose to be lower and is below z
    temp_zscore = 1
  } else {
    temp_zscore = 1 + abs(z_scoreMat[tf, sampleID])
  }

  part1 = length(targs) * temp_zscore * weights[tf]

  part2 = 0
  for(targGene in targs) {
    temp_zscore = 1
    if(z_scoreMat[targGene, sampleID] > 0 & netGenes[tf] == 1) { # if the gene suppose to be higher and is above z
      temp_zscore = 1
    } else if(z_scoreMat[targGene, sampleID] < 0 & netGenes[tf] == -1) { # if the gene is suppose to be lower and is below z
      temp_zscore = 1
    } else {
      temp_zscore = 1 + abs(z_scoreMat[targGene, sampleID])
    }

    corr_sign = sign(grnTable[(grnTable$TG == targGene & grnTable$TF == tf), "corr"])
    corr_factor = 0
    if(corr_sign == 1 & netGenes[tf] == netGenes[targGene]) { # if both are positive corr and both up or down together
      corr_factor = correlationFactor
    }
    else if(corr_sign == -1 & netGenes[tf] != netGenes[targGene]) { # negative correlation but they go different direction
      corr_factor = correlationFactor
    }
    else { # positive correlation
      corr_factor = -correlationFactor
    }

    temp_part2 = temp_zscore * weights[targGene] * corr_factor

    part2 = part2 + temp_part2
  }

  # if TF is suppose to be down regulated
  if(netGenes[tf] == -1) {
    return(-(part1 + part2))
  }
  else {
    return(part1 + part2)

  }
}

#' @title Compute mean Z scores of given genes
#' @description Compute the mean Z score of given genes in each sample
#'
#' @param genes the genes
#' @param xvals named vector
#' @param tVals mean and SD of genes in training data
#' @param ctt cell type
#'
#' @return a vector of Z scores
cn_zscoreVect<-function (genes, xvals, tVals, ctt){
  ans<-vector();
  for(gene in genes){
    ans<-append(ans, zscore(xvals[gene], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]]));
  }
  ans;
}


