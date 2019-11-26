#' @title calculate TF scores
#' @description calculate the TF that one should target to yield target cell type
#'
#' @param expQuery the logRanked expression matrix for query data
#' @param subnetName the name of the target cell type
#' @param grnAll result of running \code{\link{ccn_makeGRN}}
#' @param trainNorm result of running \code{\link{ccn_trainNorm}}
#' @param classifier_return the classifier_return list that is returned from \code{\link{broadClass_train}}
#' @param classWeight whether to take the weight of the classifier into calculation
#' @param classWeightVal the value of the classifier weight
#' @param exprWeight whether to take the weight of the expression into calculation
#' @param exprWeightVal the value of the expression weight
#' @param correlationFactor the weight of correlation direction in the network
#' @param normTFscore boolean indicate whether to normalize TF scores
#'
#' @return matrix of TF scores and query samples
#' @export
ccn_tfScores <- function(expQuery, subnetName, grnAll, trainNorm, classifier_return, classWeight=TRUE, classWeightVal = 3, exprWeight=FALSE, exprWeightVal = 3, correlationFactor = 1, prune = TRUE, normTFscore = FALSE) {
  cnProc = classifier_return$cnProc

  classList = processImportance(classifier = cnProc$classifier, xpairs = classifier_return$xpairs_list, prune = prune)


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
      if(normTFscore == TRUE) {
        temp_tfScore = temp_tfScore - calc_normTfScore(tf, targs, sampleID, netGenes, weights, grnTable, correlationFactor)

      }
      tfScore = c(tfScore, temp_tfScore)
    }

    tf_scoreMat = tf_scoreMat[names(tfScore), ]
    tf_scoreMat[, sampleID] = tfScore
  }

  return(tf_scoreMat)
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
ccn_zscoreVect<-function (genes, xvals, tVals, ctt){
  ans<-vector();
  for(gene in genes){
    ans<-append(ans, zscore(xvals[gene], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]]));
  }
  ans;
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

  # revamp this shit
  calculation_matrix = matrix(data = 1, nrow = length(targs), ncol = 8)
  rownames(calculation_matrix) = targs
  colnames(calculation_matrix) = c("TF_direction", "Target_direction","z_score", "z_mod", "weight", "corr", "correlation_factor", "total_score")

  calculation_matrix[, "TF_direction"] = netGenes[tf] # assign the TF direction
  calculation_matrix[, "Target_direction"] = netGenes[targs] # assign the target direction

  temp_zscore = z_scoreMat[targs, sampleID]
  calculation_matrix[, "z_score"] = as.vector(temp_zscore)

  # assign modification z score
  modIndex = (calculation_matrix[, "Target_direction"] == 1 & calculation_matrix[, "z_score"] < 0) # if z score is low and its suppose to be upregulated
  calculation_matrix[modIndex, "z_mod"] = abs(calculation_matrix[modIndex, "z_score"]) + 1

  modIndex = (calculation_matrix[, "Target_direction"] == -1 & calculation_matrix[, "z_score"] > 0) # if z score is high and its suppose to be downregulated
  calculation_matrix[modIndex, "z_mod"] = abs(calculation_matrix[modIndex, "z_score"]) + 1

  # assign weight
  calculation_matrix[, "weight"] = weights[targs]

  # assign correlation direction
  temp_grnTable = grnTable[grnTable$TF == tf, ]
  rownames(temp_grnTable) = as.vector(temp_grnTable$TG)
  temp_grnTable = temp_grnTable[targs, ]
  corr_sign = sign(temp_grnTable[, "corr"])
  calculation_matrix[, "corr"] = corr_sign

  # assign the correlation factor
  calculation_matrix[, "correlation_factor"] = -correlationFactor # default to be negative

  # if correlation between target gene and TF is positive
  # and if both target and TF are suppose to be upregulated or downregulated
  # it's do not penalize
  modIndex = (calculation_matrix[, "corr"] == 1 & calculation_matrix[, "Target_direction"] == calculation_matrix[, "TF_direction"]) # if z score is low and its suppose to be upregulated
  calculation_matrix[modIndex, "correlation_factor"] = abs(correlationFactor)

  # if correlation between target gene and TF is negative
  # and if both target and TF are in opposite direction
  # do not penalize
  modIndex = (calculation_matrix[, "corr"] == -1 & calculation_matrix[, "Target_direction"] != calculation_matrix[, "TF_direction"]) # if z score is low and its suppose to be upregulated
  calculation_matrix[modIndex, "correlation_factor"] = abs(correlationFactor)

  # calculate total score for individual gene
  calculation_matrix[, "total_score"] = calculation_matrix[, "correlation_factor"] * calculation_matrix[, "weight"] * calculation_matrix[, "z_mod"]

  part2 = sum(calculation_matrix[, "total_score"])

  TF_score = part1 + part2

  # if TF is suppose to be down regulated
  if(netGenes[tf] == -1) {
    return(-(TF_score))
  }
  else {
    return(TF_score)

  }

}

#' Calculate the TF normalization constant
#' Calculate the TF score given all the parameters
#'
#' @param tf TF
#' @param targs a list of targets for TF
#' @param sampleID the sample ID
#' @param netGenes genes selected for building the subnetwork
#' @param weights a vector of weights for each gene
#' @param grnTable the grn table with correlation values and signs
#' @param correlationFactor the weight of correlation between
calc_normTfScore <- function(tf, targs, sampleID, netGenes, weights, grnTable, correlationFactor = 1) {

  part1 = length(targs) * weights[tf]
  part2 = 0

  # revamp this shit
  calculation_matrix = matrix(data = 1, nrow = length(targs), ncol = 8)
  rownames(calculation_matrix) = targs
  colnames(calculation_matrix) = c("TF_direction", "Target_direction","z_score", "z_mod", "weight", "corr", "correlation_factor", "total_score")

  calculation_matrix[, "TF_direction"] = netGenes[tf] # assign the TF direction
  calculation_matrix[, "Target_direction"] = netGenes[targs] # assign the target direction

  # assign weight
  calculation_matrix[, "weight"] = weights[targs]

  # assign correlation direction
  temp_grnTable = grnTable[grnTable$TF == tf, ]
  rownames(temp_grnTable) = as.vector(temp_grnTable$TG)
  temp_grnTable = temp_grnTable[targs, ]
  corr_sign = sign(temp_grnTable[, "corr"])
  calculation_matrix[, "corr"] = corr_sign

  # assign the correlation factor
  calculation_matrix[, "correlation_factor"] = -correlationFactor # default to be negative

  # if correlation between target gene and TF is positive
  # and if both target and TF are suppose to be upregulated or downregulated
  # it's do not penalize
  modIndex = (calculation_matrix[, "corr"] == 1 & calculation_matrix[, "Target_direction"] == calculation_matrix[, "TF_direction"]) # if z score is low and its suppose to be upregulated
  calculation_matrix[modIndex, "correlation_factor"] = abs(correlationFactor)

  # if correlation between target gene and TF is negative
  # and if both target and TF are in opposite direction
  # do not penalize
  modIndex = (calculation_matrix[, "corr"] == -1 & calculation_matrix[, "Target_direction"] != calculation_matrix[, "TF_direction"]) # if z score is low and its suppose to be upregulated
  calculation_matrix[modIndex, "correlation_factor"] = abs(correlationFactor)

  # calculate total score for individual gene
  calculation_matrix[, "total_score"] = calculation_matrix[, "correlation_factor"] * calculation_matrix[, "weight"] * calculation_matrix[, "z_mod"]

  part2 = sum(calculation_matrix[, "total_score"])

  TF_score = part1 + part2

  # if TF is suppose to be down regulated
  if(netGenes[tf] == -1) {
    return(-(TF_score))
  }
  else {
    return(TF_score)
  }

}

