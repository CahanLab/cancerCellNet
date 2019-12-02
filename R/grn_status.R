#' GRN status
#'
#' Calculates the status of all GRNs in query samples as compared to training data for
#' @param expDat query expression matrix
#' @param subList of ct => genes
#' @param tVals tvals
#' @param classList classList
#' @param minVals minVals
#' @param classWeight class weight
#' @param exprWeight  expression weight
#' @return grn scores (not normalized)
ccn_netScores<-function (expDat, genes, tVals, ctt, classList=NULL, classWeight=TRUE, classWeightVal = 3, exprWeight=TRUE, exprWeightVal = 3, xmax=1e3){
  cat(ctt,"\n")
  aMat<-matrix(0, nrow=length(genes), ncol=ncol(expDat));
  rownames(aMat)<-names(genes);

  weights<-rep(1, length(genes));
  names(weights)<-names(genes);

  #otherCTs<-setdiff(names(tVals), ct)

  cat(dim(aMat),"\n")
  if(exprWeight){
    meanVect<-unlist(tVals[[ctt]][['mean']][names(genes)]);
    weights<-(exprWeightVal*meanVect)/sum(exprWeightVal*meanVect); # also arbritary value on the weight you are putting on the expression
  }

  if(classWeight){ #TODO modify this to fit the gene pairs
    #classImp<-classList[[ctt]]$importance[genes,1];

    classImp = weights
    for(gene in names(classList[[ctt]])) {
      if(gene %in% names(classImp)) {
        classImp[gene] = classList[[ctt]][gene] + classWeightVal # arbritary value
      }
    }

    ### 04-19-17
    ###classImp<-classImp/sum(classImp)
    weights<-weights*classImp;
  }

  for(gene in names(genes)){

    ### cat("***",gene,"\n")
    ###zzs<-as.matrix(cn_rawScore(expDat[gene,], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]])[1,])

    zzs<-ccn_rawScore(expDat[gene,], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]], xmax=xmax, reg_type = genes[gene])

    aMat[gene,] = zzs;
  }

  print(dim(aMat))
  xscores<-apply(aMat, 2, weighted.mean, w=weights);
  xscores;
}

#' Estimate gene expression dist in CTs
#'
#' Calculate mean and SD
#' @param expDat training data
#' @param sampTab, ### training sample table
#' @param dLevel="description1", ### column to define CTs
#' @param predictSD=FALSE ### whether to predict SD based on expression level
#'
#' @return tVals list of ct->mean->named vector of average gene expression, ->sd->named vector of gene standard deviation
ccn_make_tVals<-function (expDat, sampTab, dLevel="description1", predictSD=FALSE){

  if(predictSD){
    ans<-ccn_make_tVals_predict(expDat, sampTab, dLevel);
  }
  else{
    # Note: returns a list of dName->gene->mean, sd, where 'dName' is a ctt or lineage
    # make sure everything is lined up
    expDat<-expDat[,rownames(sampTab)];
    tVals<-list();
    dNames<-unique(as.vector(sampTab[,dLevel]));
    allGenes<-rownames(expDat);
    for(dName in dNames){
      #cat(dName,"\n");
      xx<-which(sampTab[,dLevel]==dName);
      sids<-rownames(sampTab[xx,]);
      xDat<-expDat[,sids];
      means<-apply(xDat, 1, mean);
      sds<-apply(xDat, 1, sd);
      tVals[[dName]][['mean']]<-as.list(means);
      tVals[[dName]][['sd']]<-as.list(sds);
    }
    ans<-tVals;
  }
  ans;
}

#' ccn_make_tVals_predict
#'
#' predicts SD based on mean expression using linear regression
#' @param expDat training data
#' @param sampTab training sample table
#' @param dLevel="description1" column to define CT
#'
#' @return tVals list of ct->mean->named vector of average gene expression, ->sd->named vector of gene standard deviation
ccn_make_tVals_predict<-function(expDat, sampTab, dLevel="description1"){
  # Note: returns a list of dName->gene->mean, sd, where 'dName' is a ctt or lineage
  # make sure everything is lined up
  expDat<-expDat[,rownames(sampTab)];
  tVals<-list();
  dNames<-unique(as.vector(sampTab[,dLevel]));
  allGenes<-rownames(expDat);

  # make a model to predict SD given average expression level across all samples
  sdT<-apply(expDat, 1, sd);
  mT<-apply(expDat, 1, mean);
  myModel<-lm(sdT~mT);
  for(dName in dNames){
    xx<-which(sampTab[,dLevel]==dName);
    sids<-rownames(sampTab[xx,]);
    xDat<-expDat[,sids];
    means<-apply(xDat, 1, mean);
    sds<-predict(myModel, data.frame(mT=means));
    tVals[[dName]][['mean']]<-as.list(means);
    tVals[[dName]][['sd']]<-as.list(sds);
  }
  tVals;
}

#' Min Difference
#'
#' computes mean gene A in CT 1 - mean gene A in CT 2, where CT2 has the non CT1 max value. does this for genes
#' @param tVals tVals
#' @param genes vector of gene names
#' @param ct ct to compare to
#' @return vector of differences
minDif<-function(tVals, genes, ct){
  octs<-setdiff(names(tVals), ct)
  qq<-lapply(tVals[octs], "[[", "mean")
  ##tVals[[ct]][["mean"]][[gene]]#-max(unlist(lapply(qq, "[[", gene)))
  tmpMat<-matrix(unlist(lapply(qq, "[", genes)), nrow=length(genes))
  rownames(tmpMat)<-genes
  maxes<-apply(tmpMat, 1, max)
  unlist(tVals[[ct]][["mean"]][genes])-maxes
}

#' computes the raw z score for a gene as xmax-abs(zscore).
#'
#' better values are higher.
#' @param vect a vector of gene expression values for multiple samples
#' @param mmean mean value in training data
#' @param ssd standard deviation in training data
#' @param reg_type the type of regulation (up or down) for the specific subnetwork
#' @return transformed (but not normalized) GRN score
#'
ccn_rawScore<-function(vect, mmean, ssd, xmax=1e3, reg_type){
  zcs<-zscore(vect, mmean, ssd);

  if(as.numeric(reg_type) == 1) { # if the
    zcs[zcs > 0] = 0
    return(xmax - abs(zcs * 3))
  }
  else {
    zcs[zcs < 0] = 0

    return(xmax - abs(zcs * 3))
  }
  ### xmax<-1000; # arbitrary, and corrected for later, but want some high enough that it should not be exceeded

  #xmax-abs(zcs); # this was the original scoring system.
  #zcs
}

#' GRN status
#'
#' Calculates the GRN status in query samples as compared to training data
#' @param expDat query expression matrix
#' @param subList of ct => genes
#' @param tVals list of ctt->list(means->genes, sds->genes)
#' @param classList class list
#' @param minVals  min vals
#' @param classWeight classweght
#' @param exprWeight  expression weight
#' @return GRN scores
#' @export
ccn_score<-function(expDat, subList, tVals, classList=NULL, minVals=NULL, classWeight=FALSE, exprWeight=TRUE, xmax=1e3){
  #nSubnets<-sum(sapply(subList, length));
  if(class(expDat) != "matrix") {
    expDat = as.matrix(expDat)
  }
  nSubnets<-length(subList);
  ans<-matrix(0, nrow=nSubnets, ncol=ncol(expDat));
  ctts<-names(subList);
  rnames<-vector();
  rIndex<-1;
  for(ctt in ctts){
    cat(ctt,"\n");
    genes<-subList[[ctt]];
    # 06-06-16 -- added to allow for use of GRNs defined elsewhere
    genes<-genes[intersect(names(genes), rownames(expDat))]; # only select the genes that are
    #    snNames<-names(subnets);
    #    rnames<-append(rnames, snNames);
    #    for(sName in snNames){
    ans[rIndex,]<-ccn_netScores(expDat, genes, tVals=tVals, ctt=ctt,classList=classList, classWeight=classWeight,exprWeight=exprWeight, xmax=xmax);
    rnames<-append(rnames, ctt);
    rIndex<-rIndex+1;
    #   }
  }
  rownames(ans)<-rnames;
  colnames(ans)<-colnames(expDat);
  if(!is.null(minVals)){
    minVals<-minVals[rownames(ans)];
    ans<-ans-minVals;
  }
  ans;
}

#' Normalize grn status as compared to training data
#'
#' Divide the query scores by the mean values in the training data.
#' @param ctrlScores a list of subnet->mean value, all subnets
#' @param queryScores a matrix, rownames = subnet names, all subnets
#' @param subNets a vector of subnets names to use
#'
#' @return normalized grn status matrix
#' @export
ccn_normalizeScores<-function(ctrlScores, queryScores, subNets){

  ans<-matrix(0, nrow=length(subNets), ncol=ncol(queryScores));
  rownames(ans)<-subNets
  #subNets<-rownames(queryScores);
  for(subNet in subNets){
    ### cat(subNet,"\n")
    ans[subNet,]<- queryScores[subNet,] / ctrlScores[[subNet]];
  }
  colnames(ans)<-colnames(queryScores);
  ans;
}


#' Figure out normalization factors for GRNs, and norm training data
#'
#' Exactly that.
#' @param expTrain expression matrix
#' @param stTrain  sample table
#' @param subNets named list of genes, one list per CTT, tct=>gene vector
#' @param classList list of classifiers
#' @param dLevel column name to group on
#' @param tVals seful when debugging
#' @param classWeight weight GRN status by importance of gene to classifier
#' @param exprWeight weight GRN status by expression level of gene?
#' @param sidCol sample id colname
#' @param xmax the maximum raw score that a sample could receive per gene
#' @param meanNorm normalize raw scores based on the lowest mean in a category
#'
#' @return list of trainingScores, normVals, raw_scores, minVals, tVals=tVals
#' @export
ccn_trainNorm<-function (expTrain, stTrain, subNets, classList = NULL,  dLevel = "description1", tVals=NULL, classWeight=TRUE, exprWeight=FALSE, sidCol='sample_id', xmax=1e3, meanNorm = FALSE){

  if(is.null(tVals)){
    tVals<-ccn_make_tVals(expTrain, stTrain, dLevel)
  }

  ctts<-as.vector(unique(stTrain[,dLevel]));
  scoreList<-list();
  normList<-list(); # a list of ctt->subnet->mean value
  minVect<-vector(); # a list of ctt->subnet->min value, used to shift raw grn est scores

  cat("calculating GRN scores on training data ...\n");
  tmpScores<-ccn_score(expTrain, subNets, tVals, classList, minVals=NULL, classWeight=classWeight, exprWeight=exprWeight, xmax=xmax)


  if(meanNorm == TRUE) {
    train_meanScores = meanTraining(tmpScores, stTrain, dLevel, sidCol)
    minVect<-apply(train_meanScores, 1, min);
    names(minVect)<-rownames(train_meanScores);

  } else {
    minVect<-apply(tmpScores, 1, min);
    names(minVect)<-rownames(tmpScores);

  }


  # shift the raw scores so that min=0;
  tmpScores<-tmpScores - minVect;
  cat("norm factors\n");
  for(ctt in ctts){
    # determine nomalization factors
    ##snets<-names(subNets[[ctt]]);
    snets<-ctt;

    scoreDF<-ccn_extract_SN_DF(tmpScores, stTrain, dLevel, snets, sidCol=sidCol);
    scoreDF<-ccn_reduceMatLarge(scoreDF, "score", "description", "subNet");
    xdf<-scoreDF[which(scoreDF$grp_name==ctt),];
    tmpSNS<-as.list(xdf$mean);
    names(tmpSNS)<-xdf$subNet;
    normList[names(tmpSNS)]<-tmpSNS;
  }

  # normalize training scores
  nScores<-ccn_normalizeScores(normList, tmpScores, rownames(tmpScores));

  scoreDF<-ccn_extract_SN_DF(nScores, stTrain, dLevel, sidCol=sidCol);

  scoreDF<-ccn_reduceMatLarge(scoreDF, "score", "description", "subNet");

  list(trainingScores=scoreDF,
       normVals=normList,
       raw_scores=tmpScores,
       minVals=minVect,
       tVals=tVals);
}

#' @title calculate the mean GRN scores across categories
#' @description calculate the mean GRN scores across categories
#' @param grnScores the GRN scores calculated through ccn_score
#' @param stTrain the sample table used for training
#' @param dLevel the name of the column with all the categories
#' @return an averaged GRN score matrix
meanTraining <- function(grnScores, stTrain, dLevel, sidCol) {
  rownames(stTrain) = as.vector(stTrain[, sidCol])
  # return matrix
  meanMatrix = matrix(0, nrow = nrow(grnScores), ncol = length(unique(stTrain[, dLevel])))

  rownames(meanMatrix) = rownames(grnScores)
  colnames(meanMatrix) = unique(stTrain[, dLevel])

  for(cancerName in unique(stTrain[, dLevel])) {
    tempStTrain = stTrain[stTrain[, dLevel] == cancerName, ]
    tempGRNscores = grnScores[, rownames(tempStTrain)]

    meanMatrix[, cancerName] = apply(tempGRNscores, 1, mean)
  }

  return(meanMatrix)
}
#' Z scores
#' Figure out Z score given mean and standard deviation
#' @param x query score
#' @param meanVal mean of the distribution
zscore<-function(x,meanVal,sdVal){

  (x-meanVal)/sdVal;
}


#' returns a DF of: sample_id, description, ctt, subnet_name, score
#'
#' returns a DF of: sample_id, description, ctt, subnet_name, score
#' @param scores a matrix of subNet scores
#' @param sampTab sample table
#' @param dLevel column name of sampTab to group values by
#' @param rnames rownames to extract
#' @param sidCol sample identifier column name
#'
#' @return returns a DF of: sample_id, description, ctt, subnet_name, score
ccn_extract_SN_DF<-function(scores, sampTab, dLevel, rnames=NULL, sidCol="sample_id"){

  if(is.null(rnames)){
    rnames<-rownames(scores);
    #cat("GOT NULL\n");
  }

  tss<-scores[rnames,];
  if(length(rnames)==1){
    tss<-t(as.matrix(scores[rnames,]));
    rownames(tss)<-rnames;
    #  cat(dim(tss),"\n")
  }

  nSamples<-ncol(tss);
  stTmp<-sampTab[colnames(tss),]; ####
  snNames<-rownames(tss);
  num_subnets<-length(snNames);
  snNames<-unlist(lapply(snNames, rep, times=nSamples));
  sample_ids<-rep(as.vector(stTmp[,sidCol]), num_subnets);
  descriptions<-rep(as.vector(stTmp[,dLevel]), num_subnets);
  # myCtts<-rep(ctt, length(snNames));
  scores<-as.vector(t(tss));
  data.frame(sample_id=sample_ids,
             description=descriptions,
             #         ctt=myCtts,
             subNet = snNames,
             score=scores);
  ### data.frame
}

#' reduce large matrix
#' reduce large matrix from ccn_extract_SN_DF
#'
#' @param datFrame the result from ccn_extract_SN_DF
#' @param valCol the column name of the score
#' @param cName the column name of the different groups
#' @param iterOver the column name of the subnetworks
#'
#' @return reduced large matrix
ccn_reduceMatLarge<-function (datFrame, valCol="score", cName="description", iterOver="subNet"){

  iterOvers<-unique(as.vector(datFrame[,iterOver]));
  ans<-data.frame();
  for(io in iterOvers){
    #  cat(io,"\n");
    xi<-which(datFrame[,iterOver]==io);
    dfTmp<-datFrame[xi,];
    x<- utils_reduceMat(dfTmp,valCol=valCol,cName=cName);
    x<-cbind(x, subNet=rep(io, nrow(x)));
    ans<-rbind(ans, x);
  }
  ans;
  ### ans
}

#' @title  Reduce data matrix
#' @description reduce the data.matrix values by averaging and getting st dvs
#'
#' @param datFrame the dataframe
#' @param valCol the column with scores
#' @param cName column with groups
#'
#' @return df of grp_name, mean, sd
utils_reduceMat<-function(datFrame, valCol, cName='ann'){

  mids<-unique(as.vector(datFrame[,cName]));
  means<-rep(0, length(mids));
  sds<-rep(1, length(mids));
  indexI<-rep(0, length(mids)); # index to extract other columns from original data.frame

  for(i in seq(length(mids))){
    mid<-mids[i]
    xi<-which(datFrame[,cName]==mid);
    tmpDat<-datFrame[xi,];
    means[i]<-mean(tmpDat[,valCol]);
    sds[i]<-sd(tmpDat[,valCol]);
  }
  data.frame(grp_name=mids, mean=means, stdev=sds);
}

#' @title Get Query GRN status
#' @description Get the GRN status of query samples
#'
#' @param expQuery logRanked query expression matrix
#' @param expTrain logRanked training expression matrix
#' @param stTrain sample table of training expression matrix
#' @param dLevel the name of the column with cancer types
#' @param sidCol the name of the column with sample IDs
#' @param grn_return the grn list that is returned from \code{\link{ccn_makeGRN}}
#' @param trainNorm normalization statistics from \code{\link{ccn_trainNorm}}. If you are using pre-calculated normalization statistics please make sure all the parameters are the same for applying to query and calculating normalization
#' @param classifier_return the classifier_return list that is returned from \code{\link{broadClass_train}}
#' @param classWeight boolean indicating whether to take the importance of the classification into status calculation
#' @param exprWeight boolean indicating whether to take the weight of gene expression into status calculation
#' @param prune boolean indicating whether to select exclusive genes for processing classification gene importance
#' @param predSD a parameter for calculating normalization statistics from training data
#' @return a matrix indicating the GRN status
#' @export
ccn_queryGRNstatus <- function(expQuery, expTrain, stTrain, dLevel, sidCol, grn_return, trainNorm = NULL, classifier_return,  classWeight = TRUE, exprWeight = FALSE, prune = TRUE, xmax = 1e3, predSD=FALSE) {
  cnProc = classifier_return$cnProc

  trainNorm_prior = trainNorm
  geneImportance = processImportance(classifier = cnProc$classifier, xpairs = classifier_return$xpairs_list, prune = prune)

  if(is.null(trainNorm) == TRUE) {
    trainNorm = ccn_trainNorm(expTrain, stTrain, subNets=grn_return$ctGRNs$geneLists, classList = geneImportance, dLevel = dLevel, sidCol = sidCol, classWeight = classWeight, exprWeight = exprWeight)
    cat("Finished finding normalization statistics", "\n")
  }

  status_score = ccn_score(expDat = expQuery,
                           subList = grn_return$ctGRNs$geneLists, tVals = trainNorm$tVals,
                           classList = geneImportance, minVals = trainNorm$minVals,
                           classWeight = classWeight, exprWeight = exprWeight,
                           xmax = xmax)
  normScoresQuery = ccn_normalizeScores(trainNorm$normVals, status_score, rownames(status_score))

  if(is.null(trainNorm_prior) == TRUE) {
    return(list("trainNorm" = trainNorm, "query_GRNstatus" = normScoresQuery))
  }
  else {
    return(normScoresQuery)
  }
}





