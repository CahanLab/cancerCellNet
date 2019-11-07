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
ccn_netScores<-function (expDat, genes, tVals, ctt, classList=NULL, classWeight=FALSE, exprWeight=TRUE,xmax=1e3){
  cat(ctt,"\n")
  aMat<-matrix(0, nrow=length(genes), ncol=ncol(expDat));
  rownames(aMat)<-genes;

  weights<-rep(1, length(genes));
  names(weights)<-genes;

  #otherCTs<-setdiff(names(tVals), ct)

  cat(dim(aMat),"\n")
  if(exprWeight){
    meanVect<-unlist(tVals[[ctt]][['mean']][genes]);
    weights<-(2**meanVect)/sum(2**meanVect);
  }

  if(classWeight){ #TODO modify this to fit the gene pairs
    #classImp<-classList[[ctt]]$importance[genes,1];

    classImp = weights
    for(gene in names(classList[[ctt]])) {
      if(gene %in% names(classImp)) {
        classImp[gene] = classList[[ctt]][gene] + 1
      }
    }

    ### 04-19-17
    ###classImp<-classImp/sum(classImp)
    weights<-weights*classImp;
  }

  for(gene in genes){
    ### cat("***",gene,"\n")
    ###zzs<-as.matrix(cn_rawScore(expDat[gene,], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]])[1,])

    zzs<-ccn_rawScore(expDat[gene,], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]], xmax=xmax)
    aMat[gene,]<-zzs;

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
#'
#' @return transformed (but not normalized) GRN score
#'
ccn_rawScore<-function(vect, mmean, ssd, xmax=1e3){
  zcs<-zscore(vect, mmean, ssd);
  ### xmax<-1000; # arbitrary, and corrected for later, but want some high enough that it should not be exceeded

  xmax-abs(zcs);
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
#'
ccn_score<-function(expDat, subList, tVals, classList=NULL, minVals=NULL, classWeight=FALSE, exprWeight=TRUE, xmax=1e3){
  #nSubnets<-sum(sapply(subList, length));
  nSubnets<-length(subList);
  ans<-matrix(0, nrow=nSubnets, ncol=ncol(expDat));
  ctts<-names(subList);
  rnames<-vector();
  rIndex<-1;
  for(ctt in ctts){
    cat(ctt,"\n");
    genes<-subList[[ctt]];
    # 06-06-16 -- added to allow for use of GRNs defined elsewhere
    genes<-intersect(genes, rownames(expDat));
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
#'
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
#'
#' @return list of trainingScores, normVals, raw_scores, minVals, tVals=tVals
#' @export
ccn_trainNorm<-function (expTrain, stTrain, subNets, classList = NULL,  dLevel = "description1", tVals=NULL, classWeight=FALSE, exprWeight=TRUE, sidCol='sample_id', xmax=1e3, predSD=FALSE){

  if(is.null(tVals)){
    tVals<-ccn_make_tVals(expTrain, stTrain, dLevel, predictSD=predSD)
  }

  ctts<-as.vector(unique(stTrain[,dLevel]));
  scoreList<-list();
  normList<-list(); # a list of ctt->subnet->mean value
  minVect<-vector(); # a list of ctt->subnet->min value, used to shift raw grn est scores

  cat("calculating GRN scores on training data ...\n");
  tmpScores<-ccn_score(expTrain, subNets, tVals, classList, minVals=NULL, classWeight=classWeight, exprWeight=exprWeight, xmax=xmax)


  minVect<-apply(tmpScores, 1, min);
  names(minVect)<-rownames(tmpScores);

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
