#' Network Influence Score for all GRNs
#'
#' Runs cn_nis on all GRNs
#' @param cnRes object result of running cn_apply
#' @param cnProc object result of running cn_make_processor
#' @param ctt string indicating the CT to compare against
#' @param relaWeight whether to weight by overall expression such that TFs with higher expression in ctt are more important (1=do the weighting)
#'
#' @return list of numeric matrix of TF scores
#'
#' @export
cn_nis_all<-function(cnRes, cnProc, ctt, relaWeight=1){
  snNames<-names(cnProc$ctGRNs$ctGRNs$graphLists);
  ans<-list()
  for(snName in snNames){
    cat("scoring ", snName,"\n")
    x<-cn_nis(cnRes, cnProc, snName, ctt,relaWeight);
    ans[[snName]]<-x;
  }
  ans;
}

#' network influence score
#'
#' Computes network influence score (NIS). See paper for details.
#' @param cnRes object result of running cn_apply
#' @param cnProc object result of running cn_make_processor
#' @param subnet name of subnet to evaluage
#' @param ctt string indicating the CT to compare against
#'
#' @return numeric matrix where rows are TFs in CT GRN and columns are query samples
#'
#' @export
cn_nis<-function(cnRes, cnProc, subnet, ctt, relaWeight=1){

  tfTargList<-cnProc[['ctGRNs']][['ctGRNs']][['tfTargets']];
  # return a DF of : tfs, nTargets, targetScore, tfScore, totalScore
  nTargets<-vector();
  targetScore<-vector();
  tfScore<-vector();
  totalScore<-vector();
  tfWeights<-vector();

  tfs<-names(tfTargList[[subnet]]);
  netGenes<-cnProc[['grnList']][[subnet]];
  netGenes<-intersect(netGenes, rownames(cnProc[['expTrain']])) # this should intersect with expQuery

  expDat<-cnRes[['expQuery']];
  stQuery<-cnRes[['stQuery']];
  sids<-as.vector(stQuery$sample_id);

  ans<-matrix(0, nrow=length(tfs), ncol=nrow(stQuery));
  rownames(ans)<-tfs;
  colnames(ans)<-sids;

  tVals<-cnProc[['tVals']];

  # compute a matrix of zscores.
  zzzMat<-matrix(0, nrow=length(netGenes), ncol=nrow(stQuery));

  for(i in seq(length(sids))){
    sid<-sids[i];
    #cat("computing zscores ", sid,"\n");
    xvals<-as.vector(expDat[netGenes,sid]);
    names(xvals)<-netGenes;
    zzzMat[,i]<-cn_zscoreVect(netGenes, xvals, tVals, ctt);
  }

  rownames(zzzMat)<-netGenes;
  colnames(zzzMat)<-rownames(stQuery);

  for(sid in sids){
    #cat("tf scoring ", sid,"\n");
    xvals<-as.vector(expDat[,sid]);
    names(xvals)<-rownames(expDat);


    # assign weights

    ### # meanVect<-unlist(tVals[[ctt]][['mean']][netGenes]);
    meanVect<-unlist(tVals[[subnet]][['mean']][netGenes]);
    weights<-(2**meanVect)/sum(2**meanVect); # weight of the training data

    for(i in seq(length(tfs))){

      tf<-tfs[i];

      # zscore of TF relative to target C/T
      ##      tfScore[i]<-zscore(xvals[tf], tVals[[ctt]][['mean']][[tf]], tVals[[ctt]][['sd']][[tf]]);

      tfScore[i]<-zzzMat[tf,sid];

      targs<-tfTargList[[subnet]][[tf]]; #
      targs<-intersect(targs, rownames(cnProc[['expTrain']]));

      # Zscores of TF targets, relative to C/T
      ##      tmp<-cn_zscoreVect(targs, xvals, tVals, ctt );
      tmp<-zzzMat[targs,sid];
      targetScore[i]<-sum(tmp*weights[targs]);

      ## new one:
      totalScore[i]<-targetScore[i] + (length(targs)*tfScore[i]*weights[tf]);

      if(relaWeight!=1){ # don't weight by expression
        meanW<-mean(weights)
        totalScore[i]<- sum(tmp)*meanW + (length(targs)*tfScore[i])*meanW
      }
      nTargets[i]<-length(targs) ;
      tfWeights[i]<-weights[tf];
    }
    xxx<-data.frame(tf=tfs, tfScore=tfScore, targetScore=targetScore, nTargets=nTargets,tfWeight=tfWeights, totalScore=totalScore);
    xxx<-xxx[order(xxx$totalScore),]; # puts the worst ones at top when plotting
    xxx$tf<-factor(xxx$tf, as.vector(unique(xxx$tf)));
    ans[as.vector(xxx$tf),sid]<-as.vector(xxx$totalScore);
  }
  ans;
  # returns network influence score.
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
  ### zscore vector
}
