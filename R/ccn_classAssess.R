#' Assess classifiers based on validation data
#'
#' Assess classifiers based on validation data
#' @param ct_scores matrix of classification scores, rows = classifiers, columns = samples, colnames=sampleids
#' @param stVal sample table
#' @param classLevels column name of stVal to use as ground truth to assess classifiers
#' @param resolution increment at which to evalutate classification
#' @param dLevelSID column to indicate sample id
#'
#' @return list of data frames with threshold, sens, precision
#' @export
ccn_classAssess<-function# make PR data frames for each classifier
(ct_scores,# matrix of classification scores, rows = classifiers, columns = samples, colnames=sampleids
 stVal, # sampletable
 classLevels="description2",
 dLevelSID="sample_id",
 resolution=0.005 # increment at which to evalutate classification
){
  allPRs<-list()
  evalAll<-matrix(0, nrow=nrow(ct_scores),ncol=2)
  classifications<-rownames(ct_scores)
  rownames(stVal)<-as.vector(stVal[,dLevelSID])
  i<-1;
  for(xname in classifications){
    classification<-classifications[i]
    tmpPR <- cn_eval(ct_scores[xname,],
                     stVal,
                     classLevels,
                     xname,threshs=seq(0,1, by=resolution), dLevelSID=dLevelSID)
    allPRs[[xname]]<-as.data.frame(tmpPR)
    i<-1+i
  }
  allPRs
}

#' run cn_clPerf across thresholds
#'
#' run cn_clPerf across thresholds
#' @param vect named vector
#' @param dLevel description level
#' @param classification classification matrix
#' @param threshs seq of pval cutoffs
#' @param dLevelSID column to indicate sample id
#'
#' @return return a data frame of the number of TP, FN, FP, and TN, and pval cutoff
cn_eval<-function# return a data frame of the number of TP, FN, FP, and TN, and pval cutoff
(vect, # named vector
 sampTab,
 dLevel, # description level)
 classification,
 threshs=seq(0,1,by=0.05), # pval cutoffs
 dLevelSID="sample_id"
){
  ans<-matrix(0,nrow=length(threshs), ncol=9)
  for(i in seq(length(threshs))){
    thresh<-threshs[i]
    ans[i,1:4]<-cn_clPerf(vect, sampTab, dLevel, classification, thresh, dLevelSID=dLevelSID)
  }
  ans[,5]<-threshs
  colnames(ans)<-c("TP", "FN", "FP", "TN", "thresh","FPR", "TPR", "Sens", "Prec")
  TPR<-ans[,'TP']/(ans[,'TP']+ans[,'FN'])
  FPR<-ans[,'FP']/(ans[,'TN']+ans[,'FP'])
  sens<-ans[,"TP"]/(ans[,"TP"]+ans[,"FN"])
  prec<-ans[,"TP"]/(ans[,"TP"]+ans[,"FP"])
  ans[,'TPR']<-TPR
  ans[,'FPR']<-FPR
  ans[,'Sens']<-sens
  ans[,'Prec']<-prec
  ans;
}


#' determine performance of classification at given threshold
#'
#' determine performance of classification at given threshold
#' @param vect vector of values
#' @param sampTab sample table
#' @param dLevel colname
#' @param classification actual classification
#' @param thresh threshold above which to make a call
#' @param dLevelSID column to indicate sample id
#'
#' @return vector of TP FN FP TN
cn_clPerf<-function # assumes rownames(sampTab) == sampTab identifier used as colname for vect
(vect,
 sampTab,
 dLevel,
 classification, # actual classification
 thresh,
 dLevelSID="sample_id"){
  TP<-0;
  FN<-0;
  FP<-0;
  TN<-0;
  sampIDs<-names(vect);
  classes<-as.vector(sampTab[sampIDs,dLevel]);

  ###actualPos<-as.vector(sampTab[sampTab[,dLevel]==classification,]$sample_id);#which(classes==classification));
  actualPos<-as.vector(sampTab[sampTab[,dLevel]==classification,dLevelSID])
  actualNeg<-setdiff(sampIDs, actualPos);

  calledPos<-names(which(vect>thresh));
  calledNeg<-names(which(vect<=thresh));

  TP <- length(intersect(actualPos, calledPos));
  FP <- length(intersect(actualNeg, calledPos));
  FN <- length(actualPos)-TP;
  TN <- length(actualNeg)-FP;
  c(TP, FN, FP, TN);
}
