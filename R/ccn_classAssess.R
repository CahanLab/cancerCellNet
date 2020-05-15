#' @title
#' Assessment of Classifier
#'
#' @description
#' Assess classifiers using validation data
#'
#' @param ct_scores matrix of classification scores
#' @param stVal sample table
#' @param classLevels column name of the sample table indicating the categories to be used as ground truth to assess classifiers
#' @param dLevelSID column name of the sample table that indicates sample id
#' @param resolution increment of thresholds for assessing the classifier
#'
#' @return list of data frames with threshold, sens, precision
#' @export
ccn_classAssess<-function(ct_scores, stVal, classLevels="description2", dLevelSID="sample_id", resolution=0.005) {
  allPRs = list()
  evalAll = matrix(0, nrow=nrow(ct_scores),ncol=2)
  classifications = rownames(ct_scores)
  rownames(stVal) = as.vector(stVal[,dLevelSID])
  i = 1

  for(xname in classifications){
    classification = classifications[i]
    tmpPR = cn_eval(ct_scores[xname,],
                     stVal,
                     classLevels,
                     xname,threshs=seq(0,1, by=resolution), dLevelSID=dLevelSID)
    allPRs[[xname]] = as.data.frame(tmpPR)
    i = 1+i
  }

  return(allPRs)
}

#' @title
#' determine performance of classification at given threshold
#'
#' @description
#' determine performance of classification at given threshold
#' @param vect vector of values
#' @param sampTab sample table
#' @param dLevel colname
#' @param classification actual classification
#' @param thresh threshold above which to make a call
#' @param dLevelSID column to indicate sample id
#'
#' @return return a data frame of the number of TP, FN, FP, and TN
cn_eval<-function(vect, sampTab, dLevel, classification, threshs=seq(0,1,by=0.05), dLevelSID="sample_id"){
  ans = matrix(0,nrow=length(threshs), ncol=9)
  for(i in seq(length(threshs))){
    thresh = threshs[i]
    ans[i,1:4] = cn_clPerf(vect, sampTab, dLevel, classification, thresh, dLevelSID=dLevelSID)
  }
  ans[,5] = threshs
  colnames(ans) = c("TP", "FN", "FP", "TN", "thresh","FPR", "TPR", "Sens", "Prec")
  TPR = ans[,'TP']/(ans[,'TP']+ans[,'FN'])
  FPR = ans[,'FP']/(ans[,'TN']+ans[,'FP'])
  sens = ans[,"TP"]/(ans[,"TP"]+ans[,"FN"])
  prec = ans[,"TP"]/(ans[,"TP"]+ans[,"FP"])
  ans[,'TPR'] = TPR
  ans[,'FPR'] = FPR
  ans[,'Sens'] = sens
  ans[,'Prec'] = prec

  return(ans)
}

#' run cn_clPerf across thresholds
#'
#' run cn_clPerf across thresholds
#' @param vect named vector
#' @param sampTab sample table
#' @param dLevel description level
#' @param classification classification matrix
#' @param thresh seq of pval cutoffs
#' @param dLevelSID column to indicate sample id
#'
#' @return vector of TP FN FP TN
cn_clPerf<-function(vect, sampTab, dLevel, classification, thresh, dLevelSID="sample_id"){
  TP = 0
  FN = 0
  FP = 0
  TN = 0
  sampIDs = names(vect)
  classes = as.vector(sampTab[sampIDs,dLevel])

  actualPos = as.vector(sampTab[sampTab[,dLevel]==classification,dLevelSID])
  actualNeg = setdiff(sampIDs, actualPos)

  calledPos = names(which(vect>thresh))
  calledNeg = names(which(vect<=thresh))

  TP = length(intersect(actualPos, calledPos))
  FP = length(intersect(actualNeg, calledPos))
  FN = length(actualPos)-TP
  TN = length(actualNeg)-FP

  return(c(TP, FN, FP, TN))
}
