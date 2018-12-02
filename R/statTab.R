#' @title
#' Generate a Stats Table
#' @description
#' To make a stats table with alpha, mu, overall mean, coefficient of variance, fano factor, max value, standard deviation
#' @param expDat a matrix of gene expression values generated from \code{\link{trans_prop}}
#' @param dThresh a number indication the threshold for detection
#'
#' @return a dataframe containing alpha, mu, overall mean, coefficient of variance, fano factor, max value,
#' and standard deviation of the expression matrix
#' @export
sc_statTab<-function(expDat, dThresh=0){
  statTab<-data.frame()
  # below generates various statistical values
  muAll<-sc_compMu(expDat, threshold=dThresh);
  alphaAll<-sc_compAlpha(expDat,threshold=dThresh);
  meanAll<-apply(expDat, 1, mean);
  covAll<-apply(expDat, 1, sc_cov);
  fanoAll<-apply(expDat,1, sc_fano);
  maxAll<-apply(expDat, 1, max);
  sdAll<-apply(expDat, 1, sd);

  statTabAll<-data.frame(gene=rownames(expDat), mu=muAll, alpha=alphaAll, overall_mean=meanAll, cov=covAll, fano=fanoAll, max_val=maxAll, sd=sdAll)
  statTabAll;
}

#' @title
#' compute Mu
#' @description
#'
#' compute Mu given threshold of the expression matrix
#' @param expMat a matrix with gene expressions data
#' @param threshold a number indicating the threshold at which gene expression below does not count
#'
#' @return a list of Mus for individual gene in the expression matrix
#' @export
sc_compMu<-function(expMat, threshold=0){

  afunct<-function(vector, threshold){
    mean(vector[which(vector>threshold)]);
  }

  mus<-unlist(apply(expMat, 1, afunct, threshold))
  mus[is.na(mus)]<-0;
  mus;
}


#' @title
#' Compute Alpha
#' @description
#' Compute Alpha of expression matrix given detection threshold
#'
#' @param expMat a matrix containing gene expression data
#' @param threshold a number indicating the detection threshold
#' @param pseudo logical indicating if it is going to be pseudo alpha or real alpha
#'
#' @return a list of alphas for the genes in expression matrix
#' @export
sc_compAlpha<-function(expMat, threshold=0,pseudo=FALSE){
  # identify the index of vectors greater than threshold
  lengthFunction<-function(vector, threshold){
    length(which(vector>threshold));
  }

  alphas <- apply(expMat, 1, lengthFunction, threshold)

  ans <- alphas/ncol(expMat)

  if (pseudo) {
    ans <- (alphas+1) / (ncol(expMat) + 1)
  }

  #return
  ans
}

#' @title
#' Compute Fano Factor
#' @description
#' Compute fano factor for vector
#' @param vector a vector
#' @return the fano factor
#' @export
sc_fano<-function(vector){
  var(vector)/mean(vector);
}

# compute coeef of variation on vector
#' @title Compute Coefficient of Variation
#' @description Compute the coefficient of variation of a vector
#' @param vector a vector
#'
#' @return the coefficient of variation
#' @export
sc_cov<-function(vector){
  sd(vector)/mean(vector);
}
