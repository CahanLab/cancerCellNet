#' @export
sc_testPattern<-function(pattern, expDat){
  pval<-vector();
  cval<-vector();

  geneids<-rownames(expDat); # get the gene names of the exp matrix

  llfit<-ls.print(lsfit(pattern, t(expDat)), digits=25, print=FALSE);

  xxx<-matrix( unlist(llfit$coef), ncol=8,byrow=TRUE);
  ccorr<-xxx[,6]; # t-value of X coefficient
  cval<- sqrt(as.numeric(llfit$summary[,2])) * sign(ccorr); # R squared
  pval<-as.numeric(xxx[,8]); # p-value of X coefficient

  #qval<-qvalue(pval)$qval;
  holm<-p.adjust(pval, method='holm');
  #data.frame(row.names=geneids, pval=pval, cval=cval, qval=qval, holm=holm);
  data.frame(row.names=geneids, pval=pval, cval=cval,holm=holm);
}
