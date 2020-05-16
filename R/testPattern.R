#' @title
#' Template Matching
#'
#' @description
#' Match the template created \code{\link{sc_sampR_to_pattern}} with expression data
#'
#' @param pattern the template pattern created from \code{\link{sc_sampR_to_pattern}}
#' @param expDat normalized expression data
#'
#' @return a dataframe of genes with pval, cval and holms
#' @export
sc_testPattern<-function(pattern, expDat){
  pval = vector()
  cval = vector()

  geneids = rownames(expDat)  # get the gene names of the exp matrix

  llfit = ls.print(lsfit(pattern, t(expDat)), digits=25, print=FALSE)

  xxx = matrix( unlist(llfit$coef), ncol=8,byrow=TRUE)

  # remove r squared values that are negative
  naIndex = (llfit$summary[, 2] < 0)

  ccorr = xxx[!naIndex,6]  # t-value of X coefficient
  cval =  sqrt(as.numeric(llfit$summary[!naIndex,2])) * sign(ccorr)  # R squared
  pval = as.numeric(xxx[!naIndex,8])  # p-value of X coefficient

  holm = p.adjust(pval, method='holm')


  return(data.frame(row.names=geneids[!naIndex], pval=pval, cval=cval,holm=holm))
}
