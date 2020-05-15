#' @title
#' Log scaled downsampled expression
#'
#' @description
#' Divide each column by sum of that column then scale to xFact and log it for more normalization.
#'
#' @param expDat a matrix of weighted down expression data
#' @param xFact a number representing scaling factor
#'
#' @return a matrix of downsampled read mapped to genes/transcripts
#'
#' @export
trans_prop<-function(expDat, xFact=1e5){
  ans = matrix(0, nrow=nrow(expDat), ncol=ncol(expDat))

  for(i in seq(ncol(expDat))){
    ans[,i] = expDat[,i]/sum(expDat[,i])
  }

  ans = ans*xFact
  colnames(ans) = colnames(expDat)
  rownames(ans) = rownames(expDat)

  return(log(1+ans))
}
