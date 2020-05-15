#' @title
#' Make complete gene-to-gene comparison
#' @description
#' Performs gene to gene comparison of each gene pairs in the expression table.
#' @param expDat the expression matrix
#' @param genePairs a vector with gene pairs
#' @return A matrix of the gene to gene comparison results. 1 indicates the first gene has higher expression than the second gene in the gene pair. 0 indicates the
#' second gene has higher expression than the first gene in the gene pair.
#' @export
query_transform <- function(expDat, genePairs) {

  expDat = as.matrix(expDat) # additional, in case that expDat is in dataframe which is very likely
  genes = strsplit(genePairs, "_") # split the gene pairs into seperate genes
  ans = matrix(0, nrow=length(genes), ncol=ncol(expDat))

  for(i in seq(length(genes))){
    ans[i,] = as.numeric(expDat[genes[[i]][1],]>expDat[genes[[i]][2],])
  }
  colnames(ans) = colnames(expDat)
  rownames(ans) = genePairs

  return(ans)
}
