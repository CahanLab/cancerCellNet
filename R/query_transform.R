
#' Make complete gene-to-gene comparison
#'
#' This function will compare the genes in genePairs vector sample by sample in the expression matrix.
#'
#' @param expDat the expression matrix
#' @param genePairs a vector with gene pairs
#'
#' @return A matrix indicating whether the first gene in the gene pair has a greater expression than the second gene in the gene pair
#'
#' @export
query_transform <- function(expDat, genePairs) {
  genes<-strsplit(genePairs, "_")
  ans<-matrix(0, nrow=length(genes), ncol=ncol(expDat))
  pair_index <- 1
  genes1<-vector()
  genes2<-vector()

  for(i in seq(length(genes))){
    ans[i,]<-as.numeric(expDat[genes[[i]][1],]>expDat[genes[[i]][2],])
  }
  colnames(ans)<-colnames(expDat)
  rownames(ans)<-genePairs
  ans
}
