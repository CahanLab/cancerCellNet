#' @title
#' Gene to Gene comparison
#'
#' @description
#' Make complete gene to gene comparison for all the classification worthy genes among each sample
#'
#' @param expDat raw count expression table with cgenes
#'
#' @return a dataframe with individual gene pairs on rows and samples as columns
#' @export
pair_transform<-function(expDat){
  ngenes<-nrow(expDat) #get number of genes
  genes<-rownames(expDat)

  ans<-matrix(0, nrow=ngenes*(ngenes-1)/2, ncol=ncol(expDat)) # same columns and n choose 2 combination for rows

  pair_index<-1
  genes1<-vector()
  genes2<-vector()

  for(i in 1:ngenes){
    for(j in 1:ngenes){
      if(j>i){
        genes1<-append(genes1, genes[i])
        genes2<-append(genes2, genes[j])
        ans[pair_index,]<-as.numeric(expDat[i,]>expDat[j,]) # will give a vector as long as column size of 1 and 0
        pair_index<-pair_index +1
      }
    }
  }

  colnames(ans)<-colnames(expDat)
  tList2 <- list(genes=data.frame(g1=genes1, g2=genes2), tDat=ans)

  pairNames<-paste(tList2[[1]][,1], "_",tList2[[1]][,2], sep='')

  pairDat<-tList2[[2]]
  rownames(pairDat)<-pairNames
  pairDat
}
