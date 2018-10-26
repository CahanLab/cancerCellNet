# cancerCellNet
# (C) Patrick Cahan 2016-2018


#' find best pairs
#'
#' find best pairs
#'
#' @param expDat expDat
#' @param cellLabels named vector of cell groups
#'
#' @return vector of pairs
#'
#' @export
gnrBP<-function(
  expDat,
  cellLabels,
  topX=50){

	ans<-vector()
    myPatternG<-sc_sampR_to_pattern(as.character(cellLabels))
    for(i in seq(length(myPatternG))){
    	cat(i,"\n")
    	xres<-sc_testPattern(myPatternG[[i]], expDat=expDat)
    	tmpAns<-findBestPairs(xres, topX)
    	ans<-append(ans, tmpAns)
    }
    unique(ans)
}

#' makes complete gene-to-gene comparison
#'
#' @param expDat expDat
#'
#' @return list of list(genes=data frame, expDat=binary matrix)
#'
#' @export
pair_transform<-function( # convert to a vector of length = length(vect)^2 - 1 /2
expDat){
	ngenes<-nrow(expDat)
	genes<-rownames(expDat)
	ans<-matrix(0, nrow=ngenes*(ngenes-1)/2, ncol=ncol(expDat))
	pair_index<-1
	genes1<-vector()
	genes2<-vector()
	for(i in 1:ngenes){
		for(j in 1:ngenes){
			if(j>i){
				genes1<-append(genes1, genes[i])
				genes2<-append(genes2, genes[j])
				ans[pair_index,]<-as.numeric(expDat[i,]>expDat[j,])
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







