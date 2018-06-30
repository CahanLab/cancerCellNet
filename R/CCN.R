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


#' find candidate classifier-worthy genes
#'
#' find candidate classifier-worthy genes
#'
#' @param expDat expDat
#' @param sampTab sampTab
#' @param dLevel dLevel
#' @param topX topX
#' @param dThresh dThresh
#' @param alpha1 alpha1
#' @param alpha2 alpha2
#' @param mu mu
#'
#' @return list of cgenes and grps
#' 
#' @export
findClassyGenes<-function(
	expDat,
	 sampTab,
	 dLevel,
	 topX=25,
	 dThresh=0,
	 alpha1=0.05,
	 alpha2=.001,
	 mu=2)
{
	gsTrain<-sc_statTab(expDat, dThresh=dThresh)
	ggenes<-sc_filterGenes(gsTrain, alpha1=alpha1, alpha2=alpha2, mu=mu)
	grps<-as.vector(sampTab[,dLevel])
	names(grps)<-rownames(sampTab)
	xdiff<-gnrAll(expDat[ggenes,], grps)
	cgenes<-lapply(xdiff, getClassGenes, topX=topX)
	cgenes<-unique(unlist(cgenes))
	list(cgenes=cgenes, grps=grps)
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


#' weighted subtraction from mapped reades, applied to all
#'
#' Simulate expression profile of  _total_ mapped reads
#' @param expRaw matrix of total mapped reads per gene/transcript
#' @param total numeric post transformation sum of read counts
#'
#' @return vector of downsampled read mapped to genes/transcripts
#'
#' @export
weighted_down<-function(
	expRaw,
 	total,
 	dThresh=0
 ){
    expCountDnW<-apply(expRaw, 2, downSampleW, total=total, dThresh=dThresh)
    #log(1+expCountDnW)
    expCountDnW
  }

#' divide each column by sum of that column then scale to xFact and log it
#'
#' divide each column by sum of that column then scale to xFact and log it
#' @param expDat expDat
#' @param xFact xFact
#'
#' @return vector of downsampled read mapped to genes/transcripts
#'
#' @export
trans_prop<-function(
	expDat,
 	xFact=1e5
){
  ans<-matrix(0, nrow=nrow(expDat), ncol=ncol(expDat));
  for(i in seq(ncol(expDat))){
    ans[,i]<-expDat[,i]/sum(expDat[,i]);    
  }
  ans<-ans*xFact;
  colnames(ans)<-colnames(expDat);
  rownames(ans)<-rownames(expDat);
  log(1+ans)
}


#
#' make a classifier, with a randomized class too
#'
#' rmake a classifier, with a randomized class too
#' @param expTrain training data
#' @param genes vector of genes to use as predictors
#' @param groups named vector of cells to groups or classes
#' @param nRand =50 num of randomized profiles to make
#' @param ntrees =2000 number of trees to build

#' @return RF
#' @export
#'
makeClassifier<-function(
  expTrain,
  genes,
  groups,
  nRand=50,
  ntrees=2000){


	randDat<-randomize(expTrain, num=nRand)
	expTrain<-cbind(expTrain, randDat)

	allgenes<-rownames(expTrain)

	missingGenes<-setdiff(unique(genes), allgenes)
	cat("Number of mussing genes ", length(missingGenes),"\n")
	ggenes<-intersect(unique(genes), allgenes)
	randomForest(t(expTrain[ggenes,]), as.factor(c(groups, rep("rand", ncol(randDat)))), ntree=ntrees)

}

#
#' classify samples
#'
#' classify samples
#' @param rfObj result of running sc_makeClassifier
#' @param expQuery expQuery
#' @param numRand numRand

#' @return classRes matrix
#' @export
#'
rf_classPredict<-function(
  rfObj,
  expQuery,
  numRand=50){

  	randDat<-randomize(expQuery, num=numRand)
  	expQuery<-cbind(expQuery, randDat)

    preds<-rownames(rfObj$importance)
  	xpreds<-t(predict(rfObj, t(expQuery[preds,]), type='prob'))
	colnames(xpreds)<-colnames(expQuery)
	xpreds
}


#
#' randomize data matrix 
#'
#' randomize data matrix 
#' @param expDat expDat
#' @param num number of profiles to return
#'
#' @return exp matrix random
#' @export
#'
randomize<-function(
 expDat,
 num=50){


	randDat<-t(apply(expDat, 1, sample))	
 	randDat<-apply(randDat, 2, sample)

 	randDat<-randDat[,sample(1:ncol(randDat), num)]
	colnames(randDat)<-paste0(rep("rand_", num), 1:num)
	rownames(randDat)<-rownames(expDat)
	randDat
}



#' makes complete gene-to-gene comparison
#'
#' @param expDat expDat
#' @param genePairs genePairs
#'
#' @return matrix indicating which gene of a pair is greater
#'
#' @export
query_transform<-function( # convert to a vector of length = length(vect)^2 - 1 /2
	expDat,
 	genePairs #vector of strings indicating pairs to compare
 ){
	genes<-strsplit(genePairs, "_")
    ans<-matrix(0, nrow=length(genes), ncol=ncol(expDat))
	pair_index<-1
	genes1<-vector()
	genes2<-vector()
	for(i in seq(length(genes))){
		ans[i,]<-as.numeric(expDat[genes[[i]][1],]>expDat[genes[[i]][2],])
	}
	colnames(ans)<-colnames(expDat)
	rownames(ans)<-genePairs
	ans
}



