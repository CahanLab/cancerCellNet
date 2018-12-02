# analysis and non-generalized functions for CCN project

#' Sum transcription expression estimates to gene-level expression measures
#'
#' Sums all trnascript level expression estimates to a single gene level estimate. Needs a table that maps transcript IDs to gene IDs.
#' @param expRaw salmon isoform estimates
#' @param numCores num of cores to use for parallel 
#' @param geneTabfname gene <-> transcript mapping, df that needs cols: gene_id, transcript_id
#' @param nameCol gene ann table column name to average over
#'
#' @return matrix (ansCounts)
gene_expr_sum<-function
(expRaw,
 numCores=3,
 geneTabfname="./geneToTrans_Mus_musculus.GRCm38.80.exo_Jun_02_2015.R",
 nameCol="gene_name"
){
  
  matchFunc<-function(val, vect){
    which(vect==val)
  }
  
  save_colnames = colnames(expRaw)
  expCounts<-expRaw
  if(!is.matrix(expCounts)){
    expCounts<-as.matrix(expCounts);
  }
  
  geneTab<-utils_loadObject(geneTabfname);
  
  # keep only common variables
  sameProbes<-intersect(rownames(expCounts), rownames(geneTab));
  expCounts<-as.matrix(expCounts[sameProbes,]);
  geneTab<-geneTab[sameProbes,];
  
  # all genes
  allgenes<-as.vector(geneTab[,nameCol]);
  
  # unique genes
  eids<-unique(allgenes)
  
  # make a cluster
  aClust<-parallel::makeCluster(numCores, type='SOCK')
  
  # a list of indices into all genes
  geneIndexList<-parLapply(aClust, eids, matchFunc, vect=allgenes)
  names(geneIndexList)<-eids
  uSymbols<-vector(length=length(eids))
  
  stopCluster(aClust)
  

  ansCounts<-matrix(0, nrow=length(eids), ncol=ncol(expCounts))

  for(i in seq(length(geneIndexList))){
    eid<-eids[i]
    #cat(".");
    xi <-  geneIndexList[[i]];
    
    ## desProbes <- as.vector(geneTab[xi,]$probe_id);
    desProbes <- as.character(geneTab[xi,]$transcript_id);
    if(length(xi)>1){
      ansCounts[i,]<-apply(as.matrix(expCounts[desProbes,]), 2, sum);
    }
    else{
      ansCounts[i,]<-as.matrix(expCounts[desProbes,]);
    }
    uSymbols[i]<-as.vector(geneTab[ xi[1] ,nameCol]);
  }

  rownames(ansCounts)<-uSymbols
  colnames(ansCounts) = save_colnames
  ansCounts


}