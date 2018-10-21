# cancerCellNet
# (C) Patrick Cahan 2012-2018

# commonly used or misc functions

#' find genes higher in a cluster compared to all other cells
#'
#' ind genes higher in a cluster compared to all other cells
#'
#' @param expDat a matrix of normalized gene expression taken from \code{\link(trans_prop)}
#' @param cellLabels a vector of all the named vector of cancer categories
#'
#' @return list of dataFrames containing pval, cval and holm for each gene in each cancer category
#'
#' @export
gnrAll<-function(expDat, cellLabels){

  myPatternG<-sc_sampR_to_pattern(as.character(cellLabels))
  specificSets<-lapply(myPatternG, sc_testPattern, expDat=expDat)
  cat("Done testing\n")

#  grpOrder<-myGrpSort(cellLabels)

#  specificSets[grpOrder]

  specificSets
}


#'
#' @param sampR a vector of the sample catgories
#'
#' @return a list of vectors indicating mapping of each cancer type in the sampR vector
#' @export
sc_sampR_to_pattern<-function(sampR){
  d_ids<-unique(as.vector(sampR));
  nnnc<-length(sampR);
#  ans<-matrix(nrow=length(d_ids), ncol=nnnc);
  ans<-list()
  for(d_id in d_ids){
    x<-rep(0,nnnc);
    x[which(sampR==d_id)]<-1;
    ans[[d_id]]<-x;
  }
  ans
}



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

#' @title
#' Find Genes that Pass Criteria
#'
#' @description
#' Based on idea that reliably detected genes will either be detected in many cells, or highly expressed in a small number of cells or both
#'
#' @param geneStats a matrix containing stats of genes generated from running \code{\link{sc_statTab}}
#' @param alpha1 a number representing proportion of cells in which a gene must be considered detected (as defined in geneStats)
#' @param alpha2 a number representing lower proportion of cells for genes that must have higher expression level
#' @param mu a number represeting threshold for average expression level of genes passing the lower proportion criteria
#'
#' @return a vector of gene symbols
#' @export
#'
sc_filterGenes<-function(geneStats, alpha1=0.1, alpha2=0.01, mu=2){
  passing1<-rownames(geneStats[geneStats$alpha>alpha1,])
  notPassing<-setdiff(rownames(geneStats), passing1)
  geneStats<-geneStats[notPassing,]
  c(passing1, rownames(geneStats[which(geneStats$alpha>alpha2 & geneStats$mu>mu),]))
}


#' find cells that pass criteria
#'
#' based purely on umis
#'
#' @param sampTab, which must have UMI column
#' @param minVal umis must exceed this
#' @param maxValQuant quantile to select max threshold
#'
#' @return vector rownames(sampTab) meeting criteria
#'
#' @export
#'
sc_filterCells<-function
(sampTab,
 minVal=1e3,
 maxValQuant=0.95){
  stX<-sampTab[sampTab$umis>minVal,]
  qThresh<-quantile(sampTab$umis, maxValQuant)
  rownames(stX[stX$umis<qThresh,])
}




#' @title
#' Generate a Stats Table
#' @description
#' To make a stats table with alpha, mu, overall mean, coefficient of variance, fano factor, max value, standard deviation
#' @param expDat a matrix of gene expression values generated from \code{\link{trans_prop}}
#' @param dThresh a number indication the threshold for detection
#'
#' @return a dataframe containing alpha, mu, overall mean, coefficient of variance, fano factor, max value,
#' and standard deviation of the expression matrix
#'
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
  indexFunction<-function(vector, threshold){
    names(which(vector>threshold));
  }

  indexes<-apply(expMat, 1, indexFunction, threshold);

  alphas<-unlist(lapply(indexes, length));

  ans<-alphas/ncol(expMat)

  if(pseudo){
    ans<-(alphas+1)/(ncol(expMat)+1)
  }
  ans
}

#' @title
#' compute Mu
#' @description
#' compute Mu given threshold of the expression matrix
#' @param expMat a matrix with gene expressions data
#' @param threshold a number indicating the threshold at which gene expression below does not count
#'
#' @return a list of Mus for individual gene in the expression matrix
#' @export
sc_compMu<-function(expMat, threshold=0){

  afunct<-function(vector, threshold){
    mean( vector[which(vector>threshold)] );
  }

  mus<-unlist(apply(expMat, 1, afunct, threshold))
  mus[is.na(mus)]<-0;
  mus;
}

# replavce NAs with 0
repNA<-function
(vector){
  vector[which(is.na(vector))]<-0;
  vector;
}

#' @title
#' Compute Fano Factor
#' @description
#' Compute fano factor for vector
#' @param vector a vector
#' @return the fano factor
sc_fano<-function
(vector){
  var(vector)/mean(vector);
}

# compute coeef of variation on vector
#' @title Compute Coefficient of Variation
#' @description Compute the coefficient of variation of a vector
#' @param vector a vector
#'
#' @return the coefficient of variation
sc_cov<-function
(vector){
  sd(vector)/mean(vector);
}


#' Weighted subtraction from mapped reades
#'
#' Simulate expression profile of  _total_ mapped reads as a way to normalize the counts
#'
#' @param vector a vector of total mapped reads per gene/transcript
#' @param total the sum of the read counts after transformation
#' @param dThresh the threshold at which anything lower is 0 after transformation. Usually 0
#'
#' @return vector of downsampled read mapped to genes/transcripts
#' @export
downSampleW<-function(vector, total=1e5, dThresh=0){

  totalSignal<-sum(vector) # get the sum of the vector
  wAve<-vector/totalSignal #

  resid<-totalSignal-total #num to subtract from sample
  residW<-wAve*resid # amount to substract from each gene

  ans<-vector-residW
  ans[which(ans<dThresh)]<-0
  ans
}


#' Weighted subtraction from mapped reades, apply to all columns
#'
#' Simulate expression profile of  _total_ mapped reads
#' @param expRaw matrix of total mapped reads per gene/transcript
#' @param total numeric post transformation sum of read counts
#' @param dThresh the threshold at which anything lower than that is 0
#'
#' @return matrix of downsampled read mapped to genes/transcripts
#'
#' @export
weighted_down<-function(expRaw, total=1e5, dThresh=0){
    expCountDnW<-apply(expRaw, 2, downSampleW, total=total, dThresh=dThresh)
    expCountDnW
  }


#' @title
#' Split Sample Table
#'
#' @description
#' Split a sample table into training set and validation set.
#' @param sampTab sample table (DataFrame)
#' @param ncells number of samples for training in each category (Integer)
#' @param dLevel the column name with the classification categories (String)
#' @return a list containing training sample table and validation sample table
#' @examples
#' stList<-splitCommon(stTrain, ncells=25, dLevel="description1")
#' @export
splitCommon<-function(sampTab, ncells, dLevel="description1"){

  cts<-unique(as.vector(sampTab[,dLevel])) #receive the names of the categories
  trainingids<-vector()

  for(ct in cts){
    cat(ct,": ")
    stX<-sampTab[sampTab[,dLevel]==ct,]
    ccount<-nrow(stX)-3
    ccount<-min(ccount, ncells)
    cat(nrow(stX),"\n")
    trainingids<-append(trainingids, sample(rownames(stX), ccount)) # randomly samples ccount of training samples
  }

  val_ids<-setdiff(rownames(sampTab), trainingids) # the samples that are not used to training are used for validation
  list(train=sampTab[trainingids,], val=sampTab[val_ids,])
}


#' @export
getGenesFromGO<-function# return the entrez gene ids of a given a GOID, for now assumes mouse
(GOID, # GO id to find genes for
 annList
){
  sort(as.vector(unlist(annList[['egSymbols']][annList[['goegs']][[GOID]]])));
}



#' row average (or median) based on groups
#'
#' row average (or median) based on groups
#' @param exp expression df
#' @param groupings groupings
#' @param type mean or media
#'
#' @return return a dataframe of mean or median-ed data based on given groupings.  colnames become the column name of the first sample in each group from the original data
#'
#' @export
GEP_makeMean<-function
(exp,
 groupings,
 type='mean'
){


  ans<-data.frame();
  grps<-unique(groupings);
  if(type=='mean'){
    for(grp in grps){
      gi<-which(groupings==grp);
      if(length(gi)==1){

        if(nrow(ans)==0){
          ans<-data.frame(exp[,gi]);
        }else{
          ans<-cbind(ans, exp[,gi]);
        }
      }
      else{
        xxx<-apply(exp[,gi],1,mean);
        if(nrow(ans)==0){
          ans<-data.frame(xxx);
        }
        else{
          ans<-cbind(ans, xxx);
        }
      }
    }
  }
  else{
    for(grp in grps){
      gi<-which(groupings==grp);
      xxx<-apply(exp[,gi],1,median);
      if(nrow(ans)==0){
        ans<-data.frame(xxx);
      }
      else{
        ans<-cbind(ans, xxx);
      }
    }
  }

  colnames(ans)<-grps;
  ans;
  ### data.frame of mean or median-ed data based on given groupings
}


#' find transcript factors
#'
#' find transcript factors
#' @param annotation
#' @param species defaul is 'Hs', can also be 'Mm;
#' @param ontology default is BP
#'
#' @return vector fo TF names
#' @export
#' @importFrom AnnotationDbi as.list
#'
find_genes_byGo<-function#
(annotation,
  species='Hs',
  onto="BP"
){

  cat("Loading gene annotations ...\n")
  require(GO.db);

  if(species=='Hs'){
    require(org.Hs.eg.db);
    egSymbols<-as.list(org.Hs.egSYMBOL);
    goegs<-as.list(org.Hs.egGO2ALLEGS);
  }
  else{
    require(org.Mm.eg.db);
    egSymbols<-as.list(org.Mm.egSYMBOL);
    goegs<-as.list(org.Mm.egGO2ALLEGS);
  }

  goterms<-as.list(GOTERM);
  goids<-names(goegs);
  onts<-lapply(goids, Ontology);
  bps<-onts[onts==onto];
  goids<-names(unlist(bps));

  cat("matching gene symbols and annotations")
  gobpList<-list();
  for(goid in goids){
    egs <- goegs[[ goid ]];
    goterm<-Term(goterms[[goid]]);
    genes<-sort(unique(as.vector(unlist( egSymbols[egs] ))));
    gobpList[[goterm]]<-genes;
  }

  ### newHsTRs<-gobpList[['regulation of transcription, DNA-dependent']];
  regNames<-names(gobpList)[grep(annotation, names(gobpList))];
  trs<- unique(unlist(gobpList[regNames]));
  cat(annotation, ": ", length(trs),"\n");
  sort(trs)

}




#' 1-PCC distance
#'
#' 1-PCC distance
#' @param x numeric matrix
#'
#' @return distance matrix
#'
#' @examples
#' xdist<-utils_myDist(t(expDat))
#' plot(hclust(xdist, 'ave'), hang=-1)
#'
#' @export
utils_myDist<-function
(x
){
  as.dist(1-cor(t(x)));
}

#' Load an R object
#'
#' Load an R object
#' @param fname The name of the R object
#'
#' @return A R object
#' @example
#' utils_loadObject("ccn_classifier_Jun_29_2018.rda")
#'
#' @export
utils_loadObject<-function
(fname
 ### file name
){
  x<-load(fname);
  get(x);
}

#' strip whitespace from a string
#'
#' strip whitespace from a string
#' @param string string
#'
#' @return new string
#'
#' @export
utils_stripwhite<-function
###
(string
 #### string
 ){
  gsub("^\\s+|\\s+$", "", string)
}

#' print date
#'
#' print date
#' @return string
#'
#' @export
utils_myDate<-function
###
()
{
  format(Sys.time(), "%b_%d_%Y");
}

#' reduces full path to filename
#'
#' reduces full path to filename
#' @param string
#'
#' @return something
#'
#' @export
utils_strip_fname<-function #
(str){
  a<-strsplit(str, "/")[[1]];
  a[length(a)];
}

utils_stderr<-function
### calculate standard error
(x){
  sqrt(var(x)/length(x));
  ### stderr
}

zscore<-function
### compute zscore
(x,
 ### numeric vector
 meanVal,
 ### mean of distribution to compute zscore of x against
 sdVal
 ### standard deviation of distribution to compute zscore of x agains
 ){
  (x-meanVal)/sdVal;
  ### zscore
}


zscoreVect<-function
### Compute the mean zscore of given genes in each sample
(genes,
 ### genes
 xvals,
 ### named vector
 tVals,
 ### tvals
 ctt
 ### ctt
 ){
  ans<-vector();
  for(gene in genes){
    ans<-append(ans, zscore(xvals[gene], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]]));
  }
  ans;
  ### zscore vector
}

#' make Inf and -Inf values sensible
#'
#' make Inf and -Inf values sensible
#' @param zMat zMat
#'
#' @return corrected zMat
#'
#' @export
cn_correctZmat<-function
(zmat){
  myfuncInf<-function(vect){
    xi<-which(vect=='Inf')
    if(any(xi)){
      mymax<-max(vect[-xi])
      vect[xi]<-mymax
    }
    vect
  }
  zmat<-apply(zmat,2, myfuncInf)
  zmat[is.na(zmat)]<-0
  zmat
}
