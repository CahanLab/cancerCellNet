# cancerCellNet
# (C) Patrick Cahan 2012-2018

# commonly used or misc functions

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

# replavce NAs with 0
repNA<-function
(vector){
  vector[which(is.na(vector))]<-0;
  vector;
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


#' @title
#' Convert Long Expression Table
#'
#' @description
#' Convert long gene expression/count table into wide dataframe to proceed to the next step
#' of training classifier or validating samples.
#'
#' @param longDf a dataframe of gene expression/count in long table format
#' @param SampCol the name of the column containing sample names
#' @param geneCol the name of the column containing gene names
#' @param geneExpcol the name of the column containing gene counts or expression
#' @param geneNaOmit TRUE if you want to omit the genes with NA
#'
#' @return a list containing the gene expression dataframe in wide table format and a dataframe of omitted genes if geneNaOmit is TRUE
#' Otherwise, just returns gene Expression matrix in wide table
#' @export
#'
#' @importFrom stringr str_replace
#'
convertLongTab <- function(longDf, SampCol, geneCol, geneExpCol, geneNaOmit = TRUE) {
  longDf_extract <- data.frame(longDf[, c(SampCol, geneCol, geneExpCol) ])
  exp_tab_rs <- reshape(longDf_extract, idvar = geneCol, timevar = SampCol, direction = "wide")
  colnames(exp_tab_rs) <- stringr::str_replace(colnames(exp_tab_rs), paste0(geneExpCol, "."), "")
  rownames(exp_tab_rs) = exp_tab_rs[, geneCol]
  exp_tab_rs[, geneCol] = NULL

  exp_tab_omitted <- data.frame()

  if (geneNaOmit == TRUE) {
    exp_tab_clean = na.omit(exp_tab_rs);

    exp_tab_omitted = exp_tab_rs[!(rownames(exp_tab_rs) %in% rownames(exp_tab_clean)), ]

    output = list(wideTable = exp_tab_clean, omittedGenesTable = exp_tab_omitted);

  }
  else {
    output = exp_tab_rs
  }

  output
}
