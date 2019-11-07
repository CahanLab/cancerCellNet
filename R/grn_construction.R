#' find transcript factors
#'
#' find transcript factors
#' @param species defaul is 'Hs', can also be 'Mm;
#'
#' @return vector fo TF names
#'
#' @import GO.db, org.Hs.eg.db
#' @export
#'
find_tfs<-function(species='Hs'){

  cat("Loading gene annotations ...\n")
  require(GO.db);

  if(species=='Hs'){
    require(org.Hs.eg.db);
    egSymbols<-as.list(org.Hs.egSYMBOL);
    goegs<-as.list(org.Hs.egGO2ALLEGS);
  }

  goterms<-as.list(GOTERM);
  goids<-names(goegs);
  onts<-lapply(goids, Ontology);
  bps<-onts[onts=='BP'];
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
  regNames<-names(gobpList)[grep("regulation of transcription", names(gobpList))];
  trs<- unique(unlist(gobpList[regNames]));
  cat("Regulation of transcription: ", length(trs),"\n");

  mfs<-onts[onts=='MF'];
  goidsMF<-names(unlist(mfs));

  gomfList<-list();
  for(goid in goidsMF){
    egs <- goegs[[ goid ]];
    goterm<-Term(goterms[[goid]]);
    genes<-sort(unique(as.vector(unlist( egSymbols[egs] ))));
    gomfList[[goterm]]<-genes;
  }
  dbs<-gomfList[['DNA binding']];
  cat("DNA binding: ", length(dbs),"\n");
  sort(intersect(trs, dbs));
}

#' gene-gene correlations, and round
#'
#' gene-gene correlations, and round
#' @param expDat expression matrix
#'
#' @return correlation matrix
grn_corr_round<-function(expDat) {
  corrX<-cor(t(expDat)); # pearson
  round(corrX, 3);
}

#' compute context dependent zscores
#'
#' slightly modidied from JJ Faith et al 2007
#' @param corrMat correlation matrix
#'
#' @return matrix of clr zscores
#'
mat_zscores<-function(corrMat){
  corrMat<-abs(corrMat);
  zscs_2<-round(scale(corrMat), 3);
  rm(corrMat);
  gc()
  zscs_2 + t(zscs_2);
}

#' extracts the TRs, zscores, and corr values passing thresh
#'
#' extracts the TRs, zscores, and corr values passing thresh
#' @param zscores, # zscore matrix, non-TFs already removed from columns
#' @param corrMatrix, # correlation matrix
#' @param genes, # vector of target genes
#' @param threshold # zscore threshold
#'
#' @return data.frame(target=targets, reg=regulators, zscore=zscoresX, corr=correlations);
cn_extractRegsDF<-function(zscores, corrMatrix, genes, threshold){

  targets<-vector();
  regulators=vector();
  zscoresX<-vector();
  correlations<-vector();

  targets<-rep('', 1e6);
  regulators<-rep('', 1e6);
  zscoresX<-rep(0, 1e6);
  correlations<-rep(0, 1e6);

  str<-1;
  stp<-1;
  for(target in genes){
    x<-zscores[target,];
    regs<-names(which(x>threshold));
    if(length(regs)>0){
      zzs<-x[regs];
      corrs<-corrMatrix[target,regs];
      ncount<-length(regs);
      stp<-str+ncount-1;
      targets[str:stp]<-rep(target, ncount);
      #    targets<-append(targets,rep(target, ncount));
      regulators[str:stp]<-regs;
      #regulators<-append(regulators, regs);
      #    zscoresX<-append(zscoresX, zzs);
      zscoresX[str:stp]<-zzs;
      correlations[str:stp]<-corrs;
      str<-stp+1;
    }
    #    correlations<-append(correlations, corrs);
  }
  targets<-targets[1:stp];
  regulators<-regulators[1:stp];
  zscoresX<-zscoresX[1:stp];
  correlations<-correlations[1:stp];


  data.frame(target=targets, reg=regulators, zscore=zscoresX, corr=correlations);
}

#' compute CLR-like zscores
#'
#' compute CLR-like zscores
#' @param corrVals correlation matrix
#' @param tfs vector of  transcriptional regualtor names
#'
#' @return zscore matrix
grn_zscores<-function (corrVals,tfs){
  zscs<-mat_zscores(corrVals);
  gc();
  zscs[,tfs];
}

#' convert a table to an igraph
#'
#' convert a table to an igraph. This adds an nEnts vertex attribute to count the number of entities in the sub-net
#'
#' @param grnTab table of TF, TF, maybe zscores, maybe correlations
#' @param simplify false
#' @param directed FALSE,
#' @param weights TRUE
#'
#' @return iGraph object
ig_tabToIgraph<-function(grnTab, simplify=FALSE, directed=FALSE, weights=TRUE){

  tmpAns<-as.matrix(grnTab[,c("TF", "TG")]);
  regs<-as.vector(unique(grnTab[,"TF"]));
  ###cat("Length TFs:", length(regs), "\n");
  targs<-setdiff( as.vector(grnTab[,"TG"]), regs);

  ###  cat("Length TGs:", length(targs), "\n");
  myRegs<-rep("Regulator", length=length(regs));
  myTargs<-rep("Target", length=length(targs));

  types<-c(myRegs, myTargs);
  verticies<-data.frame(name=c(regs,targs), label=c(regs,targs), type=types);

  ### iG<-graph.data.frame(tmpAns,directed=directed,v=verticies);
  iG<-igraph::graph_from_data_frame(tmpAns,directed=directed,v=verticies);

  if(weights){
    #E(iG)$weight<-grnTab$weight;
    E(iG)$weight<-grnTab$zscore;
  }

  if(simplify){
    iG<-simplify(iG);
  }
  V(iG)$nEnts<-1;
  iG;
}

#' get raw GRN from zscores, and corr
#'
#' get raw GRN from zscores, and corr
#' @param zscores zscores matrix
#' @param corrs correlation matrix
#' @param targetGenes target genes
#' @param zThresh zscore threshold
#'
#' @return list of grnTable and corresponding graph
cn_getRawGRN<-function(zscores, corrs, targetGenes, zThresh=4){

  # make a grn table
  cat("Making GRN table...\n")
  grn<-cn_extractRegsDF(zscores, corrs, targetGenes, zThresh);
  cat("Done making GRN table...\n")
  colnames(grn)[1:2]<-c("TG", "TF");

  # make an iGraph object and find communities
  cat("Make iGraph...\n")
  igTmp<-ig_tabToIgraph(grn, directed=FALSE, weights=TRUE);
  cat("done with iGraph...\n")
  list(grnTable=grn, graph=igTmp);
}

#' find genes that are preferentially expressed in specified samples
#'
#' find genes that are preferentially expressed in specified samples
#' @param expDat expression matrix
#' @param sampTab sample table
#' @param holm sig threshold
#' @param cval R thresh
#' @param dLevel annotation level to group on
#' @param prune boolean limit to genes exclusively detected as CT in one CT
#'
#' @return a list of something
#'
#' @export
#'
cn_findSpecGenes<-function (expDat, sampTab, holm=1e-50, cval=0.5, dLevel="description1", prune=FALSE){

  myPatternG<-cn_sampR_to_pattern(as.vector(sampTab[,dLevel]));
  specificSets<-apply(myPatternG, 1, cn_testPattern, expDat=expDat);

  # adaptively extract the best genes per lineage
  cvalT<-vector();
  ctGenes<-list();
  ctNames<-unique(as.vector(sampTab[,dLevel]));
  for(ctName in ctNames){
    x<-specificSets[[ctName]];
    tmp<-rownames(x[which(x$cval>cval),]);
    tmp2<-rownames(x[which(x$holm<holm),]);
    tmp<-intersect(tmp, tmp2)
    ctGenes[[ctName]]<-tmp;
    ###    cvalT<-append(cvalT, cval);
  }

  if(prune){
    # now limit to genes exclusive to each list
    specGenes<-list();
    for(ctName in ctNames){
      others<-setdiff(ctNames, ctName);
      x<-setdiff( ctGenes[[ctName]], unlist(ctGenes[others]));
      specGenes[[ctName]]<-x;
    }
    ans<-specGenes;
  }
  else{
    ans<-ctGenes;
  }
  ans;
}

#' return a pattern for use in cn_testPattern (template matching)
#'
#' return a pattern for use in cn_testPattern (template matching)
#' @param sampR vector
#'
#' @return ans
cn_sampR_to_pattern<-function (sampR){
  d_ids<-unique(as.vector(sampR));
  nnnc<-length(sampR);
  ans<-matrix(nrow=length(d_ids), ncol=nnnc);
  for(i in seq(length(d_ids))){
    x<-rep(0,nnnc);
    x[which(sampR==d_ids[i])]<-1;
    ans[i,]<-x;
  }
  colnames(ans)<-as.vector(sampR);
  rownames(ans)<-d_ids;
  ans;
}

#' template matching
#'
#' test correlation between idealized expression pattern and target gene
#' @param pattern vector of pattern
#' @param expDat expression matrix
#'
#' @return data.frame(row.names=geneids, pval=pval, cval=cval,holm=holm);
#'
cn_testPattern<-function(pattern, expDat){
  pval<-vector();
  cval<-vector();
  geneids<-rownames(expDat);
  llfit<-ls.print(lsfit(pattern, t(expDat)), digits=25, print=FALSE);
  xxx<-matrix( unlist(llfit$coef), ncol=8,byrow=TRUE);
  ccorr<-xxx[,6];
  cval<- sqrt(as.numeric(llfit$summary[,2])) * sign(ccorr);
  pval<-as.numeric(xxx[,8]);

  #qval<-qvalue(pval)$qval;
  holm<-p.adjust(pval, method='holm');
  #data.frame(row.names=geneids, pval=pval, cval=cval, qval=qval, holm=holm);
  data.frame(row.names=geneids, pval=pval, cval=cval,holm=holm);
}

#' finds general and context dependent specifc genes
#'
#' finds general and context dependent specifc genes
#' @param expDat expression matrix
#' @param sampTab sample table
#' @param holm pvalue threshold for template matching
#' @param cval template matching threshold for overall CT specific expression
#' @param cvalGK template matching threshold for developmentally shared CT specific expression
#' @param dLevel "description1",
#' @param dLevelGK "description2"
#'
#' @return list of $matcher${cell_type}->{germ_layer}$context$general${cell_type}->gene vector etc
cn_specGenesAll<-function(expDat, sampTab,holm=1e-4,cval=0.5,cvalGK=0.75, dLevel="description1", dLevelGK=NULL,prune=FALSE){
  matcher<-list();
  general<-cn_findSpecGenes(expDat, sampTab, holm=holm, cval=cval, dLevel=dLevel,prune=prune);
  ctXs<-list()# one per germlayer
  if(!is.null(dLevelGK)){

    germLayers<-unique(as.vector(sampTab[,dLevelGK]));
    for(germlayer in germLayers){
      stTmp<-sampTab[sampTab[,dLevelGK]==germlayer,];
      expTmp<-expDat[,rownames(stTmp)];
      xxx<-cn_findSpecGenes(expTmp, stTmp, holm=holm, cval=cvalGK,dLevel=dLevel, prune=prune);
      cts<-names(xxx);
      for(ct in cts){
        matcher[[ct]]<-germlayer;
        # remove general ct-specific genes from this set
        a<-general[[ct]];
        b<-xxx[[ct]];
        ba<-setdiff(b, a);
        both<-union(a,b);
        xxx[[ct]]<-ba;
      }
      ctXs[[germlayer]]<-xxx;
    }
  }
  ctXs[['general']]<-general;
  list(context=ctXs, matcher=matcher);
}

#' extract sub-networks made up of CT genes;
#'
#' extract sub-networks made up of CT genes;
#' @param rawGRNs result of running cn_getRawGRN
#' @param specGenes result of running cn_specGenesAll
#'
#' @return list(geneLists=geneLists, graphLists=graphLists, tfTargets=tfTargets)
cn_specGRNs<-function(rawGRNs, specGenes){

  # should return a list of gene lists and igraphs
  geneLists<-list();
  graphLists<-list();

  groupNames<-names(specGenes[['context']][['general']]);

  big_graph<-rawGRNs[['graph']];

  matcher<-specGenes$matcher;

  allgenes<-V(big_graph)$name;

  for(ct in groupNames){
    cat(ct,"\n")
    if(!is.null(names(matcher))){
      gll<-matcher[[ct]];
      cat(ct," ",gll,"\n");
      mygenes<-union(specGenes[['context']][['general']][[ct]], specGenes[['context']][[gll]][[ct]]);
    }
    else{
      mygenes<-specGenes[['context']][['general']][[ct]]
    }

    geneLists[[ct]]<-intersect(allgenes, mygenes);
    graphLists[[ct]]<-induced.subgraph(big_graph, geneLists[[ct]]);
  }

  tfTargets<-cn_MakeTLs(graphLists);

  list(geneLists=geneLists, graphLists=graphLists, tfTargets=tfTargets);
}

#' get targets of tFs
#'
#' get targets of tFs
#' @param graphList a list of networks represented as iGraphs
#'
#' @return list of tf=>targets
cn_MakeTLs<-function(graphList){
  tfTargs<-list();
  nnames<-names(graphList);
  for(nname in nnames){
    tfTargs[[nname]]<-cn_get_targets(graphList[[nname]]);
  }
  tfTargs;
}

#' get targets of a tf
#'
#' get targets of a tf
#' @param aGraph an iGraph
#'
#' @return target list
cn_get_targets<-function(aGraph){
  targList<-list();
  regs<-V(aGraph)$label[V(aGraph)$type=='Regulator'];
  if(length(regs)>0){
    for(reg in regs){
      targList[[reg]]<-unique(sort(V(aGraph)$label[neighbors(aGraph, reg)]));
    }
  }
  targList;
}

#' GRN status
#'
#' Calculates the status of all GRNs in query samples as compared to training data for
#' @param expDat query expression matrix
#' @param subList of ct => genes
#' @param tVals tvals
#' @param classList classList
#' @param minVals minVals
#' @param classWeight class weight
#' @param exprWeight  expression weight
#' @return grn scores (not normalized)
cn_netScores<-function (expDat, genes, tVals, ctt, classList=NULL, classWeight=FALSE, exprWeight=TRUE,xmax=1e3){
  cat(ctt,"\n")
  aMat<-matrix(0, nrow=length(genes), ncol=ncol(expDat));
  rownames(aMat)<-genes;

  weights<-rep(1, length(genes));
  names(weights)<-genes;

  #otherCTs<-setdiff(names(tVals), ct)

  cat(dim(aMat),"\n")
  if(exprWeight){
    meanVect<-unlist(tVals[[ctt]][['mean']][genes]);
    weights<-(2**meanVect)/sum(2**meanVect);
  }

  if(classWeight){ #TODO modify this to fit the gene pairs
    #classImp<-classList[[ctt]]$importance[genes,1];

    classImp = weights
    for(gene in names(classList)) {
      if(gene %in% names(classImp)) {
        classImp[gene] = classList[gene] / 4 + 1
      }
    }

    ### 04-19-17
    ###classImp<-classImp/sum(classImp)
    weights<-weights*classImp;
  }

  for(gene in genes){
    ### cat("***",gene,"\n")
    ###zzs<-as.matrix(cn_rawScore(expDat[gene,], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]])[1,])

    zzs<-cn_rawScore(expDat[gene,], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]], xmax=xmax)
    aMat[gene,]<-zzs;

  }

  print(dim(aMat))
  xscores<-apply(aMat, 2, weighted.mean, w=weights);
  xscores;
}

#' Estimate gene expression dist in CTs
#'
#' Calculate mean and SD
#' @param expDat training data
#' @param sampTab, ### training sample table
#' @param dLevel="description1", ### column to define CTs
#' @param predictSD=FALSE ### whether to predict SD based on expression level
#'
#' @return tVals list of ct->mean->named vector of average gene expression, ->sd->named vector of gene standard deviation
cn_make_tVals<-function (expDat, sampTab, dLevel="description1", predictSD=FALSE){

  if(predictSD){
    ans<-cn_make_tVals_predict(expDat, sampTab, dLevel);
  }
  else{
    # Note: returns a list of dName->gene->mean, sd, where 'dName' is a ctt or lineage
    # make sure everything is lined up
    expDat<-expDat[,rownames(sampTab)];
    tVals<-list();
    dNames<-unique(as.vector(sampTab[,dLevel]));
    allGenes<-rownames(expDat);
    for(dName in dNames){
      #cat(dName,"\n");
      xx<-which(sampTab[,dLevel]==dName);
      sids<-rownames(sampTab[xx,]);
      xDat<-expDat[,sids];
      means<-apply(xDat, 1, mean);
      sds<-apply(xDat, 1, sd);
      tVals[[dName]][['mean']]<-as.list(means);
      tVals[[dName]][['sd']]<-as.list(sds);
    }
    ans<-tVals;
  }
  ans;
}

#' cn_make_tVals_predict
#'
#' predicts SD based on mean expression
#' @param expDat training data
#' @param sampTab training sample table
#' @param dLevel="description1" column to define CT
#'
#' @return tVals list of ct->mean->named vector of average gene expression, ->sd->named vector of gene standard deviation
cn_make_tVals_predict<-function(expDat, sampTab, dLevel="description1"){
  # Note: returns a list of dName->gene->mean, sd, where 'dName' is a ctt or lineage
  # make sure everything is lined up
  expDat<-expDat[,rownames(sampTab)];
  tVals<-list();
  dNames<-unique(as.vector(sampTab[,dLevel]));
  allGenes<-rownames(expDat);

  # make a model to predict SD given average expression level across all samples
  sdT<-apply(expDat, 1, sd);
  mT<-apply(expDat, 1, mean);
  myModel<-lm(sdT~mT);
  for(dName in dNames){
    xx<-which(sampTab[,dLevel]==dName);
    sids<-rownames(sampTab[xx,]);
    xDat<-expDat[,sids];
    means<-apply(xDat, 1, mean);
    sds<-predict(myModel, data.frame(mT=means));
    tVals[[dName]][['mean']]<-as.list(means);
    tVals[[dName]][['sd']]<-as.list(sds);
  }
  tVals;
}

#' min diff
#'
#' computes mean gene A in CT 1 - mean gene A in CT 2, where CT2 has the non CT1 max value. does this for genes
#' @param tVals tVals
#' @param genes vector of gene names
#' @param ct ct to compare to
#' @return vector of differences
minDif<-function(tVals, genes, ct){
  octs<-setdiff(names(tVals), ct)
  qq<-lapply(tVals[octs], "[[", "mean")
  ##tVals[[ct]][["mean"]][[gene]]#-max(unlist(lapply(qq, "[[", gene)))
  tmpMat<-matrix(unlist(lapply(qq, "[", genes)), nrow=length(genes))
  rownames(tmpMat)<-genes
  maxes<-apply(tmpMat, 1, max)
  unlist(tVals[[ct]][["mean"]][genes])-maxes
}

#' computes the raw score for a gene as xmax-abs(zscore).
#'
#' better values are higher.
#' @param vect a vector of gene expression values for multiple samples
#' @param mmean mean value in training data
#' @param ssd standard deviation in training data
#'
#' @return transformed (but not normalized) GRN score
#'
cn_rawScore<-function(vect, mmean, ssd, xmax=1e3){
  zcs<-zscore(vect, mmean, ssd);
  ### xmax<-1000; # arbitrary, and corrected for later, but want some high enough that it should not be exceeded
  xmax-abs(zcs);
}

#' GRN status
#'
#' Calculates the GRN status in query samples as compared to training data
#' @param expDat query expression matrix
#' @param subList of ct => genes
#' @param tVals list of ctt->list(means->genes, sds->genes)
#' @param classList class list
#' @param minVals  min vals
#' @param classWeight classweght
#' @param exprWeight  expression weight
#' @return GRN scores
#'
cn_score<-function(expDat, subList, tVals, classList=NULL, minVals=NULL, classWeight=FALSE, exprWeight=TRUE, xmax=1e3){
  #nSubnets<-sum(sapply(subList, length));
  nSubnets<-length(subList);
  ans<-matrix(0, nrow=nSubnets, ncol=ncol(expDat));
  ctts<-names(subList);
  rnames<-vector();
  rIndex<-1;
  for(ctt in ctts){
    cat(ctt,"\n");
    genes<-subList[[ctt]];
    # 06-06-16 -- added to allow for use of GRNs defined elsewhere
    genes<-intersect(genes, rownames(expDat));
    #    snNames<-names(subnets);
    #    rnames<-append(rnames, snNames);
    #    for(sName in snNames){
    ans[rIndex,]<-cn_netScores(expDat, genes, tVals=tVals, ctt=ctt,classList=classList, classWeight=classWeight,exprWeight=exprWeight, xmax=xmax);
    rnames<-append(rnames, ctt);
    rIndex<-rIndex+1;
    #   }
  }
  rownames(ans)<-rnames;
  colnames(ans)<-colnames(expDat);
  if(!is.null(minVals)){
    minVals<-minVals[rownames(ans)];
    ans<-ans-minVals;
  }
  ans;
}

#' Normalize grn status as compared to training data
#'
#' Divide the query scores by the mean values in the training data.
#' @param ctrlScores a list of subnet->mean value, all subnets
#' @param queryScores a matrix, rownames = subnet names, all subnets
#' @param subNets a vector of subnets names to use
#'
#' @return normalized grn status matrix
#'
cn_normalizeScores<-function(ctrlScores, queryScores, subNets){

  ans<-matrix(0, nrow=length(subNets), ncol=ncol(queryScores));
  rownames(ans)<-subNets
  #subNets<-rownames(queryScores);
  for(subNet in subNets){
    ### cat(subNet,"\n")
    ans[subNet,]<- queryScores[subNet,] / ctrlScores[[subNet]];
  }
  colnames(ans)<-colnames(queryScores);
  ans;
}


#' Figure out normalization factors for GRNs, and norm training data
#'
#' Exactly that.
#' @param expTrain expression matrix
#' @param stTrain  sample table
#' @param subNets named list of genes, one list per CTT, tct=>gene vector
#' @param classList list of classifiers
#' @param dLevel column name to group on
#' @param tVals seful when debugging
#' @param classWeight weight GRN status by importance of gene to classifier
#' @param exprWeight weight GRN status by expression level of gene?
#' @param sidCol sample id colname
#'
#' @return list of trainingScores, normVals, raw_scores, minVals, tVals=tVals
#' @export
cn_trainNorm<-function (expTrain, stTrain, subNets, classList = NULL,  dLevel = "description1", tVals=NULL, classWeight=FALSE, exprWeight=TRUE, sidCol='sample_id', xmax=1e3, predSD=FALSE){

  if(is.null(tVals)){
    tVals<-cn_make_tVals(expTrain, stTrain, dLevel, predictSD=predSD)
  }

  ctts<-as.vector(unique(stTrain[,dLevel]));
  scoreList<-list();
  normList<-list(); # a list of ctt->subnet->mean value
  minVect<-vector(); # a list of ctt->subnet->min value, used to shift raw grn est scores

  cat("calculating GRN scores on training data ...\n");
  tmpScores<-cn_score(expTrain, subNets, tVals, classList, minVals=NULL, classWeight=classWeight, exprWeight=exprWeight, xmax=xmax)


  minVect<-apply(tmpScores, 1, min);
  names(minVect)<-rownames(tmpScores);

  # shift the raw scores so that min=0;
  tmpScores<-tmpScores - minVect;
  cat("norm factors\n");
  for(ctt in ctts){
    # determine nomalization factors
    ##snets<-names(subNets[[ctt]]);
    snets<-ctt;

    scoreDF<-cn_extract_SN_DF(tmpScores, stTrain, dLevel, snets, sidCol=sidCol);
    scoreDF<-cn_reduceMatLarge(scoreDF, "score", "description", "subNet");
    xdf<-scoreDF[which(scoreDF$grp_name==ctt),];
    tmpSNS<-as.list(xdf$mean);
    names(tmpSNS)<-xdf$subNet;
    normList[names(tmpSNS)]<-tmpSNS;
  }

  # normalize training scores
  nScores<-cn_normalizeScores(normList, tmpScores, rownames(tmpScores));

  scoreDF<-cn_extract_SN_DF(nScores, stTrain, dLevel, sidCol=sidCol);

  scoreDF<-cn_reduceMatLarge(scoreDF, "score", "description", "subNet");

  list(trainingScores=scoreDF,
       normVals=normList,
       raw_scores=tmpScores,
       minVals=minVect,
       tVals=tVals);
}

zscore<-function(x,meanVal,sdVal){
  (x-meanVal)/sdVal;
  ### zscore
}


#' returns a DF of: sample_id, description, ctt, subnet_name, score
#'
#' returns a DF of: sample_id, description, ctt, subnet_name, score
#' @param scores a matrix of subNet scores
#' @param sampTab sample table
#' @param dLevel column name of sampTab to group values by
#' @param rnames rownames to extract
#' @param sidCol sample identifier column name
#'
#' @return returns a DF of: sample_id, description, ctt, subnet_name, score
cn_extract_SN_DF<-function(scores, sampTab, dLevel, rnames=NULL, sidCol="sample_id"){

  if(is.null(rnames)){
    rnames<-rownames(scores);
    #cat("GOT NULL\n");
  }

  tss<-scores[rnames,];
  if(length(rnames)==1){
    tss<-t(as.matrix(scores[rnames,]));
    rownames(tss)<-rnames;
    #  cat(dim(tss),"\n")
  }

  nSamples<-ncol(tss);
  stTmp<-sampTab[colnames(tss),]; ####
  snNames<-rownames(tss);
  num_subnets<-length(snNames);
  snNames<-unlist(lapply(snNames, rep, times=nSamples));
  sample_ids<-rep(as.vector(stTmp[,sidCol]), num_subnets);
  descriptions<-rep(as.vector(stTmp[,dLevel]), num_subnets);
  # myCtts<-rep(ctt, length(snNames));
  scores<-as.vector(t(tss));
  data.frame(sample_id=sample_ids,
             description=descriptions,
             #         ctt=myCtts,
             subNet = snNames,
             score=scores);
  ### data.frame
}

cn_reduceMatLarge<-function (datFrame, valCol="score", cName="description", iterOver="subNet"){

  iterOvers<-unique(as.vector(datFrame[,iterOver]));
  ans<-data.frame();
  for(io in iterOvers){
    #  cat(io,"\n");
    xi<-which(datFrame[,iterOver]==io);
    dfTmp<-datFrame[xi,];
    x<- utils_reduceMat(dfTmp,valCol=valCol,cName=cName);
    x<-cbind(x, subNet=rep(io, nrow(x)));
    ans<-rbind(ans, x);
  }
  ans;
  ### ans
}

### reduce the data.matrix values by averaging and getting st dvs
utils_reduceMat<-function(datFrame, valCol, cName='ann'){

  mids<-unique(as.vector(datFrame[,cName]));
  means<-rep(0, length(mids));
  sds<-rep(1, length(mids));
  indexI<-rep(0, length(mids)); # index to extract other columns from original data.frame

  for(i in seq(length(mids))){
    mid<-mids[i]
    xi<-which(datFrame[,cName]==mid);
    tmpDat<-datFrame[xi,];
    means[i]<-mean(tmpDat[,valCol]);
    sds[i]<-sd(tmpDat[,valCol]);
  }
  data.frame(grp_name=mids, mean=means, stdev=sds);
  ### df of grp_name, mean, sd
}

###########################
# helper functions
#' @title logRank gene expression
#' @description rank the gene expression values and log10 the ranks
#' @param expTrain the expression matrix
#' @return the log10 rank of genes
logRank <- function(expTrain) {
  expTrain = apply(expTrain, FUN = rank, MARGIN = 2)

  return(log10(expTrain))
}

processImportance <- function(classifier) {
  genePairImportance = classifier$importance

  genes = unique(unlist(strsplit(x = rownames(classifier$importance), split = "_")))

  geneImportance = rep(0, length(genes))

  names(geneImportance) = genes
  genePairImportance = genePairImportance[, 1]

  for(genePair in names(genePairImportance)) {
    temp_importance = genePairImportance[genePair]
    genePairSplit = unlist(strsplit(x = genePair, split = "_"))

    if(temp_importance > geneImportance[genePairSplit[1]]) {
      geneImportance[genePairSplit[1]] = temp_importance
    }
    #if(temp_importance > geneImportance[genePairSplit[2]]) {
    #    geneImportance[genePairSplit[2]] = temp_importance
    #}
  }

  return(geneImportance)
}
