# (C) Patrick Cahan 2012-2018
# CancerCellNet


#' heatmap genes and groups
#'
#' heatmap genes and groups
#'
#' @param expDat expDat
#' @param genes genes
#' @param grps vector of cellnames -> grp label
#' @param maxPerGrp 100
#' @param cRow =FALSE,
#' @param cCol =FALSE,
#' @param limits =c(0,10),
#' @param toScale =FALSE,
#' @param fontsize_row =4
#'
#' @return pheatmap
#'
#' @export
hm_gpa_sel<-function(
  expDat,
  genes,
  grps, ## vector of cellnames -> grp label
  maxPerGrp=100,
  cRow=FALSE,
  cCol=FALSE,
  limits=c(0,10),
  toScale=FALSE,
  fontsize_row=4,
  reOrderCells=FALSE){


  allgenes<-rownames(expDat)
  missingGenes<-setdiff(genes, allgenes)
  if(length(missingGenes)>0){
    cat("Missing genes: ", paste0(missingGenes, collapse=","), "\n")
    genes<-intersect(genes, allgenes)
  }

  value<-expDat[genes,]
  if(toScale){
      value <- t(scale(t(value)))
    }

  value[value < limits[1]] <- limits[1]
  value[value > limits[2]] <- limits[2]

  groupNames<-unique(grps)
  if(reOrderCells){
    grps<-grps[order(grps)]
    groupNames<-sort(unique(grps))
  }

  cells<-names(grps)

##
 ## groupNames<-myGrpSort(grps)
##

  cells2<-vector()
  for(groupName in groupNames){
    xi<-which(grps==groupName)
    if(length(xi)>maxPerGrp){
      tmpCells<-sample(cells[xi], maxPerGrp)
    }
    else{
      tmpCells<-cells[xi]
    }
    cells2<-append(cells2, tmpCells)
  }
  value<-value[,cells2]

  xcol <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(groupNames))
    names(xcol) <- groupNames
    anno_colors <- list(group = xcol)

    xx<-data.frame(group=as.factor(grps))
    rownames(xx)<-cells

   val_col <- colorRampPalette(rev(brewer.pal(n = 12,name = "Spectral")))(25)
   #val_col <- colorRampPalette(brewer.pal(n = 12,name = "Spectral"))(100)

  pheatmap(value, cluster_rows = cRow, cluster_cols = cCol, color=val_col,
        show_colnames = FALSE, annotation_names_row = FALSE,
##        annotation_col = annTab,
        annotation_col = xx,
        annotation_names_col = FALSE, annotation_colors = anno_colors, fontsize_row=fontsize_row)
}



#' simple heatmap of the classification result
#'
#' Heatmap of the classification result.
#' @param cnRes returned from cn_sapply
#' @param isBig is this a big heatmap? TRUE or FALSE
#'
#' @return nothing
#'
#' @examples
#' cn_HmClass(cnRes, isBig=TRUE)
#'
#' @export
cn_HmClass<-function
(classRes,
 isBig=FALSE
){

  cools<-colorRampPalette(c("black", "limegreen", "yellow"))( 100 )
  bcol<-'white';
  if(isBig){
    bcol<-NA;
  }
  pheatmap(classRes,
    col=cools,
    breaks=seq(from=0, to=1, length.out=100),
    border_color=bcol,
    cluster_rows = FALSE,
    cluster_cols = FALSE)
  # classification heatmap
}


#' @title
#' Heatmap of the classification result
#' @description
#' This function generates a heatmap of the classification result for visualization
#'
#' @param classMat classification matrix generated from \code{\link{rf_classPredict}}
#' @param isBig is this a big heatmap? TRUE or FALSE
#' @param cluster_cols cluster_cols
#'
#' @return nothing
#'
#' @examples
#' ccn_HmClass(cnRes, isBig=TRUE)
#'
#' @export
ccn_hmClass<-function(
  classMat,
  grps=NULL, ## vector of cellnames -> grp label
  isBig=FALSE,
  cRow=FALSE,
  cCol=FALSE,
  fontsize_row=4,
  scale=FALSE
){

  cools<-colorRampPalette(c("black", "limegreen", "yellow"))( 100 )
  bcol<-'white';
  if(isBig){
    bcol<-NA;
  }

  if(is.null(grps)){
  	cn_HmClass(classMat,isBig=isBig)
  }
  else{

	grps<-grps[order(grps)]
  	cells<-names(grps)
  	groupNames<-sort(unique(grps))
	xcol <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(groupNames))
  	names(xcol) <- groupNames
  	anno_colors <- list(group = xcol)

 	xx<-data.frame(group=as.factor(grps))
  	rownames(xx)<-cells


  if(scale){
    mymin<-min(classMat)
    mymax<-max(classMat)
  }
  else{
    mymin<-0
    mymax<-1
  }

 	 pheatmap(classMat, col=cools, breaks=seq(from=mymin, to=mymax, length.out=100), cluster_rows = cRow, cluster_cols = cCol,
        show_colnames = FALSE, annotation_names_row = FALSE,
##        annotation_col = annTab,
        annotation_col = xx,
        annotation_names_col = FALSE, annotation_colors = anno_colors, fontsize_row=fontsize_row)
 	}

}


#' @title
#' add "random" profiles to the meta sample table to assist with plotting
#'
#' @description
#' This will add the random profiles generated from \code{\link{rf_classPredict}} to the meta sample table
#'
#' @param classRes the classification matrix generated from \code{\link{rf_classPredict}}
#' @param sampTab the meta sample table
#' @param desc the "column name" of the column containing the sample descriptions in the meta sample table
#' @param id the column name of the column containing the sample names in the meta sample table
#'
#' @return meta sample table containing the random profiles
#'
#' @examples
#' stValRand<-addRandToSampTab(classRes_val, stVal, "description2", "sample_name")
#' @export
addRandToSampTab<-function(classRes, sampTab, desc, id="cell_name") {
	cNames<-colnames(classRes)
	snames<-rownames(sampTab)

	rnames<-setdiff(cNames, snames)

	cat("number of random samples: ",length(rnames), "\n")

	stNew<-data.frame(rid=rnames, rdesc=rep("rand", length(rnames)))
	stTop<-sampTab[,c(id, desc)]
	colnames(stNew)<-c(id, desc)

	ans<-rbind(stTop, stNew)
	rownames(ans)<-colnames(classRes)

	#return
	ans
}




#' Plot results of cn_classAssess
#'
#' Plot one precision recall curve per CT
#' @param assessed result of runnung cn_classAssess
#'
#' @return ggplot pbject
#'
#' @examples
#' testAssTues<-cn_splitMakeAssess(stTrain, expTrain, ctGRNs, prop=.5)
#' plot_class_PRs(testAssTues$ROCs)
#'
#' @export
plot_class_PRs<-function
(assessed
  ){
  ctts<-names(assessed)
  df<-data.frame()
  for(ctt in ctts){
    tmp<-assessed[[ctt]]
    tmp<-cbind(tmp, ctype=ctt)
    df<-rbind(df, tmp)
  }

  ggplot(data=df, aes(x=Sens, y=Prec)) + geom_point(size = .5, alpha=.5) +  geom_path(size=.5, alpha=.75) +
  theme_bw() + xlab("Recall") + ylab("Precision") + facet_wrap( ~ ctype, ncol=4) +
  theme(axis.text = element_text(size=5)) + ggtitle("Classifier performance")
}








