#' @title
#' Heatmap of the classification result
#' @description
#' This function generates a heatmap of the classification result for visualization
#'
#' @param classMat classification matrix generated from \code{\link{rf_classPredict}}
#' @param grps
#' @param isBig TRUE if this is a big heatmap
#' @param cRow Boolean values determining if rows should be clustered
#' @param cCol Boolean values determining if columns should be clustered
#' @param fontsize_row
#' @param scale
#'
#' @return nothing
#'
#' @examples
#' ccn_HmClass(cnRes, isBig=TRUE)
#'
#' @export
ccn_hmClass<-function(classMat, grps=NULL, isBig=FALSE, cRow=FALSE, cCol=FALSE, fontsize_row=4, scale=FALSE){

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
