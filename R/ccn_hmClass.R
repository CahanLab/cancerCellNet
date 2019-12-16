#' @title
#' Heatmap of the classification result
#' @description
#' This function generates a heatmap of the classification result for visualization
#'
#' @param classMat classification matrix generated from \code{\link{rf_classPredict}}
#' @param grps a vector that maps samples to a group
#' @param isBig TRUE if this is a big heatmap
#' @param cRow TRUE if rows should be clustered
#' @param cCol TRUE if columns should be clustered
#' @param fontsize_row the font size of the row
#' @param fontsize_col the font size of the columns
#' @param main the title of the heatmap
#' @param scale FALSE if the highest value is 1 and the lowest is 0
#' @param customAnnoColor a named vector with colors and named with group names
#' @return classification heatmap
#'
#' @examples
#' ccn_HmClass(cnRes, isBig=TRUE)
#'
#' @importFrom pheatmap pheatmap
#'
#' @export
ccn_hmClass<-function(classMat, grps=NULL, isBig=FALSE, cRow=FALSE, cCol=FALSE, fontsize_row=4, fontsize_col=4, main=NA, scale=FALSE, customAnnoColor = NULL, ...){

  cools<-colorRampPalette(c("black", "limegreen", "yellow"))( 100 )
  bcol<-'white';
  if(isBig){
    bcol<-NA;
  }

  # if no groups specified, use simple heatmap
  if(is.null(grps)){
    cn_HmClass(classMat,isBig=isBig, cRow=cRow, cCol=cCol, fontsize_row=fontsize_row, fontsize_col=fontsize_col, main=main)
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

    # if there is a custom color pallete
    if(is.null(customAnnoColor) == FALSE) {

      if(!all(groupNames %in% names(customAnnoColor))) {
        stop(paste0("Not all group name has a color in custom color vetor.", "\n"))

      }
      customAnnoColor <- customAnnoColor[groupNames]
      anno_colors <- list(group = customAnnoColor)

    }

    pheatmap::pheatmap(classMat,
             col=cools,
             breaks=seq(from=mymin, to=mymax, length.out=100),
             cluster_rows = cRow,
             cluster_cols = cCol,
             show_colnames = FALSE,
             annotation_names_row = FALSE,
             ##        annotation_col = annTab,
             annotation_col = xx,
             main = main,
             annotation_names_col = FALSE, annotation_colors = anno_colors, fontsize_row=fontsize_row, fontsize_col=fontsize_col,
             ...)
  }

}

#' @title
#' Simple heatmap of the classification result
#'
#' @description
#' Heatmap of the classification result that does not have an explicit group.
#'
#' @param classMat classification matrix generated from \code{\link{rf_classPredict}}
#' @param isBig TRUE if this is a big heatmap
#' @param cRow TRUE if rows should be clustered
#' @param cCol TRUE if columns should be clustered
#' @param fontsize_row the font size of the row labels
#' @param fontsize_col the font size of the column labels
#' @param main the title of the heatmap
#' @param scale FALSE if the highest value is 1 and the lowest is 0
#'
#' @return simpler classification heatmap without groups colored
cn_HmClass<-function(classRes, isBig=FALSE, cRow=FALSE, cCol=FALSE, fontsize_row=4, fontsize_col=4, main=NA, scale=FALSE){
  cools<-colorRampPalette(c("black", "limegreen", "yellow"))( 100 )
  bcol<-'white';

  if(isBig){
    bcol<-NA;
  }

  if(scale){
  mymin<-min(classMat)
  mymax<-max(classMat)
  }else{
    mymin<-0
    mymax<-1
  }

  pheatmap::pheatmap(classRes,
           col=cools,
           breaks=seq(from=mymin, to=mymax, length.out=100),
           border_color=bcol,
           cluster_rows = cRow,
           cluster_cols = cCol,
           main=main,
           fontsize_row = fontsize_row,
           fontsize_col = fontsize_col)
  # classification heatmap
}
