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
