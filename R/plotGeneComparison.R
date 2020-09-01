#' @title
#' Gene expression plotting
#' @description
#' Plot gene expressions for visualization and comparison
#'
#' @param expDat comparison expression matrix from \code{\link{makeGeneCompareTab}}
#' @param fontsize_row row font size
#' @param cRows TRUE if cluster the rows
#' @param cCols TRUE if cluster the columns
#' @param grps the groups for the heatmap
#' @param ... pass parameters to pheamap
#' @return a heatmap of gene expressions for the comparsion expression matrix
#' @importFrom RColorBrewer brewer.pal
#' @export
plotGeneComparison<-function(expDat, fontsize_row = 6, grps = NULL, cRows = FALSE, cCols = FALSE, toScale = FALSE, ...){

  if(is.null(grps) == TRUE) {
    grps = colnames(expDat)
    names(grps) = colnames(expDat)
  }

  genes = rownames(expDat)

  value = expDat[genes,] #select the matrix with cgenes

  if(toScale){
    value = t(scale(t(value))) #scales the matrix
  }

  limits = c(0,10)

  value[value < limits[1]] = limits[1] # ensures 0 is the smallest
  value[value > limits[2]] = limits[2] # ensures 10 is the highest

  groupNames = unique(grps) # gather the names of the groups

  cells = names(grps)

  cells2 = vector()
  for(groupName in groupNames){
    xi = which(grps==groupName) #select samples that are in a certain group

    tmpCells = cells[xi] #if not over the maximum number of samples per group, then use all the samples available

    cells2 = append(cells2, tmpCells) # create a vector with all the samples selected for plotting
  }
  value = value[,cells2] # select the samples that are going to be used for plotting

  xcol = colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(groupNames))
  names(xcol) = groupNames
  anno_colors = list(group = xcol)

  xx = data.frame(group=as.factor(grps))
  rownames(xx) = cells

  val_col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 12,name = "Spectral")))(25)

  return(pheatmap(value, color=val_col, cluster_row = cRows, cluster_cols = cCols,
           show_colnames = FALSE, annotation_names_row = FALSE,
           annotation_col = xx,
           annotation_names_col = FALSE,
           annotation_colors = anno_colors,
           fontsize_row=fontsize_row, ...))
}
