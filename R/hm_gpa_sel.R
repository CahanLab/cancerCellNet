#' @title
#' Heatmap of the genes and groups
#'
#' @description
#' Generate a heatmap with classification worthy genes and cancer groups
#'
#' @param expDat normalized expression matrix from \code{\link{trans_prop}}
#' @param genes classification worth genes from \code{\link{findClassyGenes}}
#' @param grps a vector that maps sample names to a cancer category
#' @param maxPerGrp the number max samples per group
#' @param cRow TRUE if the rows should be clustered
#' @param cCol TRUE if the columns should be clustered
#' @param limits a vector of size 2 indicating the lowest and highest range for expression level
#' @param toScale TRUE if the data should be scaled
#' @param fontsize_row the font size for the row labels
#'
#' @return a heatmap of genes and their groups
#'
#' @export
hm_gpa_sel<-function(expDat, genes, grps, maxPerGrp=100, cRow=FALSE, cCol=FALSE, limits=c(0,10), toScale=FALSE, fontsize_row=4, reOrderCells=FALSE){


  allgenes<-rownames(expDat)
  missingGenes<-setdiff(genes, allgenes) # find the genes that are not classification worthy
  if(length(missingGenes)>0){
    cat("Missing genes: ", paste0(missingGenes, collapse=","), "\n")
    genes<-intersect(genes, allgenes) #this line might be redundent since all the cgenes are gathered from ths
  }

  value<-expDat[genes,] #select the dataframe with cgenes

  if(toScale){
    value <- t(scale(t(value))) #scales the dataframe
  }

  value[value < limits[1]] <- limits[1] # ensures 0 is the smallest
  value[value > limits[2]] <- limits[2] # ensures 10 is the highest

  groupNames<-unique(grps)

  if(reOrderCells){
    grps<-grps[order(grps)]
    groupNames<-sort(unique(grps))
  } # alphabetical ordering

  cells<-names(grps)

  ##
  ## groupNames<-myGrpSort(grps)
  ##

  cells2<-vector()
  for(groupName in groupNames){
    xi<-which(grps==groupName) #select samples that are in a certain group

    if(length(xi)>maxPerGrp){
      tmpCells<-sample(cells[xi], maxPerGrp) #if over the maximum number of samples per group, then sample that amount
    }
    else{
      tmpCells<-cells[xi] #if not over the maximum number of samples per group, then use all the samples available
    }
    cells2<-append(cells2, tmpCells) # create a vector with all the samples selected for plotting
  }
  value<-value[,cells2] # select the samples that are going to be used for plotting


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
