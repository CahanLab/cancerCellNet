#' @title
#' Detailed Complex heatmap for subclass classification
#' @description
#' Generates a complex heatmap for visualizing subclass classification
#'
#' @param subClassMat classification matrix of cancer sub-category
#' @param broadClassMAt classification matrix of cancer broad category
#' @param cScoreCat the broad category that users want to display the classification score
#' @param thresHold_list list of threshold for broad class classification
#' @return classification heatmap
#'
#' @import RColorBrewer ComplexHeatmap circlize grid
#'
#' @export
subClass_hmClass <- function(subClassMat, broadClassMat, subCat, cScoreCat, thresHold_list = NULL) {

  broadClassMat = broadClassMat[, colnames(subClassMat)]
  list_max = apply(broadClassMat, FUN = max, MARGIN = 2)

  model_label = c()
  if(is.null(thresHold_list) == TRUE) {
    for(colName in colnames(broadClassMat)) {
      model_label = c(model_label, names(which.max(broadClassMat[, colName])))

    }
  }
  else {
    for(colName in colnames(broadClassMat)) {
      model_label = c(model_label, findCategory(broadClassMat[, colName], thresHold_list))

    }
  }

  cools<-colorRampPalette(c("black", "limegreen", "yellow"))( 100 )

  cools = circlize::colorRamp2(seq(0, 1,length.out = 100), cools)
  xcol <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(unique(model_label)))
  subClassColor = colorRampPalette(rev(brewer.pal(n = 12,name = "Set3")))(length(unique(subCat)))

  names(xcol) = unique(model_label)

  column_ha = ComplexHeatmap::HeatmapAnnotation("Broad Classification" = model_label,
                                                "Broad C-Scores" = anno_barplot(broadClassMat[grep(cScoreCat, rownames(broadClassMat)), ],
                                                                                       axis = TRUE, axis_side = 'right',
                                                                                       ylim = c(0, 1), bar_width = 0.9,
                                                                                       gp = gpar(fill = "#66ADE5")),
                                                "Subclass Classification" = subCat,


                                                height = unit(4, "cm"),
                                                col = list("Broad Classification" = xcol),
                                                show_annotation_name = TRUE,
                                                annotation_name_offset = unit(9, "mm"), gap = unit(2, "mm"),
                                                annotation_name_gp = par(ps = 0.01, font = 4, cex=0.8))

  p = ComplexHeatmap::Heatmap(subClassMat, name = "SubClassification Score",
                              heatmap_legend_param = list(
                                title = "", at = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                legend_height = unit(6, "cm")
                              ),
                              top_annotation = column_ha, col = cools,
                              cluster_rows = FALSE,
                              cluster_columns = FALSE )
  #return
  p
}

#' @title
#' Find the broad category
#' @description
#' Find the broad cateogry based on a threshold
#'
#' @param ThisColumn the individual sample in a matrix
#' @param thresHold_list the list of threshold for each category
#' @return classification heatmap
findCategory <- function(ThisColumn, thresHold_list) {
  returnList = ""
  for(thisnames in names(thresHold_list)) {
    if(ThisColumn[thisnames] >= thresHold_list[thisnames]) {
      if (returnList == "") {
        returnList = thisnames
      }
      else {
        returnList = paste(returnList, sep = ", ", thisnames)
      }
    }
  }

  if (returnList == "") {
    returnList = "Not Classified"
  }

  # return
  returnList
}
