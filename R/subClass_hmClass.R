#' @title
#' Detailed Complex heatmap for subclass classification
#' @description
#' Generates a complex heatmap for visualizing subclass classification
#'
#' @param subClassMat classification matrix of cancer sub-category
#' @param broadClassMAt classification matrix of cancer broad category
#' @param cScoreCat the broad category that users want to display the classification score
#' @return classification heatmap
#'
#' @import RColorBrewer
#'
#' @export
subClass_hmClass <- function(subClassMat, broadClassMat, cScoreCat) {

  broadClassMat = broadClassMat[, colnames(subClassMat)]
  list_max = apply(broadClassMat, FUN = max, MARGIN = 2)

  # this gets the number of labels a model has
  model_label = c()
  for(colName in colnames(broadClassMat)) {
    model_label = c(model_label, names(which.max(broadClassMat[, colName])))

  }

  cools<-colorRampPalette(c("black", "limegreen", "yellow"))( 100 )
  xcol <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(unique(model_label)))

  names(xcol) = unique(model_label)

  column_ha = ComplexHeatmap::HeatmapAnnotation("Broad Classification" = model_label,
                                                "Broad C-Scores" = anno_barplot(broadClassMat[grep(cScoreCat, rownames(broadClassMat)), ],
                                                                                axis = TRUE, axis_side = 'right',
                                                                                ylim = c(0, 1), bar_width = 0.9,
                                                                                gp = gpar(fill = "#66ADE5")),

                                                height = unit(2, "cm"),
                                                col = list("Classification Results" = xcol),
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
