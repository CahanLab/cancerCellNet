#' @title
#' Plot Specific Gene Expression
#' @description
#' Plots a specific gene and its expression levels across the sampels
#'
#' @param GeneCompareTab the gene comparsion matrix generated from \code{\link{makeGeneCompareTab}}
#' @param sdGeneCat the standard deviation matrix generated from \code{\link{sdGeneCat}}
#' @param gene the name of the specific gene of interest
#' @param text_size the size of the ticks
#' @param title_size the size of the labels
#'
#' @return expression bar plot
#'
#' @import ggplot2
#'
#' @export
barPlotGene<-function(GeneCompareTab, sdGeneCat = NULL, gene, text_size = 10, title_size = 12) {
   # set up an error trap later

  if (is.null(sdGeneCat) == FALSE) {
    sampGeneCompareTab = GeneCompareTab[gene, ]
    sampleSDGeneCat = sdGeneCat[gene, ]

    plotMatrix = matrix(ncol = 3, nrow = length(sampGeneCompareTab))
    rownames(plotMatrix) = names(sampGeneCompareTab)
    colnames(plotMatrix) = c("Samples", "Expression_Level", "Expression_SD")

    plotMatrix[, "Samples"] = rownames(plotMatrix)
    plotMatrix[, "Expression_Level"] = sampGeneCompareTab

    for(sample in plotMatrix[, "Samples"]) {
      temp = gsub(".{3}$", "", sample)
      SDindex = grep(temp, names(sampleSDGeneCat))

      if(length(SDindex) == 1) {
        plotMatrix[sample, "Expression_SD"] = sampleSDGeneCat[SDindex]
      }
      else {
        plotMatrix[sample, "Expression_SD"] = 0
      }

    }

  }
  else {
    sampGeneCompareTab = GeneCompareTab[gene, ]
    sampleSDGeneCat = sdGeneCat[gene, ]

    plotMatrix = matrix(ncol = 3, nrow = length(sampGeneCompareTab))
    rownames(plotMatrix) = names(sampGeneCompareTab)
    colnames(plotMatrix) = c("Samples", "Expression_Level", "Expression_SD")

    plotMatrix[, "Samples"] = rownames(plotMatrix)
    plotMatrix[, "Expression_Level"] = sampGeneCompareTab
    plotMatrix[, "Expression_SD"] = 0

  }


  plotMatrix = as.data.frame(plotMatrix)

  plotMatrix[, 2] = as.numeric(as.character(plotMatrix[, 2]))
  plotMatrix[, 3] = as.numeric(as.character(plotMatrix[, 3]))

  p = ggplot2::ggplot(data=plotMatrix, aes(x=Samples, y=Expression_Level)) +
    geom_bar(stat="identity", fill = "steel blue") + theme_minimal() +
    geom_errorbar(aes(ymin=Expression_Level-Expression_SD, ymax=Expression_Level+Expression_SD), width=.2,position=position_dodge(.9))
  p = p + theme(axis.text.x.bottom = element_text(angle = 30, hjust = 1), axis.text=element_text(size = text_size), axis.title = element_text(size = title_size)) +
    ylab(label = "Normalized Expression Level")

   #return
  p
}
