#' @describeIn getPCA Given a set of genes, return a ggplot of expression
#' values.
#'
#' @param visual Type of visualization (PCA or tSNE). Default: "PCA"
#' @param x x coordinate for PCA
#' @param y y coordinate for PCA
#'
#' @return plotBiomarker(): A Biomarker plot
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' plotBiomarker(mouseBrainSubsetSCE, gene="C1qa", shape="level1class")
#'
plotBiomarker <- function(inSCE, gene, binary="Binary", visual="PCA",
                          shape="No Shape", x="PC1", y="PC2",
                          useAssay="counts", reducedDimName="PCA"){
  if (shape == "No Shape"){
    shape <- NULL
  }
  if (visual == "PCA"){
    if (is.null(SingleCellExperiment::reducedDim(inSCE, reducedDimName))) {
      inSCE <- getPCA(inSCE, useAssay = useAssay,
                      reducedDimName = reducedDimName)
    }
    axisDf <- data.frame(SingleCellExperiment::reducedDim(inSCE,
                                                           reducedDimName))
    variances <- NULL
    if (class(inSCE) == "SCtkExperiment"){
      variances <- pcaVariances(inSCE)
    }
  }
  if (visual == "tSNE"){
    if (is.null(SingleCellExperiment::reducedDim(inSCE, reducedDimName))) {
      inSCE <- getTSNE(inSCE, useAssay = useAssay,
                       reducedDimName = reducedDimName)
    }
    axisDf <- data.frame(SingleCellExperiment::reducedDim(inSCE,
                                                           reducedDimName))
  }
  if (length(gene) > 9) {
    gene <- gene[seq_len(9)]
  }
  for (i in seq_along(gene)){
    bioDf <- getBiomarker(inSCE = inSCE, gene = gene[i], binary = binary,
                          useAssay = useAssay)
    l <- axisDf
    if (!is.null(shape)){
      l$shape <- factor(SingleCellExperiment::colData(inSCE)[, shape])
    }
    geneName <- colnames(bioDf)[2]
    colnames(bioDf)[2] <- "expression"
    l$Sample <- as.character(bioDf$sample)
    l$expression <- bioDf$expression
    c <- SummarizedExperiment::assay(inSCE, useAssay)[c(geneName), ]
    percent <- round(100 * sum(c > 0) / length(c), 2)
    if (visual == "PCA"){
      if (binary == "Binary"){
        l$expression <- ifelse(l$expression, "Yes", "No")
        g <- ggplot2::ggplot(l, ggplot2::aes_string(x, y, label = "Sample",
                                                    color = "expression")) +
          ggplot2::geom_point() +
          ggplot2::scale_color_manual(limits = c("Yes", "No"),
                                      values = c("Blue", "Grey")) +
          ggplot2::labs(color = "Expression")
      }
      else if (binary == "Continuous"){
        if (min(round(l$expression, 6)) == max(round(l$expression, 6))){
          g <- ggplot2::ggplot(l, ggplot2::aes_string(x, y, label = "Sample")) +
            ggplot2::geom_point(color = "grey")
        } else{
          g <- ggplot2::ggplot(l, ggplot2::aes_string(x, y, label = "Sample",
                                                      color = "expression")) +
            ggplot2::scale_colour_gradient(limits = c(min(l$expression),
                                                      max(l$expression)),
                                           low = "grey", high = "blue") +
            ggplot2::geom_point()
        }
        g <- g + ggplot2::labs(color = "Expression")
      }
      g <- g +
        ggplot2::ggtitle(paste(geneName, " - ", percent, "%", " cells",
                               sep = "")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      if (is.null(variances)){
        g <- g + ggplot2::labs(x = x, y = y)
      } else {
        g <- g + ggplot2::labs(
          x = paste0(x, " ", toString(round(variances[x, ] * 100, 2)), "%"),
          y = paste0(y, " ", toString(round(variances[y, ] * 100, 2)), "%"))
      }
    } else if (visual == "tSNE"){
      if (binary == "Binary"){
        l$expression <- ifelse(l$expression, "Yes", "No")
        g <- ggplot2::ggplot(l, ggplot2::aes_string("X1", "X2",
                                                    label = "Sample",
                                                    color = "expression")) +
          ggplot2::geom_point() +
          ggplot2::scale_color_manual(limits = c("Yes", "No"),
                                      values = c("blue", "grey")) +
          ggplot2::labs(color = "Expression")
      }
      else if (binary == "Continuous"){
        if (min(round(l$expression, 6)) == max(round(l$expression, 6))) {
          g <- ggplot2::ggplot(l, ggplot2::aes_string("X1", "X2",
                                                      label = "Sample")) +
            ggplot2::geom_point(color = "grey")
        } else{
          g <- ggplot2::ggplot(l, ggplot2::aes_string("X1", "X2",
                                                      label = "Sample",
                                                      color = "expression")) +
            ggplot2::scale_colour_gradient(limits = c(min(l$expression),
                                                      max(l$expression)),
                                           low = "grey", high = "blue") +
            ggplot2::geom_point()
        }
        g <- g + ggplot2::labs(color = "Expression")
      }
      g <- g +
        ggplot2::ggtitle(paste(geneName, " - ", percent, "%", " cells",
                               sep = "")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
    if (!is.null(shape)){
      g <- g + ggplot2::aes_string(shape = "shape") +
        ggplot2::labs(shape = shape)
    }
    if (i == 1) {
      plist <- list(g)
    } else{
      plist <- cbind(plist, list(g))
    }
  }
  return(grid::grid.draw(gridExtra::arrangeGrob(
    grobs = plist, ncol = ceiling(sqrt(length(gene))))))
}
