#' Remove zeros
#'
#' @param expMatrix matrix
#' @return a matrix
#' @export
remove_zero <- function(expMatrix) {
  return(expMatrix[which(apply(expMatrix, 1, var) != 0), ])
}

#' Draw pca plot
#'
#' @param expMatrix matrix
#' @param is_qc boolean
#' @param scale a logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place. The default is FALSE for consistency with S, but in general scaling is advisable. Alternatively, a vector of length equal the number of columns of x can be supplied.
#' @importFrom dplyr %>%
#' @import ggplot2
#' @return ggplot2 object
#' @export
draw_pca <- function(expMatrix, is_qc = FALSE, scale = FALSE) {
  if (scale) {
    expMatrix <- remove_zero(expMatrix)
    pca_prcomp <- prcomp(t(expMatrix), scale. = TRUE)
  } else {
    pca_prcomp <- prcomp(t(expMatrix))
  }
  pcs <- predict(pca_prcomp) %>% data.frame()

  fields <- colnames(expMatrix) %>%
    strsplit("_") %>%
    as.data.frame()

  if (is_qc) {
    pcs$sample <- as.vector(as.matrix(fields[1, ]))
    pcs$group <- as.vector(as.matrix(fields[2, ]))
    shape_level <- nlevels(pcs$sample)
    if (shape_level < 15) {
      shapes <- (0:shape_level) %% 15
    } else {
      shapes <- c(0:14, c((15:shape_level) %% 110 + 18))
    }
    point <- geom_point(aes(color = group, shape = sample), size = 2.5)
    shape <- scale_shape_manual(values = shapes)
  } else {
    pcs$organ <- as.vector(as.matrix(fields[1, ]))
    pcs$group <- as.vector(as.matrix(fields[2, ]))
    pcs$sample <- colnames(expMatrix)
    print(nlevels(as.factor(pcs$organ)))
    shape_level <- nlevels(as.factor(pcs$organ))
    if (shape_level < 7) {
      shape <- scale_shape()
    } else {
      shapes <- c(0:6, c((7:shape_level) %% 110 + 18))
      shape <- scale_shape_manual(values = shapes)
    }
    point <- geom_point(aes(color = group, fill = group, shape = organ), size = 2.5)
  }

  scale.axis.x <- c(min(pcs$PC1), max(pcs$PC1))
  scale.axis.y <- c(min(pcs$PC2), max(pcs$PC2))

  ggplot(pcs, aes(x = PC1, y = PC2)) +
    point +
    shape +
    theme_few() +
    theme(
      axis.text = element_text(face = "bold", color = "gray40"),
      axis.title = element_text(face = "bold"),
      legend.text = element_text(color = "gray40"),
      plot.subtitle = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_x_continuous(limits = c(1.1 * scale.axis.x[1], 1.1 * scale.axis.x[2])) +
    scale_y_continuous(limits = c(1.1 * scale.axis.y[1], 1.1 * scale.axis.y[2])) +
    labs(
      x = sprintf("PC1 (%.0f%%)", summary(pca_prcomp)$importance[2, 1] * 100),
      y = sprintf("PC2 (%.0f%%)", summary(pca_prcomp)$importance[2, 2] * 100),
      title = "log2Expr"
    ) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    guides(shape = guide_legend(override.aes = list(size = 3)))
}

#' Draw volcano plot
#'
#' @param deg a data frame
#' @param title plot title
#' @param pCutoff cut-off for statistical significance. A horizontal line will be drawn at -log10(pCutoff).
#' @param FCcutoff cut-off for absolute log2 fold-change. Vertical lines will be drawn at the negative and positive values of log2FCcutoff.
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @return ggplot object
#' @export
draw_volcano <- function(deg, title, pCutoff = 0.05, FCcutoff = 1) {
  EnhancedVolcano(deg,
    lab = rownames(deg),
    x = "logFC",
    y = "P.Value",
    selectLab = "",
    title = title,
    cutoffLineType = "twodash",
    cutoffLineWidth = 0.8,
    pCutoff = pCutoff, ## pvalue阈值
    FCcutoff = FCcutoff, ## FC cutoff
    pointSize = 1.5,
    labSize = 6.0,
    col = c("grey30", "forestgreen", "royalblue", "red2"),
    colAlpha = 1,
    legendPosition = "bottom",
    legendLabSize = 10,
    legendIconSize = 5.0
  )
}

#' Draw barplot
#'
#' @param deg a data frame
#' @param ptitle plot title
#' @import ggplot2
#' @return ggplot object
#' @export
draw_barplot <- function(deg, ptitle) {
  up <- dim(deg[which(deg$logFC > 0), ])[1]
  down <- dim(deg[which(deg$logFC < 0), ])[1]
  max_count <- dim(deg)[1] + 100
  data <- data.frame(count = c(up, down), name = factor(c("Up-regulation", "Down-regulation"), levels = c("Up-regulation", "Down-regulation")))
  ggplot(data, aes(x = name, y = count, fill = name)) +
    geom_bar(stat = "identity", alpha = .6, width = .4) +
    ylab("Number of DEGs") +
    xlab("") +
    ggtitle(ptitle) +
    scale_y_continuous(limits = c(0, max_count), breaks = seq(0, max_count, if (max_count > 1000) 200 else 100)) +
    custom_theme()
}

#' Draw heatmap
#'
#' @param expMatrix matrix
#' @param removed_samples Which samples to remove
#' @param genes selected genes
#' @importFrom pheatmap pheatmap
#' @import ggplot2
#' @return ggplot2 object
#' @export
draw_heatmap <- function(exp_mat, removed_samples = c(), genes = c(), scale = "row") {
  if (length(removed_samples) > 0) {
    exp_mat <- exp_mat[, -which(colnames(exp_mat) %in% removed_samples)]
  }

  if (length(genes) == 0) {
    deg_exp_mat <- exp_mat
  } else {
    deg_exp_mat <- exp_mat[genes, ]
  }

  deg_exp_mat <- log2(deg_exp_mat + 0.01)
  dat <- as.matrix(deg_exp_mat)

  color <- colorRampPalette(colors = c("red", "black", "green"))(500)
  dat <- dat[apply(dat, 1, function(x) sd(x) != 0), ]
  pheatmap(dat,
    color = color, angle_col = 45, fontsize_col = 6,
    scale = scale, cluster_rows = TRUE, show_rownames = FALSE, silent = TRUE
  )
}

#' Draw volcana with ggplot2
#'
#' @param deg degs
#' @param title plot title
#' @param pCutoff cut-off for statistical significance. A horizontal line will be drawn at -log10(pCutoff).
#' @param FCcutoff cut-off for absolute log2 fold-change. Vertical lines will be drawn at the negative and positive values of log2FCcutoff.
#' @import ggplot2
#' @return ggplot2 object
#' @export
draw_volcano_ggplot <- function(deg, title, pCutoff = 0.05, FCcutoff = 1) {
  deg <- deg[which(!is.na(deg$P.Value) & !is.na(deg$logFC)), ]
  deg$threshold <- as.factor(ifelse(deg$P.Value < pCutoff & abs(deg$logFC) >= FCcutoff, ifelse(deg$logFC > FCcutoff, "Up", "Down"), "NoSignificant"))
  ggplot(data = deg, aes(x = logFC, y = -log10(P.Value), colour = threshold)) +
    geom_point(size = 1.5) +
    scale_color_manual(values = c("blue", "grey", "red")) +
    xlim(-4, 4) +
    ylim(c(0, max(-log10(deg$P.Value)) + 1)) +
    geom_vline(xintercept = c(-FCcutoff, FCcutoff), lty = 4, col = "black", lwd = 0.8) +
    geom_hline(yintercept = -log10(pCutoff), lty = 4, col = "black", lwd = 0.8) +
    labs(x = "Log2(Fold Change)", y = "-Log10(P-value)", title = title) +
    custom_theme()
}
