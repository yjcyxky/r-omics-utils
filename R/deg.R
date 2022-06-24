#' This function is for identify DEG from a given matrix
#'
#' @param expr_mat should a matrix in stead of a data frame, with column for samples and rows for signature(gene/protein/metabolics...)
#'                 For transcriptome, a matrix after log transformation is recommended.
#' @param group    should be a factor whose length is identical to the number of the columns in expr_mat,
#'                 describing the group information of each column in expr_mat
#' @importFrom edgeR DGEList
#' @importFrom edgeR filterByExpr
#' @importFrom edgeR calcNormFactors
#' @importFrom limma voom
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @importFrom data.table as.data.table
#' @export
DEGanalysis <- function(exprMat, group){
  dge <- DGEList(counts = exprMat)
  design <- model.matrix(~group)

  keep <- filterByExpr(dge, design, design,min.count = 0)
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)

  v <- voom(dge, design, plot=F)
  fit <- lmFit(v, design)

  fit <- eBayes(fit)
  result <- topTable(fit, coef=ncol(design), sort.by = 'logFC', number = Inf)
  result$gene = rownames(result)
  result$groupA =  levels(group)[1]
  result$groupB =  levels(group)[2]
  return(as.data.table(result))
}

# Testing
# d <- read.csv("~/Downloads/pm25/mouse_pm2_5_12samples_quartet_count.csv")
# rownames(d) <- d$gene_id
# d <- d[,2:dim(d)[2]]
# names <- colnames(d)[2:dim(d)[2]]
# group <- rep(c("A", "B"), c(6, 6))