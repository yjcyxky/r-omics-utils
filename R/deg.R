#' This function is for identify DEG from a given matrix
#'
#' @param expr_mat should a matrix in stead of a data frame, with column for samples and rows for signature(gene/protein/metabolics...)
#'                 For transcriptome, a matrix after log transformation is recommended.
#' @param group    should be a factor whose length is identical to the number of the columns in expr_mat,
#'                 describing the group information of each column in expr_mat
#' @param min.count numeric. Minimum count required for at least some samples.
#' @importFrom edgeR DGEList
#' @importFrom edgeR filterByExpr
#' @importFrom edgeR calcNormFactors
#' @importFrom limma voom
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @importFrom data.table as.data.table
#' @export
DEGanalysis <- function(exprMat, group, min.count = 1) {
  dge <- DGEList(counts = exprMat)
  design <- model.matrix(~group)

  keep <- filterByExpr(dge, design = design, min.count = min.count)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)

  v <- voom(dge, design, plot = F)
  fit <- lmFit(v, design)

  fit <- eBayes(fit)
  result <- topTable(fit, coef = ncol(design), sort.by = "logFC", number = Inf)
  result$gene <- rownames(result)
  return(as.data.table(result))
}

# Testing
# d <- read.csv("~/Downloads/pm25/mouse_pm2_5_12samples_quartet_count.csv")
# rownames(d) <- d$gene_id
# d <- d[,2:dim(d)[2]]
# names <- colnames(d)[2:dim(d)[2]]
# group <- rep(c("A", "B"), c(6, 6))


#' Make groups based on metadata table
#'
#' @param metadata a data frame contains several columns for metadata.
#' @param group.column a field name.
#' @return a factor
#' @export
#' @examples
#' # Make groups
#' m <- data.frame(
#'   group = rep(c("FA", "PM"), 6),
#'   month = rep(5, 12), gender = rep("male", 12),
#'   organ = rep(c("gut", "lng", "lvr"), 4)
#' )
#' make_groups(m, group_column = "group", levels = c("PM", "FA"))
make_groups <- function(metadata, group_column = "group", levels = c("PM", "FA")) {
  cols <- colnames(metadata)
  if (group_column %in% cols && all(levels %in% cols)) {
    groups <- metadata[group_column]
    return(factor(groups, levels = levels))
  } else {
    stop(paste0("Not such group.column(", group_column, ") or group values(", levels, ")."))
  }
}

#' Compute logFC value.
#'
#' @param group1 a vector which contains several values (default is log transformed).
#' @param group2 a vector which contains several values (default is log transformed).
#' @param is_log boolean, default is TRUE. if your data is not log transformed, you need to set is_log to FALSE.
#' @return a logFC value
#' @export
logFC <- function(group1, group2, is_log = TRUE) {
  if (!is_log) {
    group1 <- log2(group1)
    group2 <- log2(group2)
  }
  mean(group1) - mean(group2)
  # log2(mean(group1)/mean(group2))
}

#' Performs one and two sample t-tests on vectors of data.
#'
#' @param gene a gene symbol/id.
#' @param x a (non-empty) numeric vector of data values..
#' @param y an optional (non-empty) numeric vector of data values.
#' @return a data frame which contains gene, P.Value, statistic, and logFC columns.
#' @export
ttest <- function(gene, x, y, paired = FALSE) {
  tryCatch(
    expr = {
      result <- t.test(x, y, paired = paired)
      result <- data.frame(gene = gene, P.Value = result$p.value, statistic = result$statistic, logFC = logFC(x, y))
      rownames(result) <- gene
      return(result)
    },
    error = function(e) {
      # message(e)
      result <- data.frame(gene = gene, P.Value = NA, statistic = NA, logFC = NA)
      rownames(result) <- gene
      return(result)
    }
  )
}

#' Performs one and two sample wilcox tests on vectors of data.
#'
#' @param gene a gene symbol/id.
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @return a data frame which contains gene, P.Value, statistic, and logFC columns.
#' @export
wilcox_test <- function(gene, x, y) {
  tryCatch(
    expr = {
      result <- wilcox.test(x, y)
      result <- data.frame(gene = gene, P.Value = result$p.value, statistic = result$statistic, logFC = logFC(x, y))
      rownames(result) <- gene
      return(result)
    },
    error = function(e) {
      # message(e)
      result <- data.frame(gene = gene, P.Value = NA, statistic = NA, logFC = NA)
      rownames(result) <- gene
      return(result)
    }
  )
}

#' Remove low expression genes.
#'
#' @param expr_mat a matrix of data values.
#' @param ratio the ratio which the gene expression values equal 0 is greater than the inputed value, the related genes will be removed.
#' @return a matrix which don't contain low expression genes.
#' @export
remove_genes <- function(expr_mat, ratio = 0.6) {
  saved_genes <- apply(expr_mat, 1, function(gene) {
    sum(gene == 0) / length(gene) < ratio
  })
  expr_mat[saved_genes, ]
}

#' Compute the DEGs
#'
#' @param expr_mat a data frame which columns are genes, rows are samples. the values are counts, tpm or fpkm.
#' @param groups a factor which contains group information.
#' @param pCutoff a cutoff for p.value.
#' @param log2FCCutoff a cutoff for log2FC.
#' @param method one of the t.test, wilcox.test, limma.
#' @param min.count numeric. Minimum count required for at least some samples (when counts).
#' @importFrom purrr map_df
#' @return a data frame which contains gene, P.Value, statistic, logFC and AveExpr columns.
#' @export
deg_stat <- function(expr_mat, groups, pCutoff = 0.05, log2FCCutoff = 1,
                     method = "t.test", min.count = 0) {
  samples <- colnames(expr_mat)
  if (!is.factor(groups) || (length(samples) != length(groups))) {
    stop("groups must be a factor and the length need to be same with colnames(expr_mat).")
  }

  print(paste("Run", method, "..."))
  if (method == "t.test" || method == "wilcox.test") {
    t_exp <- t(expr_mat)
    genes <- colnames(t_exp)
    test_group <- which(groups == levels(groups)[1])
    control_group <- which(groups == levels(groups)[2])
    d <- map_df(genes, function(gene) {
      v <- t_exp[, gene]
      x <- v[test_group]
      y <- v[control_group]
      if (method == "t.test") {
        ttest(gene, x, y)
      } else if (method == "wilcox.test") {
        wilcox_test(gene, x, y)
      }
    })
  } else {
    d <- DEGanalysis(expr_mat, groups, min.count = min.count)
    d <- data.frame(gene = d$gene, P.Value = d$P.Value, statistic = d$t, logFC = d$logFC, AveExpr = d$AveExpr)
  }

  results <- d
  results$direction <- "No"
  results$direction[which(d$P.Value <= pCutoff & d$logFC >= log2FCCutoff)] <- "Up"
  results$direction[which(d$P.Value <= pCutoff & d$logFC <= -log2FCCutoff)] <- "Down"

  return(results)
}

#' Given a set of p-values, returns p-values adjusted using one of several methods.
#'
#' @param pvalues numeric vector of p-values (possibly with NAs). Any other R object is coerced by as.numeric.
#' @param method correction method, a character string. Can be BH, BY, bonferroni, holmmel, hochberg, or none.
#' @return numeric vector of corrected p-values.
#' @export
adjust_p <- function(pvalues, method = "bonferroni") {
  return(p.adjust(pvalues, method = method))
}

#' Guess database which the gene ids come from.
#'
#' @param genes string vector of gene id.
#' @return database name, ENSEMBL, ENTREZID or SYMBOL.
#' @export
which_database <- function(genes) {
  matched <- grep("^EN.*", genes)
  count <- length(genes)
  if (length(matched) != 0) {
    if (length(matched) == count) {
      return("ENSEMBL")
    } else {
      stop("Found mixed ids, you need to use the same id database.")
    }
  }

  matched <- grep("^[0-9]+", genes)
  if (length(matched) != 0) {
    if (length(matched) == count) {
      return("ENTREZID")
    } else {
      stop("Found mixed ids, you need to use the same id database.")
    }
  }

  return("SYMBOL")
}

#' Annotate gene ids.
#'
#' @param genes string vector of gene id.
#' @param species Homo.sapiens or Mus.musculus, default is Mus.musculus.
#' @import Mus.musculus
#' @import Homo.sapiens
#' @import AnnotationDbi
#' @return a data frame which contains ensembl_id, entrez_id and gene_symbol columns.
#' @export
anno_genes <- function(genes, species = "Mus.musculus") {
  database <- which_database(genes)
  results <- data.frame()
  annotated_genes <- annot_align(genes, species, database)
  results$ensembl_id <- annotated_genes$ENSEMBL
  results$entrez_id <- annotated_genes$ENTREZID
  results$gene_symbol <- annotated_genes$SYMBOL
  results
}

#' Enrich pathways.
#'
#' @param entrez_ids string vector of gene entrez id.
#' @param species Homo.sapiens or Mus.musculus, default is Mus.musculus.
#' @param universe same with enrichGO
#' @param showCategory same with dotplot
#' @param pCutoff adjusted pvalue cutoff on enrichment tests to report
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param qvalueCutoff qvalue cutoff on enrichment tests to report as significant. Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues and iii) qvalueCutoff on qvalues to be reported.
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param filename_prefix a prefix of output filename.
#' @param output_folder which directory do you want to save to.
#' @import DOSE
#' @import topGO
#' @import clusterProfiler
#' @import pathview
#' @export
enricher <- function(entrez_ids, species = "Mus.musculus", universe = NULL,
                     showCategory = 20, pCutoff = 0.05, minGSSize = 10, qvalueCutoff = 0.2,
                     pAdjustMethod = "none", filename_prefix = "enricher", output_folder = ".") {
  organism <- "mmu"
  OrgDb <- "org.Mm.eg.db"
  if (species == "Mus.musculus") {
    organism <- "mmu"
    OrgDb <- "org.Mm.eg.db"
  } else {
    organism <- "hsa"
    OrgDb <- "org.Hs.eg.db"
  }

  kk <- enrichKEGG(gene = entrez_ids, organism = organism, universe = universe, pvalueCutoff = pCutoff, minGSSize = minGSSize, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff)
  if (dim(summary(kk))[1] != 0) {
    result <- as.data.frame(kk)
    write.csv(result, file.path(output_folder, paste0(filename_prefix, "_kegg.csv")))
    p <- dotplot(kk, showCategory = showCategory, font.size = 10, color = "pvalue", title = "Enrichment KEGG Dot")
    save2pdf(p, paste0(filename_prefix, "_kegg.pdf"), output_folder, width = 8, height = 8)
    save2pdf(p, paste0(filename_prefix, "_kegg.png"), output_folder, width = 8, height = 8)
  }

  goBP <- enrichGO(gene = entrez_ids, OrgDb = OrgDb, ont = "BP", pvalueCutoff = pCutoff, minGSSize = minGSSize, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff)
  if (dim(summary(goBP))[1] != 0) {
    result <- as.data.frame(goBP)
    write.csv(result, file.path(output_folder, paste0(filename_prefix, "_goBP.csv")))
    p <- dotplot(goBP, showCategory = showCategory, font.size = 10, color = "pvalue", title = "Enrichment GO BP Dot")
    save2pdf(p, paste0(filename_prefix, "_goBP.pdf"), output_folder, width = 8, height = 8)
    save2pdf(p, paste0(filename_prefix, "_goBP.png"), output_folder, width = 8, height = 8)
  }

  goMF <- enrichGO(gene = entrez_ids, OrgDb = OrgDb, ont = "MF", pvalueCutoff = pCutoff, minGSSize = minGSSize, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff)
  if (dim(summary(goMF))[1] != 0) {
    result <- as.data.frame(goMF)
    write.csv(result, file.path(output_folder, paste0(filename_prefix, "_goMF.csv")))
    p <- dotplot(goMF, showCategory = showCategory, font.size = 10, color = "pvalue", title = "Enrichment GO MF Dot")
    save2pdf(p, paste0(filename_prefix, "_goMF.pdf"), output_folder, width = 8, height = 8)
    save2pdf(p, paste0(filename_prefix, "_goMF.png"), output_folder, width = 8, height = 8)
  }

  goCC <- enrichGO(gene = entrez_ids, OrgDb = OrgDb, ont = "CC", pvalueCutoff = pCutoff, minGSSize = minGSSize, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff)
  if (dim(summary(goCC))[1] != 0) {
    result <- as.data.frame(goCC)
    write.csv(result, file.path(output_folder, paste0(filename_prefix, "_goCC.csv")))
    p <- dotplot(goCC, showCategory = showCategory, font.size = 10, color = "pvalue", title = "Enrichment GO CC Dot")
    save2pdf(p, paste0(filename_prefix, "_goCC.pdf"), output_folder, width = 8, height = 8)
    save2pdf(p, paste0(filename_prefix, "_goCC.png"), output_folder, width = 8, height = 8)
  }
}
