# library(ggplot2)
# library(ggpubr)

#' Draw boxplot
#'
#' This function need a dataframe which contains three columns (gene_symbol, group, organ, value).
#' The organ column is optional and the gene symbol column maybe have several genes. The value column
#' has mean values of the related gene's expression value.
#' More details on
#' 1. Points: https://r-graph-gallery.com/89-box-and-scatter-plot-with-ggplot2.html
#' 2. Grouped Boxplot: https://r-graph-gallery.com/265-grouped-boxplot-with-ggplot2.html
#' 3. P-Value: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
#'
#' @param d a dataframe which contains these columns (gene_symbol, group, organ, value).
#' @param method ("t.test", "wilcox.test", "anova", "kruskal.test")
#' @return ggplot object
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 geom_jitter
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 aes
#' @importFrom ggpubr stat_compare_means
#' @export
boxplot <- function(d, method = "t.test", jitter_size = 0.4, log_scale = FALSE) {
  if (log_scale) {
    d$value <- log2(d$value)
    ytitle <- "Expression [log2(TPM + 1)]"
  } else {
    ytitle <- "Expression [TPM]"
  }

  if (!is.null(d$organ)) {
    # One gene in multiple organs
    if (length(levels(factor(d$gene_symbol))) > 1) {
      stop("Don't support multiple genes when in organ mode.")
    } else {
      ytitle <- paste(d$gene_symbol, ytitle, sep = " ")
      ggplot(d, aes(x = organ, y = value, fill = group)) +
        geom_boxplot() +
        geom_jitter(color = "black", size = jitter_size, alpha = 0.9) +
        stat_compare_means(method = method) +
        labs(x = "Organ", y = ytitle, fill = "Group") + custom_theme()
    }
  } else {
    # Multiple genes in one organ
    ggplot(d, aes(x = gene_symbol, y = value, fill = group)) +
      geom_boxplot() +
      geom_jitter(color = "black", size = jitter_size, alpha = 0.9) +
      stat_compare_means(method = method) +
      labs(x = "Gene Symbol", y = ytitle, fill = "Group") + custom_theme()
  }
}

# Testing
# source("./R/theme.R")
# # Multiple genes in one organ
# d <- data.frame(gene_symbol=rep(c("TP53", "ERBB2"), c(12, 12)),
#                 group=rep(c("Control", "Test", "Control", "Test"), c(6, 6, 6, 6)),
#                 value=round(runif(24), 1))
# boxplot(d) + custom_theme

# # One gene in multiple organs
# d <- data.frame(gene_symbol=rep("TP53", 24),
#                 organ=rep(c("OV", "ON"), c(12, 12)),
#                 group=rep(c("Control", "Test", "Control", "Test"), c(6, 6, 6, 6)),
#                 value=round(runif(24), 1))
# boxplot(d) + custom_theme


#' Query DuckDB database
#'
#' Query parquet files with SQL query string.
#' 
#' @param db_path A file path for database file, such as my-db.duckdb
#' @param query_str A query string
#' @return A data frame
#' @importFrom DBI dbConnect
#' @importFrom DBI dbExecute
#' @importFrom DBI dbGetQuery
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom duckdb duckdb
#' @export
query_db4boxplot <- function(db_path, query_str) {
  # gene_symbol, group, value, organ
  con = dbConnect(duckdb(), dbdir=db_path, read_only=TRUE)
  res = dbGetQuery(con, query_str)
  return(res)
}