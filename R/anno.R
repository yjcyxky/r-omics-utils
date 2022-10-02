#' Annotate Gene ID
#'
#' Annotate ensembl id
#'
#' @param geneids Ensembl ids
#' @param species Only support "Homo.sapiens" or "Mus.musculus"
#' @return A data frame
#' @export
annot <- function(geneids, species = "Homo.sapiens") {
  if (species == "Mus.musculus") {
    return(AnnotationDbi::select(Mus.musculus, keys = geneids, columns = c("ENSEMBL", "SYMBOL", "ENTREZID"), keytype = "ENSEMBL"))
  } else {
    return(AnnotationDbi::select(Homo.sapiens, keys = geneids, columns = c("ENSEMBL", "SYMBOL", "ENTREZID"), keytype = "ENSEMBL"))
  }
}

#' Annotate Gene ID
#'
#' Annotate ensembl id
#'
#' @param geneids Ensembl ids
#' @param species On support "Homo.sapiens" or "Mus.musculus"
#' @return A data frame
#' @export
annot_align <- function(geneids, species = "Homo.sapiens") {
  genes <- annot(geneids, species)
  duplicated_genes <- genes$ENSEMBL[duplicated(genes$ENSEMBL)]
  uniq_genes <- genes[!duplicated(genes$ENSEMBL), ]
  for (id in duplicated_genes) {
    uniq_genes$ENTREZID[which(uniq_genes$ENSEMBL == id)] <- paste0(genes[genes$ENSEMBL == id, ]$ENTREZID, collapse = ";")
    uniq_genes$SYMBOL[which(uniq_genes$ENSEMBL == id)] <- paste0(genes[genes$ENSEMBL == id, ]$SYMBOL, collapse = ";")
  }
  return(uniq_genes)
}

# library(Mus.musculus)
# fpkm <- read.csv("../../data/rawdata/pm25/mouse_pm2_5_89samples_fpkm.csv")
# geneids <- fpkm$gene_id
# genes <- annot_align(geneids)
