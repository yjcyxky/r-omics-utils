#' Annotate Gene ID
#'
#' Annotate ensembl id
#'
#' @param geneids Ensembl ids
#' @param species Only support "Homo.sapiens" or "Mus.musculus"
#' @param keytype Only support "ENSEMBL", "ENTREZID" or "SYMBOL"
#' @return A data frame
#' @export
annot <- function(geneids, species = "Homo.sapiens", keytype = "ENSEMBL") {
  if (species == "Mus.musculus") {
    return(AnnotationDbi::select(Mus.musculus, keys = geneids, columns = c("ENSEMBL", "SYMBOL", "ENTREZID"), keytype = keytype))
  } else {
    return(AnnotationDbi::select(Homo.sapiens, keys = geneids, columns = c("ENSEMBL", "SYMBOL", "ENTREZID"), keytype = keytype))
  }
}

#' Annotate Gene ID
#'
#' Annotate ensembl id
#'
#' @param geneids Ensembl ids
#' @param species On support "Homo.sapiens" or "Mus.musculus"
#' @param keytype Only support "ENSEMBL", "ENTREZID" or "SYMBOL"
#' @return A data frame
#' @export
annot_align <- function(geneids, species = "Homo.sapiens", keytype = "ENSEMBL") {
  genes <- annot(geneids, species)

  if (!(keytype %in% c("ENSEMBL", "ENTREZID", "SYMBOL"))) {
    stop("Keytype is invalid.")
  } else {
    duplicated_genes <- genes[,keytype][duplicated(genes[,keytype])]
    uniq_genes <- genes[!duplicated(genes[,keytype]), ]
    for (id in duplicated_genes) {
      uniq_genes$ENTREZID[which(uniq_genes[,keytype] == id)] <- paste0(genes[genes[,keytype] == id, ]$ENTREZID, collapse = ";")
      uniq_genes$SYMBOL[which(uniq_genes[,keytype] == id)] <- paste0(genes[genes[,keytype] == id, ]$SYMBOL, collapse = ";")
      uniq_genes$ENSEMBL[which(uniq_genes[,keytype] == id)] <- paste0(genes[genes[,keytype] == id, ]$ENSEMBL, collapse = ";")
    }
    return(uniq_genes)
  }
}

# library(Mus.musculus)
# fpkm <- read.csv("../../data/rawdata/pm25/mouse_pm2_5_89samples_fpkm.csv")
# geneids <- fpkm$gene_id
# genes <- annot_align(geneids)
