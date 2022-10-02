#' Make a pathway table
#'
#' @param species.KEGG three-letter KEGG species identifier, only support mmu and hsa.
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom purrr map
#' @importFrom limma getGeneKEGGLinks
#' @importFrom limma getKEGGPathwayNames
#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @export
make.pathwaytable <- function(species.KEGG="mmu") {
  tab <- getGeneKEGGLinks(species.KEGG=species.KEGG, convert=TRUE)
  if (species.KEGG == "mmu") {
    tab$Symbol <- mapIds(org.Mm.eg.db, keys=tab$GeneID,
                         column="SYMBOL", keytype="ENTREZID")
    tab$ENSEMBL_ID <- mapIds(org.Mm.eg.db, keys=tab$GeneID,
                             column="ENSEMBL", keytype="ENTREZID")
  } else if (species.KEGG == "hsa") {
    tab$Symbol <- mapIds(org.Hs.eg.db, keys=tab$GeneID,
                         column="SYMBOL", keytype="ENTREZID")
    tab$ENSEMBL_ID <- mapIds(org.Hs.eg.db, keys=tab$GeneID,
                             column="ENSEMBL", keytype="ENTREZID")    
  } else {
    stop("Don't know the species.")
  }
  
  pathway_names <- getKEGGPathwayNames(species=species.KEGG)
  ordered_pathway_names <- unlist(map(tab$PathwayID, function(id) return(pathway_names[which(pathway_names$PathwayID == id),]$Description)))
  tab$PathwayName <- ordered_pathway_names
  return(tab)
}

#' Read pathway file which contains two columns pathway_id and pathway_name.
#'
#' @param file file path
#' @export
read.pathways <- function(file) {
  pathways <- read.csv(file, header=TRUE, sep=",", quote="\"")
  if (!(length(pathways$pathway_name) > 0 && length(pathways$pathway_id) > 0)) {
    stop("Wrong format, need to contain pathway_id and pathway_name columns.")
  }
  return(pathways)
}

#' Get the genes in a pathway.
#'
#' @param pathways a data frame which contains two columns pathway_id and pathway_name
#' @importFrom stringr str_starts
#' @importFrom dplyr bind_rows
#' @export
match.genelists <- function(pathways) {
  if (all(str_starts(pathways$pathway_id, "hsa"))) {
    species.KEGG <- "hsa"
  } else if (all(str_starts(pathways$pathway_id, "mmu"))) {
    species.KEGG <- "mmu"
  } else {
    stop("Not valid pathway id or all pathway ids are not from the same species.")
  }

  pathway_table <- make.pathwaytable(species.KEGG=species.KEGG)
  return(bind_rows(map(pathways$pathway_id, function(id) return(pathway_table[which(pathway_table$PathwayID == paste("path:", id, sep="")),]))))
}