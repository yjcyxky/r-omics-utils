run_limma <- function(data_file, metadata, dest_dir = ".",
                      remove_low_expr_genes_fn = NULL,
                      group_column = "group", levels = c("PM", "FA"),
                      ...) {
  method <- "limma"
  method_abbr <- "counts_limma"
  suffix <- "_counts_limma"

  if (!file.exists(data_file)) {
    stop(paste0()("No such file: ", data_file))
  } else {
    raw_expr_mat <- read.csv(data_file, row.names = 1)
  }

  if (class(metadata)[1] != "Metadata") {
    stop("metadata is not a valid Metadata class.")
  }

  subsets <- unique(metadata$subsets)
  if (is.null(subsets)) {
    stop("No subsets in metadata, you need to run select_samples firstly.")
  }

  for (s in subsets) {
    sample_ids <- get_sample_ids(metadata, subset_value = s)
    expr_mat <- subset(raw_expr_mat, select = sample_ids)
    
    if (!is.null(remove_low_expr_genes_fn)) {
      print(paste("Remove low expression genes with function: "))
      print(remove_low_expr_genes_fn)
      expr_mat <- remove_low_expr_genes_fn(expr_mat)
    }
    
    rootdir <- file.path(dest_dir, s, method_abbr)
    folder <- file.path(rootdir, "degs")
    dir.create(file.path(rootdir, "degs"), recursive = TRUE)
    dir.create(file.path(rootdir, "barplot"), recursive = TRUE)
    dir.create(file.path(rootdir, "enrichment"), recursive = TRUE)
    dir.create(file.path(rootdir, "heatmap"), recursive = TRUE)
    dir.create(file.path(rootdir, "volcano"), recursive = TRUE)

    groups <- make_groups(metadata, s, group_column = group_column, levels = levels)
    r <- analyze_degs(expr_mat,
      groups = groups, method = method, suffix = suffix,
      folder = folder, ...
    )

    write.csv(r$all, file.path(folder, paste0(s, suffix, "_all", ".csv")), row.names = FALSE, quote = FALSE)
    write.csv(r$degs, file.path(folder, paste0(s, suffix, ".csv")), row.names = FALSE, quote = FALSE)

    # Barplot
    input_folder <- file.path(rootdir, "degs")
    output_folder <- file.path(rootdir, "barplot")
    deg <- read.csv(file.path(input_folder, paste0(s, suffix, ".csv")), row.names = 1)
    p <- draw_barplot(deg, s)
    save2pdf(p, paste0(s, suffix, "_barplot.pdf"), output_folder, width = 8, height = 8)

    # Heatmap
    input_folder <- file.path(rootdir, "degs")
    output_folder <- file.path(rootdir, "heatmap")
    deg <- read.csv(file.path(input_folder, paste0(s, suffix, ".csv")), row.names = 1)
    genes <- rownames(deg)
    p <- draw_heatmap(expr_mat, c(), genes)
    save2pdf(p, paste0(s, suffix, "_heatmap.pdf"), output_folder, width = 8, height = 8)

    # Volcano
    input_folder <- file.path(rootdir, "degs")
    output_folder <- file.path(rootdir, "volcano")
    deg <- read.csv(file.path(input_folder, paste0(s, suffix, "_all.csv")), row.names = 1)
    p <- draw_volcano(deg, paste0("Volcano of ", s))
    save2pdf(p, paste0(s, suffix, "_volcano.pdf"), output_folder, width = 8, height = 8)

    deg <- read.csv(file.path(input_folder, paste0(s, suffix, "_all.csv")), row.names = 1)
    p <- draw_volcano_ggplot(deg, paste0("Volcano of ", s))
    save2pdf(p, paste0(s, suffix, "_volcano_ggplot.pdf"), output_folder, width = 8, height = 8)

    # Enrichment
    input_folder <- file.path(rootdir, "degs")
    output_folder <- file.path(rootdir, "enrichment")
    deg <- read.csv(file.path(input_folder, paste0(s, suffix, ".csv")), row.names = 1)
    validGenes <- deg$entrez_id[which(!is.na(deg$entrez_id))]
    validGenes <- unlist(lapply(as.character(validGenes), function(id) {
      unlist(strsplit(id, split = ";"))
    }))
    enricher(validGenes, output_folder = output_folder, filename_prefix = paste(s, suffix, sep = ""))
  }
}


run_ttest <- function(data_file, metadata, dest_dir = ".",
                      remove_low_expr_genes_fn = NULL, data_type = "fpkm",
                      plus_num = 0.01, enable_log2 = TRUE,
                      group_column = "group", levels = c("PM", "FA"),
                      ...) {
  method <- "t.test"
  method_abbr <- paste(data_type, "_ttest", sep = "")
  suffix <- paste0(c("_", data_type, "_ttest"), collapse = "")

  if (!file.exists(data_file)) {
    stop(paste0()("No such file: ", data_file))
  } else {
    raw_expr_mat <- read.csv(data_file, row.names = 1)
  }

  if (class(metadata)[1] != "Metadata") {
    stop("metadata is not a valid Metadata class.")
  }

  subsets <- unique(metadata$subsets)
  if (is.null(subsets)) {
    stop("No subsets in metadata, you need to run select_samples firstly.")
  }

  for (s in subsets) {
    sample_ids <- get_sample_ids(metadata, subset_value = s)
    expr_mat <- subset(raw_expr_mat, select = sample_ids)
    
    if (!is.null(remove_low_expr_genes_fn)) {
      print(paste("Remove low expression genes with function: "))
      print(remove_low_expr_genes_fn)
      expr_mat <- remove_low_expr_genes_fn(expr_mat)
    }
    
    rootdir <- file.path(dest_dir, s, method_abbr)
    folder <- file.path(rootdir, "degs")
    dir.create(file.path(rootdir, "degs"), recursive = TRUE)
    dir.create(file.path(rootdir, "barplot"), recursive = TRUE)
    dir.create(file.path(rootdir, "enrichment"), recursive = TRUE)
    dir.create(file.path(rootdir, "heatmap"), recursive = TRUE)
    dir.create(file.path(rootdir, "volcano"), recursive = TRUE)

    expr_mat <- expr_mat + plus_num

    if (enable_log2) {
      expr_mat <- log2(expr_mat)
    }

    groups <- make_groups(metadata, s, group_column = group_column, levels = levels)
    r <- analyze_degs(expr_mat,
      groups = groups, method = method, suffix = suffix,
      folder = folder, ...
    )

    write.csv(r$all, file.path(folder, paste0(s, suffix, "_all", ".csv")), row.names = FALSE, quote = FALSE)
    write.csv(r$degs, file.path(folder, paste0(s, suffix, ".csv")), row.names = FALSE, quote = FALSE)

    # Barplot
    input_folder <- file.path(rootdir, "degs")
    output_folder <- file.path(rootdir, "barplot")
    deg <- read.csv(file.path(input_folder, paste0(s, suffix, ".csv")), row.names = 1)
    p <- draw_barplot(deg, s)
    save2pdf(p, paste0(s, suffix, "_barplot.pdf"), output_folder, width = 8, height = 8)

    # Heatmap
    input_folder <- file.path(rootdir, "degs")
    output_folder <- file.path(rootdir, "heatmap")
    deg <- read.csv(file.path(input_folder, paste0(s, suffix, ".csv")), row.names = 1)
    genes <- rownames(deg)
    if (enable_log2) {
      p <- draw_heatmap(2 ^ expr_mat, c(), genes, scale = "column")
    } else {
      p <- draw_heatmap(expr_mat, c(), genes, scale = "column")
    }
    save2pdf(p, paste0(s, suffix, "_heatmap.pdf"), output_folder, width = 8, height = 8)

    # Volcano
    input_folder <- file.path(rootdir, "degs")
    output_folder <- file.path(rootdir, "volcano")
    deg <- read.csv(file.path(input_folder, paste0(s, suffix, "_all.csv")), row.names = 1)
    p <- draw_volcano(deg, paste0("Volcano of ", s))
    save2pdf(p, paste0(s, suffix, "_volcano.pdf"), output_folder, width = 8, height = 8)

    deg <- read.csv(file.path(input_folder, paste0(s, suffix, "_all.csv")), row.names = 1)
    p <- draw_volcano_ggplot(deg, paste0("Volcano of ", s))
    save2pdf(p, paste0(s, suffix, "_volcano_ggplot.pdf"), output_folder, width = 8, height = 8)

    # Enrichment
    input_folder <- file.path(rootdir, "degs")
    output_folder <- file.path(rootdir, "enrichment")
    deg <- read.csv(file.path(input_folder, paste0(s, suffix, ".csv")), row.names = 1)
    validGenes <- deg$entrez_id[which(!is.na(deg$entrez_id))]
    validGenes <- unlist(lapply(as.character(validGenes), function(id) {
      unlist(strsplit(id, split = ";"))
    }))
    enricher(validGenes, output_folder = output_folder, filename_prefix = paste(s, suffix, sep = ""))
  }
}

run_wilcox <- function(data_file, metadata, dest_dir = ".",
                       remove_low_expr_genes_fn = NULL, data_type = "fpkm",
                       plus_num = 0.01, enable_log2 = TRUE,
                       group_column = "group", levels = c("PM", "FA"),
                       ...) {
  method <- "wilcox.test"
  method_abbr <- paste(data_type, "_wilcox", sep = "")
  suffix <- paste0(c("_", data_type, "_wilcox"), collapse = "")

  if (!file.exists(data_file)) {
    stop(paste0()("No such file: ", data_file))
  } else {
    raw_expr_mat <- read.csv(data_file, row.names = 1)
  }

  if (class(metadata)[1] != "Metadata") {
    stop("metadata is not a valid Metadata class.")
  }

  subsets <- unique(metadata$subsets)
  if (is.null(subsets)) {
    stop("No subsets in metadata, you need to run select_samples firstly.")
  }

  for (s in subsets) {
    sample_ids <- get_sample_ids(metadata, subset_value = s)
    expr_mat <- subset(raw_expr_mat, select = sample_ids)
    
    if (!is.null(remove_low_expr_genes_fn)) {
      print(paste("Remove low expression genes with function: "))
      print(remove_low_expr_genes_fn)
      expr_mat <- remove_low_expr_genes_fn(expr_mat)
    }
    
    rootdir <- file.path(dest_dir, s, method_abbr)
    folder <- file.path(rootdir, "degs")
    dir.create(file.path(rootdir, "degs"), recursive = TRUE)
    dir.create(file.path(rootdir, "barplot"), recursive = TRUE)
    dir.create(file.path(rootdir, "enrichment"), recursive = TRUE)
    dir.create(file.path(rootdir, "heatmap"), recursive = TRUE)
    dir.create(file.path(rootdir, "volcano"), recursive = TRUE)

    expr_mat <- expr_mat + plus_num

    if (enable_log2) {
      expr_mat <- log2(expr_mat)
    }

    groups <- make_groups(metadata, s, group_column = group_column, levels = levels)
    r <- analyze_degs(expr_mat,
      groups = groups, method = method, suffix = suffix,
      folder = folder, ...
    )

    write.csv(r$all, file.path(folder, paste0(s, suffix, "_all", ".csv")), row.names = FALSE, quote = FALSE)
    write.csv(r$degs, file.path(folder, paste0(s, suffix, ".csv")), row.names = FALSE, quote = FALSE)

    # Barplot
    input_folder <- file.path(rootdir, "degs")
    output_folder <- file.path(rootdir, "barplot")
    deg <- read.csv(file.path(input_folder, paste0(s, suffix, ".csv")), row.names = 1)
    p <- draw_barplot(deg, s)
    save2pdf(p, paste0(s, suffix, "_barplot.pdf"), output_folder, width = 8, height = 8)

    # Heatmap
    input_folder <- file.path(rootdir, "degs")
    output_folder <- file.path(rootdir, "heatmap")
    deg <- read.csv(file.path(input_folder, paste0(s, suffix, ".csv")), row.names = 1)
    genes <- rownames(deg)
    
    if (enable_log2) {
      p <- draw_heatmap(2 ^ expr_mat, c(), genes, scale = "column")
    } else {
      p <- draw_heatmap(expr_mat, c(), genes, scale = "column")
    }
    
    save2pdf(p, paste0(s, suffix, "_heatmap.pdf"), output_folder, width = 8, height = 8)

    # Volcano
    input_folder <- file.path(rootdir, "degs")
    output_folder <- file.path(rootdir, "volcano")
    deg <- read.csv(file.path(input_folder, paste0(s, suffix, "_all.csv")), row.names = 1)
    p <- draw_volcano(deg, paste0("Volcano of ", s))
    save2pdf(p, paste0(s, suffix, "_volcano.pdf"), output_folder, width = 8, height = 8)

    deg <- read.csv(file.path(input_folder, paste0(s, suffix, "_all.csv")), row.names = 1)
    p <- draw_volcano_ggplot(deg, paste0("Volcano of ", s))
    save2pdf(p, paste0(s, suffix, "_volcano_ggplot.pdf"), output_folder, width = 8, height = 8)

    # Enrichment
    input_folder <- file.path(rootdir, "degs")
    output_folder <- file.path(rootdir, "enrichment")
    deg <- read.csv(file.path(input_folder, paste0(s, suffix, ".csv")), row.names = 1)
    validGenes <- deg$entrez_id[which(!is.na(deg$entrez_id))]
    validGenes <- unlist(lapply(as.character(validGenes), function(id) {
      unlist(strsplit(id, split = ";"))
    }))
    enricher(validGenes, output_folder = output_folder, filename_prefix = paste(s, suffix, sep = ""))
  }
}
