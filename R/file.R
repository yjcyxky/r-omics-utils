#' Write ggplot object to a file
#'
#' @param plot ggplot object.
#' @param filename a file name.
#' @param filepath a directory.
#' @param width image width.
#' @param height image height.
#' @param units Plot size in units ("in", "cm", "mm", or "px"). If not supplied, uses the size of current graphics device.
#' @importFrom ggplot2 ggsave
#' @export
save2pdf <- function(plot, filename, filepath, width = 4, height = 4, units = "in") {
  ggsave(filename = filename, plot = plot, path = filepath, width = width, height = height, units = units)
}

#' Metadata S3 class
#'
#' @param metadata_file a csv/tsv file
#' @importFrom purrr partial
#' @return a s3 class
#' @export
#' @examples 
#' # List all example files
#' system.file("extdata", package = "rOmicsUtils") |> list.files()
#' 
#' metadata <- rOmicsUtils::Metadata(file.path(system.file("extdata", package = "rOmicsUtils")[1], "GSE146508.csv"))
#' select_samples(metadata, subset_fields = c("organ", "name"))
#' make_groups(metadata, group_column = "group_name", subset_value = "Bld_Mutlu")
#' get_subset(metadata, subset = "Bld_Mutlu")
#' 
#' select_samples(metadata, sql_query_str = "organ = 'Bld' AND name = 'Biswal'", subset_fields = c("organ", "name"))
Metadata <- function(metadata_file) {
  if (grep(".*.csv", metadata_file)) {
    metadata <- read.csv(metadata_file, header=TRUE)
  } else if (grep(".*.tsv", metadata_file)) {
    metadata <- read.table(metadata_file, sep="\t", header=TRUE)
  } else {
    stop("Only support csv or tsv file.")
  }

  structure(local({
    metadata = metadata
    filtered_metadata = NULL
    subsets = NULL
    subset_fields = NULL
    environment()
  }), class = c("Metadata", "environment"))
}

get_sample_ids <- function(object, ...) {
  UseMethod("get_sample_ids", object)
}

get_sample_ids.default <- function(object, ...) {
  print("The class of this object can not be found.")
}

#' Get sample ids
#'
#' @param metadata a metadata object.
#' @return a vector
#' @export
get_sample_ids.Metadata <- function(m, subset_value=NULL) {
  if (is.null(m$filtered_metadata)) {
    warning("Use the full metadata table. if you want to use subset, you need to run select_samples firstly.")
    metadata <- m$metadata
  } else {
    metadata <- m$filtered_metadata
  }

  if (!is.null(subset_value)) {
    # Select a subset
    metadata[which(metadata$subset == subset_value),]$sample_id  
  } else {
    metadata$sample_id
  }
}

get_subset <- function(object, ...) {
  UseMethod("get_subsets", object)
}

get_subset.default <- function(object, ...) {
  print("The class of this object can not be found.")
}

#' Get a subset of metadata
#'
#' @param metadata a metadata object.
#' @param subset which subset value in subset column.
#' @return a data frame
#' @export
get_subset.Metadata <- function(m, subset_value=NULL) {
  if (is.null(m$filtered_metadata)) {
    warning("Use the full metadata table. if you want to use subset, you need to run select_samples firstly.")
    metadata <- m$metadata
  } else {
    metadata <- m$filtered_metadata
  }
  
  if (!is.null(subset_value)) {
    # Select a subset
    metadata[which(metadata$subset == subset_value),]
  } else {
    metadata
  }
}

make_groups <- function(object, ...) {
  UseMethod("make_groups", object)
}

make_groups.default <- function(object, ...) {
  print("The class of this object can not be found.")
}

#' Make groups based on metadata table
#'
#' @param metadata a data frame contains several columns for metadata.
#' @param group_column a field name.
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
make_groups.Metadata <- function(m, subset_value, group_column = "group", levels = c("PM", "FA")) {
  if (length(levels) != 2) {
    stop("The number of levels is not equal to 2, only support two groups")
  }

  if (is.null(m$filtered_metadata)) {
    warning("Use the full metadata table. if you want to use subset, you need to run select_samples firstly.")
    metadata <- m$metadata
  } else {
    metadata <- m$filtered_metadata
  }

  # Select a subset
  metadata <- metadata[which(metadata$subset == subset_value),]
  cols <- colnames(metadata)
  if (group_column %in% cols) {
    groups <- metadata[,group_column]
    if (all(levels %in% groups)) {
      if (length(unique(groups)) != 2) {
        stop("The number of groups is not equal to 2, only support two groups")
      }

      return(factor(groups, levels = levels))
    } else {
      stop(paste0("No such group values: ", paste0(levels, collapse = ' or ')))
    }
  } else {
    stop(paste0("No such group_column: ", group_column))
  }
}

select_samples <- function(object, ...) {
  UseMethod("select_samples", object)
}

select_samples.default <- function(object, ...) {
  print("The class of this object cannot be found.")
}

#' Select samples
#'
#' @param metadata a data frame.
#' @param sql_query_str a query string.
#' @param subset_fields which fields do you want to filter.
#' @importFrom sqldf sqldf
#' @return a data frame
#' @export
select_samples.Metadata <- function(m, subset_fields, sql_query_str = NULL) {
  metadata <- m$metadata

  if (!(c("sample_id") %in% colnames(metadata))) {
    stop("No such column: sample_id")
  }

  if (length(unique(metadata$sample_id)) != length(metadata$sample_id)) {
    stop("sample_id need to be unique.")
  }

  if (!is.vector(subset_fields)) {
    stop("subset_fields must be a vector.")
  } else {
    if (length(subset_fields) == 0) {
      stop("subset_fields must have one element at least.")
    }
  }

  if (!all(subset_fields %in% colnames(metadata))) {
    stop(paste0("No such columns: ", paste0(subset_fields, collapse = " or ")))
  }

  subsets <- apply(subset(metadata, select=subset_fields), 1, paste0, collapse = "_")
  m$metadata$subset <- subsets
  m$subsets <- unique(subsets)
  m$subset_fields <- subset_fields

  if (is.null(sql_query_str)) {
    m$filtered_metadata <- m$metadata
  } else {
    metadata <- m$metadata
    query_str <- paste0(c("SELECT * FROM metadata WHERE", sql_query_str), collapse = " ")
    m$filtered_metadata <- sqldf(query_str)
  }

  print(paste0(c(length(unique(m$filtered_metadata$subset)), 
               "subsets are selected by", 
               paste0(subset_fields, collapse = " & ")),
               collapse = " "))
  return(m$filtered_metadata)
}
