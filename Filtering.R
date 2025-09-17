

create_folders <- function(project_paths) {


  lapply(project_paths, function (x) dir.create(x, showWarnings = FALSE))
  dir.create(
    paste0(project_paths[3], "/MDS"),
    showWarnings = FALSE
  )

  keyfile <- read_csv(project_paths[4])
  
  return(keyfile)
}

file_paths <- function(project_folder, analysis) {
  paths <- vector(mode = "character", length = 0)

  paths <- c(paths, paste0(
    project_folder,
    "/",
    analysis
  ))

  time_stamp <- format(
    Sys.time(),
    "%Y-%m-%d-%H%M"
  )

  paths <- c(paths, readline(
    "Enter the location of your featureCounts results here: "
  ))

  paths <- c(paths, paste0(
    paths[1],
    "/DE_analysis_",
    time_stamp
  ))

  paths <- c(paths, readline("Enter the location of the keyfile here: "))
  paths <- c(paths, paste0(paths[1], "/Mapping_stats"))

  return(paths)
}

## Summary of reads assigned to features by featureCounts (after mapping)

install_libraries <- function(libraries, installer) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  lapply(libraries, function(x) {
    if (!require(x, character.only = TRUE)) do.call(installer, list(x))
    library(x, character.only = TRUE)
  })
}
assignment_summary <- function(paths, keyfile) {

  print("Reading in files")
  files <- dir(
    paste0(paths[2], "/"),
    pattern = "*.counts.summary"
  )
  count_files <- paste0(paths[2], "/", files)
  fcs         <- do.call(
    "cbind",
    lapply(
      count_files,
      function(fn) data.frame(read.delim(fn, row.names = 1))
    )
  )

  print("Creating mapping plots")

  fcs            <- data.frame(t(fcs))
  row.names(fcs) <- as.data.frame(keyfile)[, 1]
  fcs_t          <- t(fcs)
  fcs_plot       <- melt(fcs_t)

  prop           <- prop.table(fcs_t, margin = 2)
  prop_t         <- t(prop)
  prop_plot      <- melt(prop_t)

  plot_data      <- prop_plot

  plot_colours   <- rev(
    c("black", "grey",
      viridis_pal(option = "C")(n = length(levels(plot_data$Var2))))
  )

  print("Proportional mapped plot")

  pdf(paste0(paths[5], "/percent_mapped_plot.pdf"), w = 8, h = 4)
  print(
    ggplot(plot_data, aes(x = Var1, y = rev(value), fill = rev(Var2))) +
      geom_bar(stat = "identity")                                      +
      xlab("Sample")                                                   +
      ylab("Percentage of reads")                                      +
      scale_fill_manual(values = plot_colours, name = "Mapping cat")   +
      text_size_theme_8_labels                                         +
      theme_classic()                                                  +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  dev.off()

  ###Assigned reads as total read number

  print("Total mapped plot")

  pdf(paste0(paths[5], "/total_mapped_plot.pdf"), w = 6, h = 4)
  plot_data <- fcs_plot
  print(
    ggplot(plot_data, aes(x = Var2, y = value, fill = Var1))         +
      geom_bar(stat = "identity")                                    +
      xlab("Sample")                                                 +
      ylab("Number of reads")                                        +
      scale_fill_manual(values = plot_colours, name = "Mapping cat") +
      theme_classic()                                                +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
}

read_data <- function(keyfile, group, paths) {
  ## Read in data
  keyfile <- as.data.frame(keyfile)
  samples <- keyfile[, 1]
  sample_groups <- keyfile[, which(colnames(keyfile) == group)]

  count_files <- paste0(paths[2], "/", samples, ".counts")

  tmp <- read.delim(count_files[[1]], header = TRUE, skip = 1)
  gene_lengths        <- as.numeric(tmp[, "Length"])
  names(gene_lengths) <- tmp[, "Geneid"]
  rm(tmp)

  dge <- readDGE(count_files, skip = 1,
                 columns = c(1, 7),
                 group = sample_groups,
                 labels = as.character(samples))
  dge_raw <- dge
  gene_names <- as.character(rownames(dge$counts))
  gene_names_raw <- as.character(rownames(dge_raw$counts))

  return(dge)
}

filter.wrapper <-  function(keyfile,
                            group,
                            paths,
                            annotation,
                            func_path,
                            func_focus) {
  dge        <- read_data(keyfile, group, paths)
  gene_names <- as.character(rownames(dge$counts))

  dge_raw        <- dge
  gene_names_raw <- as.character(rownames(dge_raw$counts))
  fcs_folder <- paste0(paths[3], "/GEO_featureCountsSummaries")
  dir.create(fcs_folder, showWarnings = FALSE)

  raw_counts <- data.frame(dge$counts)
  raw_counts$AGI <- row.names(raw_counts)
  write.csv(raw_counts,
            paste0(fcs_folder, "/raw_counts_matrix.csv"),
            row.names = FALSE)

  old_dge        <- dge
  old_gene_names <- gene_names

  dge <- filtering_step(raw_counts, old_dge)
  gene_names <- as.character(rownames(dge$counts))

  if (annotation) {
    if (grepl(".tsv", func_path)) {
      funcs <- read.table(func_path, sep = "\t", header = TRUE)
    } else if (grepl(".csv", func_path)) {
      funcs <- read.csv(file = func_path)
    } else {
      errorCondition("File format not recognised")
    }

    gene_match <- match(row.names(dge$counts),
                        funcs[, which(colnames(funcs) == func_focus)])
    gene_description_matched <- funcs[gene_match, ]

    dge$genes <- gene_description_matched
  }
  return(dge)
}

filtering_step <- function(raw_counts, dge) {

  sample_num <- seq_along((dge$samples[, 1]))
  boxplot(colSums(raw_counts[, sample_num]))

  #42604 genes have at least 1 count
  print("Table of genes with at least one count:")
  print(table(rowSums(raw_counts[, sample_num]) > 0))

  #39728 genes have at least 1 count in 2 samples
  print("Table of genes with at least 1 count in 2 samples")
  print(table(rowSums((raw_counts[, sample_num]) > 0) > 2))

  #histogram of number of samples having 1 or more counts per gene
  hist(rowSums((raw_counts[, sample_num]) > 0), breaks = 100)

  #histogram of reads per gene
  hist(log10(rowSums(raw_counts[, sample_num])), breaks = 100)

  min_reads <- readline("Set number of minimum Count per Million: ")
  min_samples <- readline("Set number of minimum samples with minimum reads: ")

  min_reads <- as.numeric(min_reads)
  min_samples <- as.numeric(min_samples)

  loci_to_keep <- rowSums(cpm(dge) > min_reads) > min_samples

  print("Table of genes that passed or failed the filter")
  print(table(loci_to_keep))
  hist(
    log10(
      rowSums(
        raw_counts[
          which(rownames(raw_counts) %in%
                  names(which(loci_to_keep == TRUE))),
          sample_num
        ]
      )
    ),
    breaks = 100
  )

  confirm      <- readline("Accept these filters? y/n ")
  filtered_dge <- continue_filter(confirm, raw_counts, dge, loci_to_keep)
  return(filtered_dge)

}

continue_filter <- function(confirm, raw_counts, dge, loci_to_keep) {
  if (confirm == "y") {
    old_dge <- dge
    dge <- old_dge[loci_to_keep, ]
    dge$samples$lib.size <- colSums(dge$counts)

    #make sure gene_names get updated

    return(dge)
  } else if (confirm == "n") {
    filtering_step(raw_counts, dge)
  } else {
    new_confirm <- readline("Invalid response. Accept these filters? y/n ")
    continue_filter(new_confirm, raw_counts, dge)
  }
}
