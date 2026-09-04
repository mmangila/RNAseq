lfc_suffixes <- data.frame(Level  =  c(0, 1.5, 2),
                           Suffix =  c("_detags.csv",
                                       "_detags_1point5FC.csv",
                                       "_detags_2FC.csv"))

analyse_rnaseq <-  function(project_folder,
                            analysis,
                            group,
                            padj,
                            mapman_focus,
                            annotation = FALSE, fc_shrink = FALSE,
                            surrogate_variable = FALSE,
                            batch_vars = vector(mode = "character"),
                            de_logfc = 1.5, go_pval = 0.05) {

  bioconductor_libs <-  c("edgeR",
                          "limma",
                          "pcaMethods",
                          "DESeq2",
                          "biomaRt",
                          "Rgraphviz",
                          "topGO",
                          "EnhancedVolcano",
                          "sva")
  base_libs         <-  c("seqinr",
                          "car",
                          "data.table",
                          "devtools",
                          "ggalluvial",
                          "ggplot2",
                          "ggrepel",
                          "ggthemes",
                          "GO.db",
                          "gplots",
                          "grid",
                          "lattice",
                          "parallel",
                          "pheatmap",
                          "plyr",
                          "RColorBrewer",
                          "reshape2",
                          "scales",
                          "scatterplot3d",
                          "statmod",
                          "tidytext",
                          "tidyverse",
                          "viridis",
                          "wesanderson",
                          "devtools",
                          "data.table")

  urls <- c(
    "https://github.com/mmangila/RNAseq/raw/main/CombinedCode.R",
    "https://github.com/mmangila/RNAseq/raw/main/DESEQ2Code.R",
    "https://github.com/mmangila/RNAseq/raw/main/edgeRCode.R",
    "https://github.com/mmangila/RNAseq/raw/main/Filtering.R",
    "https://github.com/mmangila/RNAseq/raw/main/MAplotGeneSetLimma.R"
  )

  sapply(urls, devtools::source_url)

  install_libraries(bioconductor_libs, BiocManager::install, "https://cran.csiro.au/")
  install_libraries(base_libs,         install.packages,     "https://cran.csiro.au/")

  project_paths <- file_paths(project_folder, analysis)
  keyfile       <- create_folders(project_paths)
  keyfile[sapply(keyfile, is.character)] <- lapply(
    keyfile[sapply(keyfile, is.character)],
    as.factor
  )
  keyfile[sapply(keyfile, is.numeric)] <- lapply(
    keyfile[sapply(keyfile, is.numeric)],
    as.factor
  )
  colnames(keyfile)[1] <- "Sample_ID"

  keyfile <- as.data.frame(keyfile)

  keyfile$sample_group <- keyfile[, group]
  keyfile <- as_tibble(keyfile)
  combos  <- combn(as.data.frame(keyfile %>% distinct(sample_group))[,1], 2)

  if (annotation) {
    func_path <- readline(paste(
      "Enter the location of the genome functional annotation here",
      "(Valid format: .csv): "
    ))
    funcs <- read_files(func_path)
    func_focus <- readline(paste(
      "Which column of the functional annotation",
      "is reflected in the featureCounts results? "
    ))
  } else {
    funcs <- NULL
    func_focus <- "X"
  }

  print("Begin analysis")
  assignment_summary(project_paths, keyfile)

  print("Filter data")
  old_dge   <- filter_wrapper(keyfile,
                              group,
                              project_paths,
                              annotation,
                              funcs, func_focus)
  dge_deseq <- read_data(keyfile, group, project_paths)

  batch_design <- create_design(group, batch_vars)

  # edgeRCode
  print("Running edgeR")
  find_de_edger(
    old_dge,
    group, batch_design, keyfile,
    project_paths, analysis,
    surrogate_variable, padj
  )

  # DESEQ2Code
  print("Running DESeq2")
  find_de_deseq(dge_deseq,
                keyfile,
                group,
                batch_design,
                padj,
                project_paths,
                fc_shrink,
                surrogate_variable)
  # Find the union
  print("Combining the analyses")

  find_combined_de(keyfile,
                   group,
                   combos,
                   lfc_suffixes,
                   padj,
                   funcs,
                   func_focus,
                   project_paths,
                   go,
                   project_folder,
                   de_logfc,
                   go_pval)

  # Find GO terms if they exist

  if (!is.null(funcs)) {
    print("Run GO analysis?")
    go_choice <- readline("Choose Yes/No:")

    if (tolower(go_choice) %in% c("yes", "yeah", "ye", "yea", "y", "agree")) {
      print("Begin GO analysis")
      devtools::source_url(
        "https://github.com/mmangila/RNAseq/raw/main/topGO_functions.R"
      )

      de_logfc <- readline("Input minimum gene logFC for GO analysis: ")
      go_pval  <- readline("Input GO term p-value threshold: ") 

      
      analyse_go(funcs, func_focus,
                 project_folder, combos, project_paths,
                 de_logfc, go_pval)
    }
  }

  print("Analysis finished")
}

create_design <- function (group, batch_vars) {
  all_vars <- c(group, batch_vars)
  
  all_design_terms <- unlist(lapply(seq_along(all_vars), function (size) {
    subsets <- combn(all_vars, size, simplify = FALSE)
    design_terms <- sapply(subsets, function (x) {paste(x, collapse = ":")})
    return(design_terms)
  }))

  batch_design <- paste0("~0+", paste(all_design_terms, collapse = "+"))
  return(as.formula(batch_design))
}

read_files <- function (file_to_read) {
  if (grepl(".tsv", file_to_read)) {
      funcs <- read.table(file_to_read, sep = "\t", header = TRUE)
    } else if (grepl(".csv", file_to_read)) {
      funcs <- read.csv(file = file_to_read)
    } else {
      warning("File format not recognised. Continuing without annotation.")
    }
}
