require("ggplot2")
library(ggplot2)
options(ggrepel.max.overlaps = Inf)

text_size_theme_8 <- theme(axis.text    = element_text(size  = 8),
                           axis.title   = element_text(size  = 8),
                           legend.title = element_text(size  = 8),
                           legend.text  = element_text(size  = 8),
                           axis.text.x  = element_text(angle = 45, hjust = 1))

transparent_element <- element_rect(fill = "transparent")

# text sizes
text_size_theme_8_labels <- theme(
  axis.text    = element_text(size = 8),
  axis.title   = element_text(size = 8),
  legend.title = element_text(size = 8),
  legend.text  = element_text(size = 8),
  axis.text.x  = element_text(angle = 45,
                              hjust = 1),

  panel.grid.major      = element_blank(),
  panel.grid.minor      = element_blank(),
  legend.background     = transparent_element,
  legend.box.background = transparent_element,
  panel.background      = transparent_element,
  plot.background       = element_rect(fill = "transparent",
                                       color = NA)
)

# magic geom_text conversion ratio
# https://stackoverflow.com/questions/17311917/ggplot2-the-unit-of-size
label_size <- 25.4 / 72 * 8

lfc_suffixes <- data.frame(
  Level = c(0, 1.5, 2),
  Suffix = c("_detags.csv",
             "_detags_1point5FC.csv",
             "_detags_2FC.csv")
)

analyse_rnaseq <-  function(project_folder,
                            analysis,
                            group,
                            padj,
                            mapman_focus,
                            annotation = FALSE, go = FALSE, fc_shrink = FALSE,
                            surrogate_variable = FALSE) {

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
    "https://github.com/mmangila/RNAseq/raw/main/rewrittenHonoursCode.R",
    "https://github.com/pedrocrisp/NGS-pipelines/raw/master/R_functions/MAplotGeneSetLimma.R"
  )

  sapply(urls, devtools::source_url)

  install_libraries(bioconductor_libs, BiocManager::install)
  install_libraries(base_libs,         install.packages)

  project_paths <- file_paths(project_folder, analysis)
  keyfile       <- create.folders(project_paths)
  keyfile[sapply(keyfile, is.character)] <- lapply(
    keyfile[sapply(keyfile, is.character)],
    as.factor
  )
  keyfile[sapply(keyfile, is.numeric)] <- lapply(
    keyfile[sapply(keyfile, is.numeric)],
    as.factor
  )
  colnames(keyfile)[1] <- "Sample_ID"

  if (annotation) {
    func_path <- readline(
      "Enter the location of the genome functional annotation here (Valid format: .csv): "
    )
    func_focus <- readline(
      "Which column of the functional annotation is reflected in the featureCounts results? "
    )
  } else {
    func_path <- "None"
    func_focus <- "X"
  }

  print("Begin analysis")
  assignment.summary(project_paths, keyfile)

  print("Filter data")
  old_dge   <- filter.wrapper(keyfile, group, project_paths, annotation, func_path, func_focus)
  dge_deseq <- read.data(keyfile, group, project_paths)

  # DESEQ2Code
  print("Running DESeq2")
  find_de_deseq(
    dge_deseq, keyfile, group,
    padj, project_paths, fc_shrink, surrogate_variable
  )

  #edgeRCode
  print("Running edgeR")
  find_de_edger(old_dge, group, keyfile, project_paths, analysis,
                surrogate_variable, padj)



  # Find the union
  print("Combining the analyses")
  find_combined_de(keyfile,
                   group,
                   lfc_suffixes,
                   annotation,
                   func_path,
                   func_focus,
                   project_paths,
                   go)

  print("Analysis finished")



}
