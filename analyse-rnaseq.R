require("ggplot2")
library(ggplot2)
options(ggrepel.max.overlaps = Inf)

text_size_theme_8 <- theme(axis.text=element_text(size=8),
    axis.title=element_text(size=8),
    axis.text.x=element_text(angle = 45, hjust = 1),
    legend.title=element_text(size=8),
    legend.text=element_text(size=8))

# text sizes
text_size_theme_8_labels <- theme(axis.text=element_text(size=8),
       axis.title=element_text(size=8),
       axis.text.x=element_text(angle = 45, hjust = 1),
       legend.title=element_text(size=8),
       legend.text=element_text(size=8),
       panel.background = element_rect(fill = "transparent") # bg of the panel
   , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
   , panel.grid.major = element_blank() # get rid of major grid
   , panel.grid.minor = element_blank() # get rid of minor grid
   , legend.background = element_rect(fill = "transparent") # get rid of legend bg
   , legend.box.background = element_rect(fill = "transparent")) # get rid of legend panel bg
# magic geom_text conversion ratio
# https://stackoverflow.com/questions/17311917/ggplot2-the-unit-of-size
label_size = 25.4/72 * 8

lfc.suffixes <- data.frame(
  Level = c(0,1.5,2),
  Suffix = c("_detags.csv",
             "_detags_1point5FC.csv",
             "_detags_2FC.csv"
  )
)

analyse.RNAseq <- function (project.folder,
                            analysis,
                            group,
                            padj,
                            mapman_focus,
                            annotation = FALSE, go = FALSE, fcShrink = FALSE) {

    biocManagerLibs <- c("edgeR",
                         "limma",
                         "pcaMethods",
                         "DESeq2",
                         "biomaRt",
                         "Rgraphviz",
                         "topGO",
                         "EnhancedVolcano")
    libs            <- c("seqinr",
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

    devtools::source_url(
      "https://github.com/mmangila/RNAseq/raw/main/CombinedCode.R"
    )
    devtools::source_url(
      "https://github.com/mmangila/RNAseq/raw/main/DESEQ2Code.R"
      )
    devtools::source_url(
      "https://github.com/mmangila/RNAseq/raw/main/edgeRCode.R"
      )
    devtools::source_url(
      "https://github.com/mmangila/RNAseq/raw/main/rewrittenHonoursCode.R"
      )
    devtools::source_url(
      "https://github.com/pedrocrisp/NGS-pipelines/raw/master/R_functions/MAplotGeneSetLimma.R"
      )

    install_libraries(biocManagerLibs,BiocManager::install)
    install_libraries(libs,install.packages)

    project_paths <- file_paths(project.folder,analysis)
    keyfile <- create.folders(project_paths)
    keyfile[sapply(keyfile, is.character)] <- lapply(keyfile[sapply(keyfile, is.character)],
                                       as.factor)
    keyfile[sapply(keyfile, is.numeric)] <- lapply(keyfile[sapply(keyfile, is.numeric)],
                                       as.factor)
    colnames(keyfile)[1] <- "Sample_ID"

    if (annotation) {
      func_path <- readline("Enter the location of the genome functional annotation here (Valid format: .csv): ")
      func_focus <- readline("Which column of the functional annotation is reflected in the featureCounts results? ")
    } else {
      func_path <- NULL
      func_focus <- "X"
    }

    print("Begin analysis")
    assignment.summary(project_paths,keyfile)

    print("Filter data")
    old.dge <- filter.wrapper(keyfile,group,project_paths)
    dge.DESeq <- read.data(keyfile,group,project_paths)

    #edgeRCode
    print("Running edgeR")
    find_de_edger(old.dge, group, keyfile, project_paths, analysis, padj)

    #DESEQ2Code
    print("Running DESeq2")
    find_de_deseq(dge.DESeq, keyfile, group, padj, project_paths, fcShrink)

    # Find the union
    print("Combining the analyses")
    find_combined_de(keyfile,
                     group,
                     lfc.suffixes,
                     annotation,
                     func_path,
                     func_focus,
                     project_paths,
                     go,
                     project.folder)
    
    print("Generate MapMan files")
    dir.create(paste0(project_paths[3], "/Combined/MapMan"))
    sapply(1:length(combos[1, ]), function (x) {
      test.name <- paste0(as.character(combos[1, x]),
                          ".vs.",
                          as.character(combos[2, x]))
      comb_res <- read.csv(paste0(project_paths[3],
                                  "/Combined/DE_tables/",
                                  test.name, "/",
                                  test.name,
                                  "_detags_1point5FC.csv"))
      write.table(comb_res[, which(colnames(comb_res) %in% c(mapman_focus, "edger_logFC"))],
                  file = paste0(project_paths[3],
                                "/Combined/MapMan/",
                                test.name,
                                "_detags_1point5FC_mapman.txt"),
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE)
})
    
    print("Analysis finished")



}
