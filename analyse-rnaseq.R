

analyse.RNAseq <- function (project.folder,analysis, group, func_path, script_path) {

    biocManagerLibs <- c("edgeR",
                         "limma",
                         "pcaMethods",
                         "DESeq2",
                         "biomaRt",
                         "Rgraphviz",
                         "topGO")
    libs            <- c("seqinr",
                         "car",
                         "data.table",
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
    source(paste0(script_path,"rewrittenHonoursCode.R"))
    source(paste0(script_path,"CombinedCode.R"))
    source(paste0(script_path,"DESEQ2Code.R"))
    source(paste0(script_path,"edgeRCode.R"))

    install_libraries(biocManagerLibs,BiocManager::install)
    install_libraries(libs,install.packages)

    project_paths <- file_paths(project.folder,analysis)
    keyfile <- create.folders(project_paths)
    assignment.summary(project_paths,keyfile)
    old.dge <- filter.wrapper(keyfile,group,project_paths)
    dge.DESeq <- old.dge

    #edgeRCode

    find_de_edger(old.dge, group, keyfile, project_paths)

    #DESEQ2Code

    find_de_deseq(dge.DESeq, keyfile, group, project_paths)

    # Find the union

    find_combined_de(keyfile, group, lfc.suffixes, func_path, project_paths)



}
