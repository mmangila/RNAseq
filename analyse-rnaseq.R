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

lfc.suffixes <- data.frame(
 Level = c(0,1.5,2),
 Suffix = c("_detags.csv",
            "_detags_1point5FC.csv",
            "_detags_2FC.csv"
 )
)

source("rewrittenHonoursCode.R")
source("CombinedCode.R")
source("DESEQ2Code.R")
source("edgeRCode.R")

analyse.RNAseq <- function (project.folder,analysis, group, func_path) {

    install_libraries(biocManagerLibs,BiocManager::install)
    install_libraries(libs,install.packages)

    project_paths <- file_paths(project.folder,analysis,func_path)
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
