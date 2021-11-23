find_de_deseq <- function (dge.DESeq, keyfile, group, padj, paths) {

  print("Preparing DESEQ2 data.")
  dge.DESeq$samples$lib.size <- colSums(
    dge.DESeq$counts
  )

  dds <- eval(parse(text = paste0(
    "DESeqDataSetFromMatrix(",
    "countData = dge.DESeq$counts,",
    "colData = keyfile,",
    "design= ~ ",
    group,")"
  )))

  dds <- DESeq(dds)

  resultsNames(dds) # lists the coefficients

  print("Begin PCA")
  de_deseq_pca(dds, group, paths)
  print("Creating gene tables.")
  de_deseq_tables(keyfile, group, dds, padj, paths)

}

de_deseq_tables <- function (keyfile, group, dds, padj, paths) {

  combos <- eval(
    parse(
      text = paste0(
        "combn(as.data.frame(keyfile %>% distinct(",
        group,
        "))[,1],2)"
      )
    )
  )

  out.base <- paste0(paths[3],
                     "/DESEQ2/DE_tables/"
  )
  ##
  dir.create(out.base,
             recursive=T,
             showWarnings = F
  )
  lfc.suffixes <- data.frame(
    Level = c(0,1.5,2),
    Suffix = c("_detags.csv",
               "_detags_1point5FC.csv",
               "_detags_2FC.csv"
    )
  )

  sapply(1:length(combos[1,]), function (x) {
    test.name <- paste0(
      combos[1,x],
      ".vs.",
      combos[2,x]
    )
    test.base.dir <- paste0(out.base,
                            test.name,
                            "/")
    dir.create(test.base.dir, showWarnings=F)
    de.genes <- lfcShrink(dds,
                        coef=paste0(group,
                                    "_",
                                    as.character(combos[1,x]),
                                    "_",
                                    as.character(combos[2,x])),
                        type="apeglm")
    de.genes  <- results(
      dds,
      contrast = c(group,
                   as.character(combos[1,x]),
                   as.character(combos[2,x])
      ),
      alpha = 0.99999
    )
    de.genes <- lfcShrink(dds,
                          contrast = c(group,
                                       as.character(combos[1,x]),
                                       as.character(combos[2,x])
                          ),
                          res = de.genes,
                          type="apeglm")
    de.genes$X <- rownames(de.genes)

    write_csv(as.data.frame(de.genes),
              file = paste0(test.base.dir,
                            test.name,
                            "_alltags.csv"
              )
    )

    sapply(1:3, function (x) {
      de.genes.sigs <- filter.de.set(
        de.genes,
        lfc.suffixes[x,1],
        padj
      )
      write_csv(as.data.frame(de.genes.sigs),
                file = paste0(test.base.dir,
                              test.name,
                              lfc.suffixes[x,2]))
    })
  })
}

de_deseq_pca <- function (dds, group, paths) {
  vsd <- rlog(dds, blind=FALSE)
  print(head(assay(vsd), 3))
  pcaData <- plotPCA(vsd, intgroup=c(group), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  print(percentVar)

  pdf(paste0(paths[3],"/MDS/deseq2_e-counts_PCA_all_labels.pdf"), width = 5, height = 3.5)
  print(ggplot(
    pcaData,
    aes(PC1, PC2, color=group)
    ) +
    geom_point(size=3) +
    xlab(
      paste0("PC1: ",percentVar[1],"% variance")
      ) +
    ylab(
      paste0("PC2: ",percentVar[2],"% variance")
      ) +
    geom_text_repel(aes(label = name)))
  dev.off()

  pdf(paste0(paths[3],"/MDS/deseq2_e-counts_PCA.pdf"), width = 5, height = 3.5)
  print(ggplot(
    pcaData,
    aes(PC1, PC2, color=group)
    ) +
    geom_point(size=3) +
    xlab(
      paste0("PC1: ",percentVar[1],"% variance")
      ) +
    ylab(
      paste0("PC2: ",percentVar[2],"% variance")
      ))
  dev.off()
}


filter.de.set <- function(deset, lfc = 0, padj = 0.05) {
  filtered.de.set <- deset[
    which(
      deset$log2FoldChange > lfc &
        deset$padj < padj
    ),
    ]

  return(filtered.de.set)
}
