find_de_deseq <- function (dge.DESeq, keyfile, group, paths) {

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



  print("Creating gene tables.")
  de_deseq_tables(keyfile, group, dds, paths)

}

de_deseq_tables <- function (keyfile, group, dds, paths) {

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
    de.genes  <- results(
      dds,
      contrast = c(group,
                   as.character(combos[1,x]),
                   as.character(combos[2,x])
      ),
      alpha = 0.99999
    )
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
        0.05
      )
      write_csv(as.data.frame(de.genes.sigs),
                file = paste0(test.base.dir,
                              test.name,
                              lfc.suffixes[x,2]))
    })
  })
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
