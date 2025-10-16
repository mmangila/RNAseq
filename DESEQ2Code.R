find_de_deseq <- function(dge_deseq,
                          keyfile,
                          group,
                          batch_design,
                          padj,
                          paths,
                          fc_shrink,
                          surrogate_variable) {

  print("Preparing DESEQ2 data.")
  dge_deseq$samples$lib.size <- colSums(
    dge_deseq$counts
  )

  dds <- eval(parse(text = paste0(
    "DESeqDataSetFromMatrix(",
    "countData = dge_deseq$counts,",
    "colData   = keyfile,",
    "design    = ~ 1 + ",
    batch_design, ")"
  )))

  print("Estimating size factors")

  dds <- estimateSizeFactors(dds)

  if (surrogate_variable) {
    print("Running SVAseq")
    print(paste0(
      "model.matrix(~ 1 + ",
      batch_design,
      ", colData(dds))"
    ))

    dat <- counts(dds, normalized = TRUE)
    mod <- eval(parse(text = paste0(
      "model.matrix(~ 1 + ",
      batch_design,
      ", colData(dds))"
    )))
    mod0 <- model.matrix(~ 1, colData(dds))

    svs <- sva::num.sv(dat, mod, method = "leek")
    svseq <- run_svaseq(dat, mod, mod0, svs)

    if(svseq == 0) {
      print("Run again without SVAseq")
      find_de_deseq(dge_deseq,
                    keyfile,
                    group,
                    batch_design,
                    padj,
                    paths,
                    fc_shrink,
                    surrogate_variable = FALSE)
    }

    ddssva <- dds

    sapply(seq_along(svseq[1, ]), function(sv) {
      eval(parse(text = paste0(
        "ddsva$SV",
        sv,
        " <- svseq$sv[, ]",
        sv
      )))
    })

    design(ddssva) <- eval(parse(text = paste("~",
                                              paste0("SV", 1:(svs - 1),
                                                     collapse = " + "),
                                              "+", group)))
    dds <- ddssva
  }

  print("Run DESeq2")

  dds <- DESeq2::DESeq(dds)

  print(DESeq2::resultsNames(dds)) # lists the coefficients

  print("Begin PCA")
  de_deseq_pca(dds, group, paths)
  print("Creating gene tables")
  de_deseq_tables(keyfile, group, dds, padj, paths, fc_shrink)

}

de_deseq_tables <- function(keyfile, group, dds, padj, paths, fc_shrink) {

  combos <- eval(
    parse(
      text = paste0(
        "combn(as.data.frame(keyfile %>% distinct(",
        group,
        "))[,1],2)"
      )
    )
  )

  out_base <- paste0(paths[3],
                     "/DESEQ2/DE_tables/")
  ##
  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
  lfc_suffixes <- data.frame(
    Level = c(0, 1.5, 2),
    Suffix = c("_detags.csv",
               "_detags_1point5FC.csv",
               "_detags_2FC.csv")
  )

  print(lfc_suffixes)

  sapply(seq_along(combos[1, ]), function(x) {
    test_name <- paste0(
      combos[1, x],
      ".vs.",
      combos[2, x]
    )
    test_base_dir <- paste0(out_base,
                            test_name,
                            "/")
    dir.create(test_base_dir, showWarnings = FALSE)
    de_genes  <- results(
      dds,
      contrast = c(group,
                   as.character(combos[1, x]),
                   as.character(combos[2, x])),
      alpha = 0.99999
    )

    if (fc_shrink == TRUE) {
      de_genes <- lfc_shrink(dds,
                             contrast = c(group,
                                          as.character(combos[1, x]),
                                          as.character(combos[2, x])),
                             res      = de_genes,
                             type     = "ashr")
    }
    de_genes$X <- rownames(de_genes)

    write_csv(as.data.frame(de_genes),
              file = paste0(test_base_dir,
                            test_name,
                            "_alltags.csv"))
    print("Creating volcano plot")
    pdf(paste0(test_base_dir,test_name,"_volcano.pdf"), width = 5, height = 3.5)
    print(EnhancedVolcano(de_genes,
                          lab = rownames(de_genes),
                          x = "log2FoldChange",
                          y = "pvalue"))
    dev.off()

    sapply(seq_along(lfc_suffixes$Level), function(x) {
      print(print("Get significantly DE genes above", lfc_suffixes$Level[x], "FC"))
      de_genes_sigs <- filter_de_set(
        de_genes,
        lfc_suffixes$Level[x],
        padj
      )
      write_csv(as.data.frame(de_genes_sigs),
                file = paste0(test_base_dir,
                              test_name,
                              lfc_suffixes$Suffix[x]))
    })
  })
}

de_deseq_pca <- function(dds, group, paths) {
  vsd <- rlog(dds, blind = FALSE)
  print(head(assay(vsd), 3))
  pca_data <- plotPCA(vsd, intgroup = c(group), returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percent_var"))
  print(percent_var)

  pdf(paste0(paths[3], "/MDS/deseq2_e-counts_PCA_all_labels.pdf"),
      width = 5, height = 3.5)
  print(ggplot(pca_data,
               aes(PC1, PC2, color = group)) +
    geom_point(size = 3) +
    xlab(
      paste0("PC1: ", percent_var[1], "% variance")
    ) +
    ylab(
      paste0("PC2: ",percent_var[2],"% variance")
    ) +
    geom_text_repel(aes(label = name)))
  dev.off()

  pdf(paste0(paths[3], "/MDS/deseq2_e-counts_PCA.pdf"),
      width = 5, height = 3.5)
  print(ggplot(pca_data,
               aes(PC1, PC2, color = group)) +
    geom_point(size = 3) +
    xlab(
      paste0("PC1: ", percent_var[1], "% variance")
    ) +
    ylab(
      paste0("PC2: ", percent_var[2], "% variance")
    ))
  dev.off()
}


filter_de_set <- function(deset, lfc = 0, padj = 0.05) {
  filtered_de_set <- deset[
    which(
      deset$log2FoldChange > lfc &
        deset$padj < padj
    ),
  ]

  return(filtered_de_set)
}

run_svaseq <- function (dat, mod, mod0, svs) {
  if(svs == 0) {
    t <- 0
  } else {
    t <- try(sva::svaseq(dat, mod, mod0, n.sv = svs))
  }
  if("try-error" %in% class(t)) t <- run_svaseq(dat, mod, mod0, svs - 1)
  return(t)
}