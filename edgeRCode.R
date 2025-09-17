find_de_edger <- function(old_dge,
                          group,
                          batch_design,
                          keyfile,
                          paths,
                          analysis,
                          surrogate_variable,
                          padj) {

  print("Calculating normalisation factors")
  dge <- edgeR::calcNormFactors(old_dge, method = "TMM")

  print(batch_design)
  design <- eval(parse(text = paste0(
    "model.matrix(~ 0 + ",
    batch_design,
    ", data = keyfile)"
  )))

  group_levels <- eval(parse(text = paste0("levels(keyfile$)", group)))
  print(group_levels)

  colnames(design)[seq_along(group_levels)] <- eval(parse(
    text = paste0(
      "levels(as.factor(keyfile$",
      group,
      "))"
    )
  ))

  if (surrogate_variable) {
    print("Running SV analysis")
    dat    <- cpm(dge)
    mod    <- eval(parse(text = paste0(
      "model.matrix(~ ",
      batch_design,
      ", data = keyfile)"
    )))
    mod0   <- model.matrix(~ 1, data = keyfile)
    svseq  <- sva::svaseq(dat, mod, mod0, n.sv = 2)
    colnames(svseq$sv) <- c("SV1", "SV2")
    design <- cbind(design, svseq$sv)
  }

  print("Running limma and voom")
  v   <- edgeR::voomLmFit(dge,
                          design,
                          plot = TRUE, sample.weights = surrogate_variable)
  fit <- limma::eBayes(v)

  # print("Begin MDS")
  # de_edger_mds(dge, group, keyfile, paths, v)

  print("Begin PCA")
  de_edger_pca(v, group, keyfile, paths)
  print("Begin DE gene table generation")
  gene_names <- as.character(rownames(dge$counts))
  de_edger_tables(keyfile,
                  group,
                  fit,
                  paths,
                  design,
                  gene_names,
                  analysis,
                  v,
                  padj)
}

de_edger_mds <- function(dge, group, keyfile, paths, v) {
  sample_colours <- ggthemes::ptol_pal()(length(dge$samples$lib.size) / 2)
  if (length(dge$samples$lib.size) %% 2 == 1) {
    ddd_sample_colours <- c(rep(sample_colours, each = 2),
                            ggthemes::ptol_pal()(1))
  } else {
    ddd_sample_colours <- c(rep(sample_colours, each = 2))
  }

  mds2 <- plotmds_invisible(v, ndim = 3)

  # Make MDS plots
  par(mfrow = c(2, 2), ps = 10)
  sapply(1:3, function(x) {
    with(mds2,
         eval(parse(text = paste0("plot(V",
                                  combn(1:3, 2)[1, x],
                                  ", V",
                                  combn(1:3, 2)[2, x],
                                  ", pch=16,",
                                  " col=sample_colours, ",
                                  " cex=1,",
                                  "main='MDS plot',",
                                  "xlab='Leading logFC dim",
                                  combn(1:3, 2)[1, x],
                                  "',",
                                  "ylab='Leading logFC dim ",
                                  combn(1:3, 2)[2, x],
                                  "')"))))

  })

  mds2$group <- keyfile[, which(colnames(keyfile) == group)]

  with(mds2, {
    s3d <-  scatterplot3d(V1, V3, V2,
                          color       = ddd_sample_colours,
                          pch         = 19,
                          cex.symbols = 1.5,
                          type        = "h",
                          main        = "3D MDS plot",
                          xlab        = "Leading logFC dim 1",
                          ylab        = "Leading logFC dim 3",
                          zlab        = "Leading logFC dim 2")
    s3d_coords <- s3d$xyz.convert(V1, V3, V2)
    text(s3d_coords$x,
         s3d_coords$y,
         labels = row.names(mds2),
         cex = 0.5, pos = 4)
  })

  par(mfrow = c(1, 1))

  pdf(paste0(paths[3], "/MDS/voom_mds_3d.pdf"), width = 6.5, height = 6.5)
  with(mds2, {
    s3d        <- scatterplot3d(V1, V3, V2,
                                pch         = 19,
                                cex.axis    = 0.5,
                                cex.symbols = 1,
                                cex.lab     = 0.5,
                                type        = "h",
                                xlab        = "Leading logFC dim 1",
                                ylab        = "Leading logFC dim 3",
                                zlab        = "Leading logFC dim 2")
    s3d_coords <- s3d$xyz.convert(V1, V3, V2)
    text(s3d_coords$x,
         s3d_coords$y,
         labels = row.names(mds2),
         cex = 0.25, pos = 4)
  })
  dev.off()
}

de_edger_pca <- function(v, group, keyfile, paths) {
  voom_matrix <- v$EList$E

  cpm_tbl <- tibble::as_tibble(voom_matrix, rownames = "Gene")

  # replace NaN with NA? Nope there arent any...
  which(is.na(cpm_tbl))

  #mds_table <- limma::plotMDS(cpm_tbl)

  pdf(paste0(paths[3],"/MDS/voom_e-counts_PCA_Treatment.pdf"),
      width = 5, height = 3.5)
  make_PCA_plots(Timepoint = "ALL", dot_colour = group, keyfile, cpm_tbl)
  dev.off()

  pdf(paste0(paths[3],"/MDS/voom_e-counts_PCA_Treatment_labels.pdf"),
      width = 10, height = 7)
  make_PCA_plots_large_labels(
    dot_colour = group,
    cpm_tbl,
    keyfile
  )
  dev.off()
}

de_edger_tables <- function(keyfile,
                            group,
                            fit,
                            paths,
                            design,
                            gene_names,
                            analysis,
                            v,
                            padj) {
  combos <- eval(
    parse(
      text = paste0(
        "combn(as.data.frame(keyfile %>% distinct(",
        group,
        "))[,1],2)"
      )
    )
  )

  inside <- vector(mode = "character")
  inside2 <- sapply(seq_along(combos[1, ]), function(x) {
    inside <- c(
      inside,
      paste0(
        combos[1, x],
        ".vs.",
        combos[2, x],
        " = ",
        combos[1, x],
        " - ",
        combos[2, x]
      )
    )
  })

  if (length(inside2) > 1) {
    contrast_matrix <- eval(
      parse(
        text = paste0("makeContrasts(",
                      paste(inside2, collapse = ","),
                      ", levels = design)")
      )
    )
  } else if (length(inside2) == 1) {
    contrast_matrix <- eval(
      parse(
        text = paste0(
          "makeContrasts(",
          inside2,
          ", levels = design)"
        )
      )
    )
  }

  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2)

  results            <- limma::decideTests(fit2)
  results_2fc        <- limma::decideTests(fit2, lfc = log2(2))
  results_1point5_fc <- limma::decideTests(fit2, lfc = log2(1.5))

  print("Differentially expressed genes")
  print(summary(results))
  print("Differentially expressed genes with >1.5 fold change")
  print(summary(results_1point5_fc))
  print("Differentially expressed genes with >2 fold change")
  print(summary(results_2fc))

  de_table_all <- limma::topTable(fit = fit2, number = Inf)
  de_table_2fc <- limma::topTable(fit = fit2, number = Inf, lfc = 1)

  coefs <- colnames(fit2$coefficients)
  out_base <- paste0(paths[3], "/edgeR/DE_tables/")
  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

  for (comparison in coefs){
    test_name <- comparison
    test_base_dir <- paste0(out_base, test_name, "/")
    dir.create(test_base_dir, showWarnings = FALSE)
    print(test_name)

    pdf(paste0(test_base_dir, test_name, "_volcano.pdf"))
    plot_table <-  topTable(fit     = fit2,
                            number  = Inf, coef = comparison,
                            sort.by = "logFC")
    with(plot_table, plot(logFC, -(log(adj.P.Val)), pch = 16, cex = 0.3))
    title(main = paste0("Foldchange vs FDR (adjusted p-value)",
                        comparison), cex.main = 0.5)
    plot_table <-  topTable(fit     = fit2,
                            number  = Inf, coef = comparison,
                            sort.by = "logFC",
                            p.value = padj, lfc = log2(2))
    # writ e if statement to accou nt for when there are no DE genes
    if (length(plot_table) > 0) {
      with(plot_table, points(logFC, -(log(adj.P.Val)), pch = 16,
                              col = "red", cex = 0.5))
      dev.off()
    } else {
      dev.off()
    }

    results            <- decideTests(fit2)
    results_2fc        <- decideTests(fit2, lfc = log2(2))
    results_1point5_fc <- decideTests(fit2, lfc = log2(1.5))

    pdf(paste0(test_base_dir, test_name, "_smear.pdf"))
    MAplotGeneSetLimma(
      MArrayLMobject  = fit2,
      resultsMatrix   = results,
      geneSetList     = comparison,
      geneSetListName = comparison,
      inputList       = comparison,
      inputListName   = paste0("Up-regulated genes", comparison)
    )
    dev.off()

    pdf(paste0(test_base_dir, test_name, "_smear1point5FC.pdf"))
    MAplotGeneSetLimma(
      MArrayLMobject  = fit2,
      resultsMatrix   = results_1point5_fc,
      geneSetList     = comparison,
      geneSetListName = comparison,
      inputList       = comparison,
      inputListName   = paste0("Up-regulated genes (1.5FC)", comparison)
    )
    dev.off()

    pdf(paste0(test_base_dir, test_name, "_smear2FC.pdf"))
    MAplotGeneSetLimma(
      MArrayLMobject  = fit2,
      resultsMatrix   = results_2fc,
      geneSetList     = comparison,
      geneSetListName = comparison,
      inputList       = comparison,
      inputListName   = paste0("Up-regulated genes (2FC)", comparison)
    )
    dev.off()

    tt <- topTable(fit = fit2, number = Inf, coef = comparison)
    write.csv(tt, paste0(test_base_dir, test_name, "_alltags.csv"))
    tt <- topTable(fit = fit2, number = Inf, coef = comparison,
                   p.value = padj)
    write.csv(tt, paste0(test_base_dir, test_name, "_detags.csv"))
    tt <- topTable(fit = fit2, number = Inf, coef = comparison,
                   p.value = padj, lfc = log2(1.5))
    write.csv(tt, paste0(test_base_dir, test_name, "_detags_1point5FC.csv"))
    tt <- topTable(fit = fit2, number = Inf, coef = comparison,
                   p.value = padj, lfc = log2(2))
    write.csv(tt, paste0(test_base_dir, test_name, "_detags_2FC.csv"))
  }

  tests <- mclapply(coefs,
                    function(x) {
                      topTable(fit     = fit2,
                               number  = Inf,
                               coef    = x,
                               sort.by = "none")
                    })
  names(tests) <- coefs
  fc_matrix <- as.data.frame(sapply(tests,
                                    function(t) t$logFC, simplify = "array"))
  fdr_matrix <- sapply(tests, function(t) t$adj.P.Val, simplify = "array")
  #be careful that gene_names is in the correct order...
  rownames(fc_matrix)  <- gene_names
  rownames(fdr_matrix) <- gene_names
  write.csv(fc_matrix,  file = paste0(out_base, analysis, "_fc.csv"))
  write.csv(fdr_matrix, file = paste0(out_base, analysis, "_fdr.csv"))

  tests_sig_2fc <- mclapply(coefs,
                            function(x) {
                              topTable(fit     = fit2,
                                       p.value = padj, lfc = log2(2),
                                       number  = Inf, coef = x,
                                       sort.by = "none")
                            })
  names(tests_sig_2fc) <- coefs

  tests_sig_1point5fc <- mclapply(coefs,
                                  function(x) {
                                    topTable(fit     = fit2,
                                             p.value = padj, lfc = log2(1.5),
                                             number  = Inf, coef = x,
                                             sort.by = "none")
                                  })
  names(tests_sig_1point5fc) <- coefs

  cpm_matrix <- cpm(v)
  write.csv(fdr_matrix, file = paste0(out_base, analysis, "_fdr.csv"))
}

plotmds_invisible <- function(...) {
  ff <- tempfile()
  png(filename = ff)
  mds <- plotMDS(...)
  dev.off()
  unlink(ff)
  mds2 <- as.data.frame(mds$cmdscale.out)
  mds2
}


######### PCA 1 FUNCTION - JUST DOTS - SMALL PDF

make_PCA_plots <- function(Timepoint = "ALL", dot_colour, keyfile, cpm_tbl){
  # Timepoint = "ALL"
  # dot_colour = "Time"

  # If else to subset data by timepoint
  if (Timepoint == "ALL") {

    tc_pca <- pca(cpm_tbl[,-1],
                  scale = "uv",
                  center = T,
                  nPcs = 3,
                  method = "nipals")

  } else {

    cpm_table_group <- cpm_tbl %>%
      gather(key = Sample_ID, value = log2CPM, -Gene) %>%
      left_join(keyfile, by = "Sample_ID")

    cpm_table_group_spread <- cpm_table_group %>%
      select(Gene, Sample_ID, log2CPM) %>%
      spread(key = Sample_ID, value = log2CPM)

    pca_input_table <- cpm_table_group_spread[, -1]

    tc_pca <- pca(pca_input_table,
                  scale = "uv",
                  center = TRUE,
                  nPcs = 3,
                  method = "nipals")
  }

  # to get the variance (R2) explained by each axis
  summary(tc_pca)
  pc1_v <- round(pull(as_tibble(summary(tc_pca)), PC1)[1] * 100, digits = 1)
  pc2_v <- round(pull(as_tibble(summary(tc_pca)), PC2)[1] * 100, digits = 1)
  pc3_v <- round(pull(as_tibble(summary(tc_pca)), PC3)[1] * 100, digits = 1)

  pca_loadings <- as_tibble(loadings(tc_pca)) %>%
    mutate(Sample_ID = row.names(loadings(tc_pca))) %>%
    left_join(keyfile, by = "Sample_ID")

  ### plots PC1 vs PC2
  g <- ggplot(pca_loadings,
              aes(x = PC1, y = PC2,
                  colour = as.character(!!as.name(dot_colour)))) +
        geom_point() +
        scale_colour_ptol() +
        theme_classic() +
        labs(title = paste0("Timepoint ", Timepoint, " hr"),
             x     = paste0("PC1 (", pc1_v, "%)"),
             y     = paste0("PC2 (", pc2_v, "%)"))
  print(g)

  # add text labels
  g_labs <- g + geom_text_repel(aes(label = Sample_ID), show.legend = FALSE)

  g <- ggplot(pca_loadings,
              aes(x = PC1, y = PC2,
                  colour = as.character(!!as.name(dot_colour)))) +
         geom_point() +
         theme_classic() +
         scale_colour_ptol() +
         labs(title = paste0("Timepoint ", Timepoint, " hr - coord prop"),
              x = paste0("PC1 (", pc1_v, "%)"),
              y = paste0("PC2 (", pc2_v, "%)")) +
         coord_equal()
  print(g)

  ### plots PC1 vs PC3
  g <- ggplot(
    pca_loadings,
    aes(x = PC1, y = PC3, colour = as.character(!!as.name(dot_colour)))
    ) +
    geom_point() +
    scale_colour_ptol() +
    theme_classic() +
    labs(title = paste0("Timepoint ", Timepoint, " hr"),
         x = paste0("PC1 (", pc1_v, "%)"),
         y = paste0("PC3 (", pc3_v, "%)"))
  print(g)

  # add text labels
  g_labs <- g + geom_text_repel(aes(label=Sample_ID), show.legend = F)
  # print(g_labs)

  g <- ggplot(
    pca_loadings,
    aes(x = PC1, y = PC3, colour = as.character(!!as.name(dot_colour)))
    ) +
    geom_point() +
    theme_classic() +
    scale_colour_ptol() +
    labs(title = paste0("Timepoint ", Timepoint, " hr - coord prop"),
         x = paste0("PC1 (", pc1_v, "%)"),
         y = paste0("PC3 (", pc3_v, "%)")) +
    coord_equal()
  print(g)

  ### plots PC2 vs PC3
  g <- ggplot(pca_loadings, aes(x = PC2, y = PC3, colour = as.character(!!as.name(dot_colour)))) +
    geom_point() +
    scale_colour_ptol() +
    theme_classic() +
    labs(title = paste0("Timepoint ", Timepoint, " hr"),
         x = paste0("PC2 (", pc2_v, "%)"),
         y = paste0("PC3 (", pc3_v, "%)"))
  print(g)

  # add text labels
  g_labs <- g + geom_text_repel(aes(label=Sample_ID), show.legend = F)
  # print(g_labs)

  g <- ggplot(pca_loadings, aes(x = PC2, y = PC3, colour = as.character(!!as.name(dot_colour)))) +
    geom_point() +
    theme_classic() +
    scale_colour_ptol() +
    labs(title = paste0("Timepoint ", Timepoint, " hr - coord prop"),
         x = paste0("PC2 (", pc2_v, "%)"),
         y = paste0("PC3 (", pc3_v, "%)")) +
    coord_equal()
  print(g)
}


######### ######### ######### ######### #########
######### PCA 2 FUNCTION - labels - LARGE PDF

make_PCA_plots_large_labels <- function(Timepoint = "ALL", dot_colour, keyfile, cpm_tbl){

  if(Timepoint == "ALL"){

    tc_pca <- pca(cpm_tbl[,-1], scale = "uv", center = T, nPcs = 3, method = "nipals")

  }else{

    cpm_table_group <- cpm_tbl %>%
      gather(key = Sample_ID, value = log2CPM, -Gene) %>%
      left_join(keyfile, by = "Sample_ID") %>%
      filter(Time == Timepoint)

    cpm_table_group_spread <- cpm_table_group %>%
      select(Gene, Sample_ID, log2CPM) %>%
      spread(key = Sample_ID, value = log2CPM)

    pca_input_table <- cpm_table_group_spread[,-1]

    tc_pca <- pca(pca_input_table, scale = "uv", center = T, nPcs = 3, method = "nipals")
  }

  slplot(tc_pca, scoresLoadings = c(T,T))

  # this plot displays the cumulative variance explained by PC1 and PC2
  plotPcs(tc_pca, type = c("loadings"))

  # scores(tc_pca)

  # to get the variance (R2) explained by each axis
  summary(tc_pca)
  pc1_v <- round(pull(as_tibble(summary(tc_pca)), PC1)[1]*100, digits = 1)
  pc2_v <- round(pull(as_tibble(summary(tc_pca)), PC2)[1]*100, digits = 1)
  pc3_v <- round(pull(as_tibble(summary(tc_pca)), PC3)[1]*100, digits = 1)

  pca_loadings <- as_tibble(loadings(tc_pca)) %>%
    mutate(Sample_ID = row.names(loadings(tc_pca))) %>%
    left_join(keyfile, by = "Sample_ID")

  ### set colour scheme if present in the sample.key
  # cat_colour_names <- pull(pca_loadings, colour)
  # names(cat_colour_names) <- pull(pca_loadings, Descriptive_ID)

  ### plots PC1 vs PC2
  g <- ggplot(pca_loadings, aes(x = PC1, y = PC2, colour = as.character(!!as.name(dot_colour)))) +
    geom_point() +
    scale_colour_ptol() +
    theme_classic() +
    labs(title = paste0("Timepoint ", Timepoint, " hr"),
         x = paste0("PC1 (", pc1_v, "%)"),
         y = paste0("PC2 (", pc2_v, "%)"))
  print(g)

  # add text labels
  g_labs <- g + geom_text_repel(aes(label=Sample_ID), show.legend = F)
  print(g_labs)

  g <- ggplot(pca_loadings, aes(x = PC1, y = PC2, colour = as.character(!!as.name(dot_colour)))) +
    geom_point() +
    theme_classic() +
    scale_colour_ptol() +
    labs(title = paste0("Timepoint ", Timepoint, " hr - coord prop"),
         x = paste0("PC1 (", pc1_v, "%)"),
         y = paste0("PC2 (", pc2_v, "%)")) +
    coord_equal()
  print(g)

  ### plots PC1 vs PC3
  g <- ggplot(pca_loadings, aes(x = PC1, y = PC3, colour = as.character(!!as.name(dot_colour)))) +
    geom_point() +
    scale_colour_ptol() +
    theme_classic() +
    labs(title = paste0("Timepoint ", Timepoint, " hr"),
         x = paste0("PC1 (", pc1_v, "%)"),
         y = paste0("PC3 (", pc3_v, "%)"))
  print(g)

  # add text labels
  g_labs <- g + geom_text_repel(aes(label=Sample_ID), show.legend = F)
  print(g_labs)

  g <- ggplot(pca_loadings, aes(x = PC1, y = PC3, colour = as.character(!!as.name(dot_colour)))) +
    geom_point() +
    theme_classic() +
    scale_colour_ptol() +
    labs(title = paste0("Timepoint ", Timepoint, " hr - coord prop"),
         x = paste0("PC1 (", pc1_v, "%)"),
         y = paste0("PC3 (", pc3_v, "%)")) +
    coord_equal()
  print(g)

  ### plots PC2 vs PC3
  g <- ggplot(pca_loadings, aes(x = PC2, y = PC3, colour = as.character(!!as.name(dot_colour)))) +
    geom_point() +
    scale_colour_ptol() +
    theme_classic() +
    labs(title = paste0("Timepoint ", Timepoint, " hr"),
         x = paste0("PC2 (", pc2_v, "%)"),
         y = paste0("PC3 (", pc3_v, "%)"))
  print(g)

  # add text labels
  g_labs <- g + geom_text_repel(aes(label=Sample_ID), show.legend = F)
  print(g_labs)

  g <- ggplot(pca_loadings, aes(x = PC2, y = PC3, colour = as.character(!!as.name(dot_colour)))) +
    geom_point() +
    theme_classic() +
    scale_colour_ptol() +
    labs(title = paste0("Timepoint ", Timepoint, " hr - coord prop"),
         x = paste0("PC2 (", pc2_v, "%)"),
         y = paste0("PC3 (", pc3_v, "%)")) +
    coord_equal()
  print(g)
}

make_PCA_plots_large_labels <- function (dot_colour,
                                         cpm_tbls = cpm_tbl,
                                         keyfile) {

  tc_pca <- pca(
    cpm_tbls[,-1],
    scale = "uv",
    center = T,
    nPcs = 3,
    method = "nipals"
  )

  slplot(tc_pca, scoresLoadings = c(T,T))

  # this plot displays the cumulative variance explained by PC1 and PC2
  plotPcs(tc_pca, type = c("loadings"))

  # to get the variance (R2) explained by each axis
  summary(tc_pca)
  pc1_v <- round(
    pull(
      as_tibble(summary(tc_pca)),
      PC1
    )[1] * 100,
    digits = 1
  )
  pc2_v <- round(
    pull(
      as_tibble(summary(tc_pca)),
      PC2
    )[1] * 100,
    digits = 1
  )
  pc3_v <- round(
    pull(
      as_tibble(summary(tc_pca)),
      PC3
    )[1] * 100,
    digits = 1
  )

  pca_loadings <- as_tibble(loadings(tc_pca)) %>%
    mutate(
      Sample_ID = row.names(loadings(tc_pca))
    ) %>%
    left_join(keyfile, by = "Sample_ID")

  ### plots PC1 vs PC2
  g <- ggplot(
    pca_loadings,
    aes(
      x = PC1,
      y = PC2,
      colour = as.character(!!as.name(dot_colour))
    )
  ) +
    geom_point() +
    scale_colour_ptol() +
    theme_classic() +
    labs(
      x = paste0("PC1 (", pc1_v, "%)"),
      y = paste0("PC2 (", pc2_v, "%)")
    )
  print(g)

  # add text labels
  g_labs <- g +
    geom_text_repel(
      aes(label=Sample_ID),
      show.legend = F
    )
  print(g_labs)

  g <- ggplot(
    pca_loadings,
    aes(
      x = PC1,
      y = PC2,
      colour = as.character(!!as.name(dot_colour))
    )
  ) +
    geom_point() +
    theme_classic() +
    scale_colour_ptol() +
    labs(
      title = "Coordinate prop",
      x = paste0("PC1 (", pc1_v, "%)"),
      y = paste0("PC2 (", pc2_v, "%)")
    ) +
    coord_equal()
  print(g)

  ### plots PC1 vs PC3
  g <- ggplot(
    pca_loadings,
    aes(
      x = PC1,
      y = PC3,
      colour = as.character(!!as.name(dot_colour))
    )
  ) +
    geom_point() +
    scale_colour_ptol() +
    theme_classic() +
    labs(
      x = paste0("PC1 (", pc1_v, "%)"),
      y = paste0("PC3 (", pc3_v, "%)")
    )
  print(g)

  # add text labels
  g_labs <- g +
    geom_text_repel(
      aes(label=Sample_ID),
      show.legend = F
    )
  print(g_labs)

  g <- ggplot(
    pca_loadings,
    aes(
      x = PC1,
      y = PC3,
      colour = as.character(!!as.name(dot_colour))
    )
  ) +
    geom_point() +
    theme_classic() +
    scale_colour_ptol() +
    labs(
      title = "Coordinate prop",
      x = paste0("PC1 (", pc1_v, "%)"),
      y = paste0("PC3 (", pc3_v, "%)")
    ) +
    coord_equal()
  print(g)

  ### plots PC2 vs PC3
  g <- ggplot(
    pca_loadings,
    aes(
      x = PC2,
      y = PC3,
      colour = as.character(!!as.name(dot_colour))
    )
  ) +
    geom_point() +
    scale_colour_ptol() +
    theme_classic() +
    labs(
      x = paste0("PC2 (", pc2_v, "%)"),
      y = paste0("PC3 (", pc3_v, "%)")
    )
  print(g)

  # add text labels
  g_labs <- g +
    geom_text_repel(
      aes(label=Sample_ID),
      show.legend = F
    )
  print(g_labs)

  g <- ggplot(
    pca_loadings,
    aes(
      x = PC2,
      y = PC3,
      colour = as.character(!!as.name(dot_colour))
    )
  ) +
    geom_point() +
    theme_classic() +
    scale_colour_ptol() +
    labs(
      title = "Coordinate prop",
      x = paste0("PC2 (", pc2_v, "%)"),
      y = paste0("PC3 (", pc3_v, "%)")
    ) +
    coord_equal()
  print(g)
}
