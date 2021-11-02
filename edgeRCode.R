find_de_edger <- function (old.dge, group, keyfile, paths, analysis) {
  dge <- calcNormFactors(old.dge, method="TMM")
  design <- eval(parse(text = paste0('model.matrix(~0 + ',
                                     group,
                                     ', data = keyfile)')))
  colnames(design) <- eval(parse(text = paste0("levels(as.factor(keyfile$",group,"))")))
  v <- voom(dge,design,plot=TRUE)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)

  print("Begin MDS")
  de_edger_mds(dge, group, keyfile, paths, v)
  print("Begin PCA")
  de_edger_pca(v, group, keyfile, paths)
  print("Begin DE gene table generation")
  gene.names <- as.character(rownames(dge$counts))
  de_edger_tables(keyfile, group, fit, paths, design, gene.names, analysis)
}

de_edger_mds <- function (dge, group, keyfile, paths, v) {
  SampleColours <- ptol_pal()(length(dge$samples$lib.size)/2)
  if (length(dge$samples$lib.size) %% 2 == 1) {
    dddSampleColours <- c(rep(SampleColours, each=2), ptol_pal()(1))
  } else {
    dddSampleColours <-c(rep(SampleColours, each=2))
  }

  mds2 <- plotMDS.invisible(v, ndim=3)

  # Make MDS plots
  par(mfrow=c(2,2), ps=10)
  sapply(1:length(combn(1:3,2)[1,]), function (x) {
    with(mds2,
         eval(parse(text = paste0("plot(V",
                                  combn(1:3,2)[1, x],
                                  ", V",
                                  combn(1:3,2)[2, x],
                                  ", pch=16,",
                                  " col=SampleColours, ",
                                  " cex=1,",
                                  "main='MDS plot',",
                                  "xlab='Leading logFC dim",
                                  combn(1:3,2)[1, x],
                                  "',",
                                  "ylab='Leading logFC dim ",
                                  combn(1:3,2)[2, x],
                                  "')"))))

  })

  mds2$group <- keyfile[, which(colnames(keyfile) == group)]

  with(mds2, {
    s3d <- scatterplot3d(V1, V3, V2,
                         color = dddSampleColours,
                         pch=19,
                         cex.symbols = 1.5,
                         type="h",
                         main="3D MDS plot",
                         xlab="Leading logFC dim 1",
                         ylab="Leading logFC dim 3",
                         zlab="Leading logFC dim 2")
    s3d.coords <- s3d$xyz.convert(V1, V3, V2)
    text(s3d.coords$x,
         s3d.coords$y,
         labels=row.names(mds2),
         cex=.5, pos=4)
  })

  par(mfrow=c(1,1))

  pdf(paste0(paths[3],"/MDS/voom_mds_3d.pdf"), width = 6.5, height = 6.5)
  with(mds2, {
    s3d <- scatterplot3d(V1, V3, V2,
                         pch=19,
                         cex.axis = 0.5,
                         cex.symbols = 1,
                         cex.lab = 0.5,
                         type="h",
                         xlab="Leading logFC dim 1",
                         ylab="Leading logFC dim 3",
                         zlab="Leading logFC dim 2")
    s3d.coords <- s3d$xyz.convert(V1, V3, V2)
    text(s3d.coords$x,
         s3d.coords$y,
         labels=row.names(mds2),
         cex=.25, pos=4)
  })
  dev.off()
}

de_edger_pca <- function (v, group, keyfile, paths) {
  voom_matrix <- v$E

  CPM_tbl <- as_tibble(voom_matrix, rownames = "Gene")

  # replace NaN with NA? Nope there arent any...
  which(is.na(CPM_tbl))
  # CPM_tbl[is.na(CPM_tbl)] <- NA

  MDS_table <- plotMDS(CPM_tbl[,-1], plot = F, ndim = 5, top = 500)
  # MDS_table <- plotMDS(v, plot = F, ndim = 5)

  cmdscale_out <- as_tibble(MDS_table$cmdscale.out)

  cmdscale_out <- cmdscale_out %>%
    mutate(Sample_ID = colnames(v)) %>%
    left_join(keyfile, by = "Sample_ID")

  pdf(paste0(paths[3],"/MDS/voom_e-counts_PCA_Treatment.pdf"), width = 5, height = 3.5)
  make_PCA_plots(Timepoint = "ALL", dot_colour = group, keyfile, CPM_tbl)
  dev.off()

  pdf(paste0(paths[3],"/MDS/voom_e-counts_PCA_Treatment_labels.pdf"), width = 10, height = 7)
  make_PCA_plots_large_labels(
    dot_colour = group,
    CPM_tbl,
    keyfile
    )
  dev.off()
}

de_edger_tables <- function (keyfile, group, fit, paths, design, gene.names, analysis) {
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
  inside2 <- sapply(1:length(combos[1,]), function (x) {
    inside <- c(
      inside,
      paste0(
        combos[1,x],
        ".vs.",
        combos[2,x],
        " = ",
        combos[1,x],
        " - ",
        combos[2,x]
        )
      )
    return(inside)
    })
  contrast.matrix <- eval(
    parse(
      text = paste0(
        "makeContrasts(",
        paste(inside2, collapse = ","),
        ", levels = design)"
        )
      )
    )

  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  results <- decideTests(fit2)
  results2FC <- decideTests(fit2, lfc = log2(2))
  results1.5FC <- decideTests(fit2, lfc = log2(1.5))

  print("summary(results)")
  print(summary(results))
  print("summary(results1.5FC)")
  print(summary(results1.5FC))
  print("summary(results2FC)")
  print(summary(results2FC))

  DETableAll <- topTable(fit=fit2, number=Inf)
  DETable2FC <- topTable(fit=fit2, number=Inf, lfc = 1)

  coefs <- colnames(fit2$coefficients)
  out.base <- paste0(paths[3],"/edgeR/DE_tables/")
  dir.create(out.base, recursive=T, showWarnings = F)

  for (comparison in coefs){
    test.name <- comparison
    test.base.dir <- paste0(out.base, test.name, "/")
    dir.create(test.base.dir, showWarnings=F)
    print(test.name)

    pdf(paste0(test.base.dir, test.name, "_volcano.pdf"))
    plotTable <- topTable(fit=fit2, number=Inf, coef = comparison, sort.by = "logFC")
    with(plotTable, plot(logFC, -(log(adj.P.Val)), pch = 16, cex = 0.3))
    title(main=paste0("Foldchange vs FDR (adjusted p-value)", comparison), cex.main=0.5)
    plotTable <- topTable(fit=fit2, number=Inf, coef=comparison, sort.by = "logFC", p.value = 0.05, lfc = log2(2))
    # write if statement to account for when there are no DE genes
    if(length(plotTable) > 0){
      with(plotTable, points(logFC, -(log(adj.P.Val)), pch = 16,
                             col = "red", cex = 0.5))
      dev.off()
    } else {
      dev.off()
    }

    results <- decideTests(fit2)
    results2FC <- decideTests(fit2, lfc = log2(2))
    results1.5FC <- decideTests(fit2, lfc = log2(1.5))

    pdf(paste0(test.base.dir, test.name, "_smear.pdf"))
    MAplotGeneSetLimma(
      MArrayLMobject = fit2,
      resultsMatrix=results,
      geneSetList = comparison,
      geneSetListName = comparison,
      inputList = comparison,
      inputListName = paste0("Up-regulated genes", comparison)
    )
    dev.off()

    pdf(paste0(test.base.dir, test.name, "_smear1point5FC.pdf"))
    MAplotGeneSetLimma(
      MArrayLMobject = fit2,
      resultsMatrix=results1.5FC,
      geneSetList = comparison,
      geneSetListName = comparison,
      inputList = comparison,
      inputListName = paste0("Up-regulated genes (1.5FC)", comparison)
    )
    dev.off()

    pdf(paste0(test.base.dir, test.name, "_smear2FC.pdf"))
    MAplotGeneSetLimma(
      MArrayLMobject = fit2,
      resultsMatrix=results2FC,
      geneSetList = comparison,
      geneSetListName = comparison,
      inputList = comparison,
      inputListName = paste0("Up-regulated genes (2FC)", comparison)
    )
    dev.off()

    tt <- topTable(fit=fit2, number=Inf, coef=comparison)
    write.csv(tt, paste0(test.base.dir, test.name, "_alltags.csv"))
    tt <- topTable(fit=fit2, number=Inf, coef=comparison, p.value = 0.05)
    write.csv(tt, paste0(test.base.dir, test.name, "_detags.csv"))
    tt <- topTable(fit=fit2, number=Inf, coef=comparison, p.value = 0.05, lfc = log2(1.5))
    write.csv(tt, paste0(test.base.dir, test.name, "_detags_1point5FC.csv"))
    tt <- topTable(fit=fit2, number=Inf, coef=comparison, p.value = 0.05, lfc = log2(2))
    write.csv(tt, paste0(test.base.dir, test.name, "_detags_2FC.csv"))
  }

  tests <- mclapply(coefs, function (x) topTable(fit=fit2, number=Inf, coef=x, sort.by = "none"))
  names(tests) <- coefs
  fc.matrix <- as.data.frame(sapply(tests, function (t) t$logFC, simplify = "array"))
  fdr.matrix <- sapply(tests, function (t) t$adj.P.Val, simplify = "array")
  #be careful that gene.names is in the correct order...
  rownames(fc.matrix) <- gene.names
  rownames(fdr.matrix) <- gene.names
  write.csv(fc.matrix, file=paste0(out.base, analysis, "_fc.csv"))
  write.csv(fdr.matrix, file=paste0(out.base, analysis, "_fdr.csv"))

  tests.sig.2FC <- mclapply(coefs, function (x) topTable(fit=fit2, p.value = 0.05, lfc = log2(2), number=Inf, coef=x, sort.by = "none"))
  names(tests.sig.2FC) <- coefs

  tests.sig.1.5FC <- mclapply(coefs, function (x) topTable(fit=fit2, p.value = 0.05, lfc = log2(1.5), number=Inf, coef=x, sort.by = "none"))
  names(tests.sig.1.5FC) <- coefs

  cpm.matrix <- cpm(v)
  write.csv(fdr.matrix, file=paste0(out.base, analysis, "_fdr.csv"))
}

plotMDS.invisible <- function(...){
  ff <- tempfile()
  png(filename=ff)
  mds <- plotMDS(...)
  dev.off()
  unlink(ff)
  mds2 <- as.data.frame(mds$cmdscale.out)
  mds2
}


######### PCA 1 FUNCTION - JUST DOTS - SMALL PDF

make_PCA_plots <- function(Timepoint = "ALL", dot_colour, keyfile, CPM_tbl){
  # Timepoint = "ALL"
  # dot_colour = "Time"

  # If else to subset data by timepoint
  if(Timepoint == "ALL"){

    TC_PCA <- pca(CPM_tbl[,-1], scale = "uv", center = T, nPcs = 3, method = "nipals")

  }else{

    CPM_table_group <- CPM_tbl %>%
      gather(key = Sample_ID, value = log2CPM, -Gene) %>%
      left_join(keyfile, by = "Sample_ID")

    CPM_table_group_spread <- CPM_table_group %>%
      select(Gene, Sample_ID, log2CPM) %>%
      spread(key = Sample_ID, value = log2CPM)

    PCA_input_table <- CPM_table_group_spread[,-1]

    TC_PCA <- pca(PCA_input_table, scale = "uv", center = T, nPcs = 3, method = "nipals")
  }

  # to get the variance (R2) explained by each axis
  summary(TC_PCA)
  PC1_v <- round(pull(as.tibble(summary(TC_PCA)), PC1)[1]*100, digits = 1)
  PC2_v <- round(pull(as.tibble(summary(TC_PCA)), PC2)[1]*100, digits = 1)
  PC3_v <- round(pull(as.tibble(summary(TC_PCA)), PC3)[1]*100, digits = 1)

  PCA_loadings <- as.tibble(loadings(TC_PCA)) %>%
    mutate(Sample_ID = row.names(loadings(TC_PCA))) %>%
    left_join(keyfile, by = "Sample_ID")

  ### plots PC1 vs PC2
  g <- ggplot(
    PCA_loadings,
    aes(x = PC1, y = PC2, colour = as.character(!!as.name(dot_colour)))
    ) +
    geom_point() +
    scale_colour_ptol() +
    theme_classic() +
    labs(title = paste0("Timepoint ", Timepoint, " hr"),
         x = paste0("PC1 (", PC1_v, "%)"),
         y = paste0("PC2 (", PC2_v, "%)"))
  print(g)

  # add text labels
  g_labs <- g + geom_text_repel(aes(label=Sample_ID), show.legend = F)
  # print(g_labs)

  g <- ggplot(
    PCA_loadings,
    aes(x = PC1, y = PC2, colour = as.character(!!as.name(dot_colour)))
    ) +
    geom_point() +
    theme_classic() +
    scale_colour_ptol() +
    labs(title = paste0("Timepoint ", Timepoint, " hr - coord prop"),
         x = paste0("PC1 (", PC1_v, "%)"),
         y = paste0("PC2 (", PC2_v, "%)")) +
    coord_equal()
  print(g)

  ### plots PC1 vs PC3
  g <- ggplot(
    PCA_loadings,
    aes(x = PC1, y = PC3, colour = as.character(!!as.name(dot_colour)))
    ) +
    geom_point() +
    scale_colour_ptol() +
    theme_classic() +
    labs(title = paste0("Timepoint ", Timepoint, " hr"),
         x = paste0("PC1 (", PC1_v, "%)"),
         y = paste0("PC3 (", PC3_v, "%)"))
  print(g)

  # add text labels
  g_labs <- g + geom_text_repel(aes(label=Sample_ID), show.legend = F)
  # print(g_labs)

  g <- ggplot(
    PCA_loadings,
    aes(x = PC1, y = PC3, colour = as.character(!!as.name(dot_colour)))
    ) +
    geom_point() +
    theme_classic() +
    scale_colour_ptol() +
    labs(title = paste0("Timepoint ", Timepoint, " hr - coord prop"),
         x = paste0("PC1 (", PC1_v, "%)"),
         y = paste0("PC3 (", PC3_v, "%)")) +
    coord_equal()
  print(g)

  ### plots PC2 vs PC3
  g <- ggplot(PCA_loadings, aes(x = PC2, y = PC3, colour = as.character(!!as.name(dot_colour)))) +
    geom_point() +
    scale_colour_ptol() +
    theme_classic() +
    labs(title = paste0("Timepoint ", Timepoint, " hr"),
         x = paste0("PC2 (", PC2_v, "%)"),
         y = paste0("PC3 (", PC3_v, "%)"))
  print(g)

  # add text labels
  g_labs <- g + geom_text_repel(aes(label=Sample_ID), show.legend = F)
  # print(g_labs)

  g <- ggplot(PCA_loadings, aes(x = PC2, y = PC3, colour = as.character(!!as.name(dot_colour)))) +
    geom_point() +
    theme_classic() +
    scale_colour_ptol() +
    labs(title = paste0("Timepoint ", Timepoint, " hr - coord prop"),
         x = paste0("PC2 (", PC2_v, "%)"),
         y = paste0("PC3 (", PC3_v, "%)")) +
    coord_equal()
  print(g)
}


######### ######### ######### ######### #########
######### PCA 2 FUNCTION - labels - LARGE PDF

make_PCA_plots_large_labels <- function(Timepoint = "ALL", dot_colour, keyfile, CPM_tbl){

  if(Timepoint == "ALL"){

    TC_PCA <- pca(CPM_tbl[,-1], scale = "uv", center = T, nPcs = 3, method = "nipals")

  }else{

    CPM_table_group <- CPM_tbl %>%
      gather(key = Sample_ID, value = log2CPM, -Gene) %>%
      left_join(keyfile, by = "Sample_ID") %>%
      filter(Time == Timepoint)

    CPM_table_group_spread <- CPM_table_group %>%
      select(Gene, Sample_ID, log2CPM) %>%
      spread(key = Sample_ID, value = log2CPM)

    PCA_input_table <- CPM_table_group_spread[,-1]

    TC_PCA <- pca(PCA_input_table, scale = "uv", center = T, nPcs = 3, method = "nipals")
  }

  slplot(TC_PCA, scoresLoadings = c(T,T))

  # this plot displays the cumulative variance explained by PC1 and PC2
  plotPcs(TC_PCA, type = c("loadings"))

  # scores(TC_PCA)

  # to get the variance (R2) explained by each axis
  summary(TC_PCA)
  PC1_v <- round(pull(as.tibble(summary(TC_PCA)), PC1)[1]*100, digits = 1)
  PC2_v <- round(pull(as.tibble(summary(TC_PCA)), PC2)[1]*100, digits = 1)
  PC3_v <- round(pull(as.tibble(summary(TC_PCA)), PC3)[1]*100, digits = 1)

  PCA_loadings <- as.tibble(loadings(TC_PCA)) %>%
    mutate(Sample_ID = row.names(loadings(TC_PCA))) %>%
    left_join(keyfile, by = "Sample_ID")

  ### set colour scheme if present in the sample.key
  # cat_colour_names <- pull(PCA_loadings, colour)
  # names(cat_colour_names) <- pull(PCA_loadings, Descriptive_ID)

  ### plots PC1 vs PC2
  g <- ggplot(PCA_loadings, aes(x = PC1, y = PC2, colour = as.character(!!as.name(dot_colour)))) +
    geom_point() +
    scale_colour_ptol() +
    theme_classic() +
    labs(title = paste0("Timepoint ", Timepoint, " hr"),
         x = paste0("PC1 (", PC1_v, "%)"),
         y = paste0("PC2 (", PC2_v, "%)"))
  print(g)

  # add text labels
  g_labs <- g + geom_text_repel(aes(label=Sample_ID), show.legend = F)
  print(g_labs)

  g <- ggplot(PCA_loadings, aes(x = PC1, y = PC2, colour = as.character(!!as.name(dot_colour)))) +
    geom_point() +
    theme_classic() +
    scale_colour_ptol() +
    labs(title = paste0("Timepoint ", Timepoint, " hr - coord prop"),
         x = paste0("PC1 (", PC1_v, "%)"),
         y = paste0("PC2 (", PC2_v, "%)")) +
    coord_equal()
  print(g)

  ### plots PC1 vs PC3
  g <- ggplot(PCA_loadings, aes(x = PC1, y = PC3, colour = as.character(!!as.name(dot_colour)))) +
    geom_point() +
    scale_colour_ptol() +
    theme_classic() +
    labs(title = paste0("Timepoint ", Timepoint, " hr"),
         x = paste0("PC1 (", PC1_v, "%)"),
         y = paste0("PC3 (", PC3_v, "%)"))
  print(g)

  # add text labels
  g_labs <- g + geom_text_repel(aes(label=Sample_ID), show.legend = F)
  print(g_labs)

  g <- ggplot(PCA_loadings, aes(x = PC1, y = PC3, colour = as.character(!!as.name(dot_colour)))) +
    geom_point() +
    theme_classic() +
    scale_colour_ptol() +
    labs(title = paste0("Timepoint ", Timepoint, " hr - coord prop"),
         x = paste0("PC1 (", PC1_v, "%)"),
         y = paste0("PC3 (", PC3_v, "%)")) +
    coord_equal()
  print(g)

  ### plots PC2 vs PC3
  g <- ggplot(PCA_loadings, aes(x = PC2, y = PC3, colour = as.character(!!as.name(dot_colour)))) +
    geom_point() +
    scale_colour_ptol() +
    theme_classic() +
    labs(title = paste0("Timepoint ", Timepoint, " hr"),
         x = paste0("PC2 (", PC2_v, "%)"),
         y = paste0("PC3 (", PC3_v, "%)"))
  print(g)

  # add text labels
  g_labs <- g + geom_text_repel(aes(label=Sample_ID), show.legend = F)
  print(g_labs)

  g <- ggplot(PCA_loadings, aes(x = PC2, y = PC3, colour = as.character(!!as.name(dot_colour)))) +
    geom_point() +
    theme_classic() +
    scale_colour_ptol() +
    labs(title = paste0("Timepoint ", Timepoint, " hr - coord prop"),
         x = paste0("PC2 (", PC2_v, "%)"),
         y = paste0("PC3 (", PC3_v, "%)")) +
    coord_equal()
  print(g)
}

make_PCA_plots_large_labels <- function (dot_colour,
                                         CPM_tbls = CPM_tbl,
                                         keyfile) {

  TC_PCA <- pca(
    CPM_tbls[,-1],
    scale = "uv",
    center = T,
    nPcs = 3,
    method = "nipals"
  )

  slplot(TC_PCA, scoresLoadings = c(T,T))

  # this plot displays the cumulative variance explained by PC1 and PC2
  plotPcs(TC_PCA, type = c("loadings"))

  # to get the variance (R2) explained by each axis
  summary(TC_PCA)
  PC1_v <- round(
    pull(
      as.tibble(summary(TC_PCA)),
      PC1
    )[1] * 100,
    digits = 1
  )
  PC2_v <- round(
    pull(
      as.tibble(summary(TC_PCA)),
      PC2
    )[1] * 100,
    digits = 1
  )
  PC3_v <- round(
    pull(
      as.tibble(summary(TC_PCA)),
      PC3
    )[1] * 100,
    digits = 1
  )

  PCA_loadings <- as.tibble(loadings(TC_PCA)) %>%
    mutate(
      Sample_ID = row.names(loadings(TC_PCA))
    ) %>%
    left_join(keyfile, by = "Sample_ID")

  ### plots PC1 vs PC2
  g <- ggplot(
    PCA_loadings,
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
      x = paste0("PC1 (", PC1_v, "%)"),
      y = paste0("PC2 (", PC2_v, "%)")
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
    PCA_loadings,
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
      x = paste0("PC1 (", PC1_v, "%)"),
      y = paste0("PC2 (", PC2_v, "%)")
    ) +
    coord_equal()
  print(g)

  ### plots PC1 vs PC3
  g <- ggplot(
    PCA_loadings,
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
      x = paste0("PC1 (", PC1_v, "%)"),
      y = paste0("PC3 (", PC3_v, "%)")
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
    PCA_loadings,
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
      x = paste0("PC1 (", PC1_v, "%)"),
      y = paste0("PC3 (", PC3_v, "%)")
    ) +
    coord_equal()
  print(g)

  ### plots PC2 vs PC3
  g <- ggplot(
    PCA_loadings,
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
      x = paste0("PC2 (", PC2_v, "%)"),
      y = paste0("PC3 (", PC3_v, "%)")
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
    PCA_loadings,
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
      x = paste0("PC2 (", PC2_v, "%)"),
      y = paste0("PC3 (", PC3_v, "%)")
    ) +
    coord_equal()
  print(g)
}
