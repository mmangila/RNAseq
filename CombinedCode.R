combine_desets <-  function(combos,
                            lfc_suffixes,
                            combined_folder,
                            comp_num,
                            pval,
                            funcs,
                            func_focus,
                            paths) {

  de_genes_summary <- combos[, comp_num]
  test_name <- paste(
    de_genes_summary[1],
    "vs",
    de_genes_summary[2],
    sep = "."
  )

  print(paste("Combining",
              de_genes_summary[1], "vs",
              de_genes_summary[2], "datasets"))

  print("Get edgeR dataset")
  edger_deset <- read.csv(paste0(paths[3],
                                 "/edgeR/DE_tables/",
                                 test_name,
                                 "/",
                                 test_name,
                                 "_alltags.csv"))
  row.names(edger_deset) <- edger_deset$X

  print("Get DESeq2 dataset")
  deseq_deset <- read.csv(paste0(paths[3],
                                 "/DESEQ2/DE_tables/",
                                 test_name,
                                 "/",
                                 test_name,
                                 "_alltags.csv"))
  row.names(deseq_deset) <- deseq_deset$X

  print("Processing datasets")
  union_deset <- merge(deseq_deset, edger_deset, all = TRUE)
  union_deset <- union_deset %>%
    mutate(logFC           = ifelse(is.na(logFC), 0, logFC),
           log2FoldChange  = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
           adj.P.Val       = ifelse(is.na(adj.P.Val), 1, adj.P.Val),
           padj            = ifelse(is.na(padj), 1, padj),
           Union_logFC     = ifelse(abs(logFC) > abs(log2FoldChange),
                                    logFC, log2FoldChange),
           Union_padj      = pmin(adj.P.Val, padj),
           Intersect_logFC = ifelse(abs(logFC) < abs(log2FoldChange),
                                    logFC, log2FoldChange),
           Intersect_padj  = pmax(padj, adj.P.Val))

  union_deset[, func_focus] <- union_deset$X

  if (!is.null(funcs)) {
    annotated_deset <- merge(funcs, union_deset)
  } else {
    annotated_deset <- union_deset
  }

  annotated_deset <- annotated_deset[, colnames(annotated_deset) != "X"]

  sapply(c("Union", "Intersect"), function (combo_mode) {
  
    combo_vars      <- !grepl(combo_mode, names(annotated_deset))
    padj_col        <- sum(combo_vars)
    logfc_col       <- padj_col - 1
    combined_allset <- annotated_deset[, combo_vars]
  
    combined_allset$logFC <- combined_allset[, logfc_col]
    combined_allset$padj  <- combined_allset[, padj_col]
    combined_allset       <- combined_allset[, c(-logfc_col, -padj_col)]
  
    combined_deset  <- combined_allset[combined_allset$padj < pval, ]
    combined_1.5set <- combined_deset[combined_deset$logFC  < 1.5,  ]
    combined_2set   <- combined_deset[combined_deset$logFC  < 2,    ]

    print(paste("Get", combo_mode, "of edgeR and DESeq2"))

    de_table_path <- paste(combined_folder, combo_mode, "DE_tables", test_name, sep = "/")
    dir.create(de_table_path,
               showWarnings = FALSE, recursive = TRUE)
    write.csv(combined_allset,
              paste(de_table_path, paste0(test_name, "_alltags.csv"), sep = "/"))
    sapply(seq_along(lfc_suffixes$Level), function (x) {
      write.csv(combined_deset[combined_deset$logFC < lfc_suffixes$Level[x], ],
                paste(de_table_path, paste0(test_name, lfc_suffixes$Suffix[x]), sep = "/"))
    })

    print(paste(combo_mode, "found"))
  })
}


find_combined_de <-  function(keyfile,
                              group,
                              lfc_suffixes,
                              pval,
                              funcs,
                              func_focus,
                              paths,
                              go,
                              project_folder,
                              de_logfc,
                              go_pval) {

  combined_folder <- paste0(paths[3], "/Combined")
  
  sapply(c("Union", "Intersect"), function (de_combo) {
    dir.create(paste(combined_folder, de_combo, "GO_tests"),
               showWarnings = FALSE, recursive = TRUE)
  })

  combos <- combn(as.data.frame(keyfile %>% distinct(sample_group))[,1], 2)

  sapply(seq_along(combos[1, ]), function (comparison) {
    combine_desets(combos,
                   lfc_suffixes,
                   combined_folder,
                   comparison,
                   pval,
                   funcs,
                   func_focus,
                   paths)
  })
}
