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

  used_loci <- union(deseq_deset$X, edger_deset$X)

  tmp_data  <- sapply(seq_along(used_loci), function (gene_num) {
  
    gene        <- used_loci[gene_num]
    edger_logfc <- edger_deset[gene, "logFC"]
    deseq_logfc <- deseq_deset[gene, "log2FoldChange"]
    edger_pval  <- edger_deset[gene, "adj.P.Val"]
    deseq_pval  <- deseq_deset[gene, "padj"]

    if (is.na(edger_logfc)) edger_logfc <- 0
    if (is.na(deseq_logfc)) deseq_logfc <- 0
    if (is.na(edger_pval))  edger_pval  <- 1
    if (is.na(deseq_pval))  deseq_pval  <- 1

    return(c(max(edger_logfc,deseq_logfc), min(edger_pval, deseq_pval),
             min(edger_logfc,deseq_logfc), max(edger_pval, deseq_pval)))

    print(paste("Processed", paste0(gene_num, "/", length(used_loci)), genes))
  })

  union_deset <- data.frame(X               = used_loci,
                            Union_logFC     = tmp_data[1, ],
                            Union_padj      = tmp_data[2, ],
                            Intersect_logFC = tmp_data[3, ],
                            Intersect_padj  = tmp_data[4, ])

  union_deset[, func_focus] <- union_deset$X

  if (!is.null(funcs)) {
    annotated_deset <- merge(funcs[funcs[, func_focus] %in% used_loci, ], union_deset,
                             by = func_focus)
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
      write.csv(combined_deset[combined_deset$logFC < lfc_Suffixes$Level[x], ],
                paste(de_table_path, paste0(test_name, lfc_Suffixes$Suffix[x]), sep = "/"))
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
