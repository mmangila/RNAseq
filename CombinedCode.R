lfc_suffixes <- data.frame(
  Level  =  c(0, 1.5, 2),
  Suffix =  c("_detags.csv",
              "_detags_1point5FC.csv",
              "_detags_2FC.csv")
)

find_de_combined <-  function(combos,
                              lfc_suffixes,
                              combined_folder,
                              annotation,
                              funcs,
                              func_focus,
                              paths) {

  de_genes_summary <- t(combos)
  rownames(de_genes_summary) <- paste(
    de_genes_summary[,1],
    "vs",
    de_genes_summary[,2],
    sep = "."
  )

  sapply(1:length(combos[1,]), function (x) {
    test_name <- paste0(
      combos[1,x],
      ".vs.",
      combos[2,x]
    )

    dir.create(paste0(combined_folder,
                      "/DE_tables/",
                      test_name), showWarnings = FALSE)

    edger_deset <- read.csv(paste0(
      paths[3],
      "/edgeR/DE_tables/",
      test_name,
      "/",
      test_name,
      "_alltags.csv"
    ))


    colnames(edger_deset) <-  paste("edger",
                                    colnames(edger_deset),
                                    sep = "_")
    colnames(edger_deset)[which(colnames(edger_deset) == "edger_X")] <- "X"

    deseq_deset <- read.csv(paste0(
      paths[3],
      "/DESEQ2/DE_tables/",
      test_name,
      "/",
      test_name,
      "_alltags.csv"
    ))
    colnames(deseq_deset) <-  paste("deseq",
                                    colnames(deseq_deset),
                                    sep = "_")
    colnames(deseq_deset)[which(colnames(deseq_deset) == "deseq_X")] = "X"

    sapply(1:3, function (y) {


      print(paste0(
        "Finding differentially expressed genes in ",
        test_name,
        " with a fold change greater than ",
        lfc_suffixes[y,1], "."
      ))

      edger_genes <-  try(read.csv(paste0(paths[3],
                                          "/edgeR/DE_tables/",
                                          test_name,
                                          "/",
                                          test_name,
                                          lfc_suffixes[y,2]))$X)

      if ("try-error" %in% class(edger_genes)) {
        edger_genes <- vector()
      }

      deseq_genes <-  try(read.csv(paste0(paths[3],
                                          "/DESEQ2/DE_tables/",
                                          test_name,
                                          "/",
                                          test_name,
                                          lfc_suffixes[y, 2]))$X)

      if ("try-error" %in% class(deseq_genes)) {
        deseq_genes <- vector()
      }

      if (length(edger_genes) == 0 && length(deseq_genes) == 0) {
        print(paste0(
          "Found no differentially expressed genes in ",
          test_name,
          " with a fold change greater than ",
          lfc_suffixes[y,1], "."
        ))
      } else if (length(edger_genes) == 0) {
        comb_genes <- deseq_genes
      } else if (length(deseq_genes) == 0) {
        comb_genes <- edger_genes
      } else {
        comb_genes <- union(edger_genes, deseq_genes)
      }
      edger_table <- edger_deset[which(edger_deset$X %in% comb_genes), ]
      deseq_table <- deseq_deset[which(deseq_deset$X %in% comb_genes), ]

      gene_table  <- merge(data.table(edger_table,
                                      key = names(edger_table)),
                           data.table(deseq_table,
                                      key = names(deseq_table)))
      colnames(gene_table)[which(colnames(gene_table) == "X")] <- func_focus

      if (annotation) {
        gene_funcs  <- funcs[which(
          funcs[, which(colnames(funcs) == func_focus)] %in%
            comb_genes
        ), ]
        final_table <- try(merge(data.table(gene_funcs,
                                            key = names(gene_funcs)),
                                 data.table(gene_table,
                                            key = names(gene_table))))

        if ("try-error" %in% class(final_table)) {
          final_table <- data.frame(X         = character(),
                                    logFC     = numeric(),
                                    AveExpr   = numeric(),
                                    t         = numeric(),
                                    P.Value   = numeric(),
                                    adj.P.Val = numeric(),
                                    B         = numeric())
        }
      } else {
        final_table <- data.table(gene_table, key = names(gene_table))
      }

      if (length(final_table[, 1]) > 0) {
        write_csv(final_table,
                  file = paste0(combined_folder,
                                "/DE_tables/",
                                test_name,
                                "/",
                                test_name,
                                lfc_suffixes[y, 2]))

        print(table(sign(final_table$edger_logFC)))
      }
    })
    gene_table <- try(merge(data.table(edger_deset,
                                       key = names(edger_deset)),
                            data.table(deseq_deset,
                                       key = names(deseq_deset))))
    colnames(gene_table)[which(colnames(gene_table) == "X")] <- func_focus

    if (annotation) {
      gene_funcs  <-  funcs[which(funcs[,
                                        which(
                                          colnames(funcs) == func_focus
                                        )]
      %in% union(deseq_deset$X, edger_deset$X)), ]

      final_table <- try(merge(data.table(gene_funcs,
                                          key = names(gene_funcs)),
                               data.table(gene_table,
                                          key = names(gene_table))))
    } else {
      final_table <- data.table(gene_table, key = names(gene_table))
    }
    write_csv(
      final_table,
      file = paste0(combined_folder,
                    "/DE_tables/",
                    test_name,
                    "/",
                    test_name,
                    "_alltags.csv"))

  })
}


find_combined_de <-  function(keyfile,
                              group,
                              lfc_suffixes,
                              annotation,
                              func_path,
                              func_focus,
                              paths,
                              go,
                              project_folder) {

  combined_folder <- paste0(
    paths[3],
    "/Combined"
  )
  dir.create(
    combined_folder,
    showWarnings = FALSE
  )
  dir.create(
    paste0(combined_folder, "/DE_tables/"),
    showWarnings = FALSE
  )

  combos <- eval(
    parse(
      text = paste0(
        "combn(as.data.frame(keyfile %>% distinct(",
        group,
        "))[,1],2)"
      )
    )
  )

  if (annotation) {
    if (grepl(".tsv", func_path)) {
      funcs <- read.table(func_path, sep = "\t", header = TRUE)
    } else if (grepl(".csv", func_path)) {
      funcs <- read.csv(file = func_path)
    } else {
      errorCondition("File format not recognised")
    }
  }

  sink(file = paste0(combined_folder, "/DE_tables/de_genes_summary.txt"))
  find_de_combined(combos,
                   lfc_suffixes,
                   combined_folder,
                   annotation,
                   funcs,
                   func_focus,
                   paths)
  sink()

  print("Begin GO analysis")

  if (go == TRUE) {
    devtools::source_url(
      "https://github.com/mmangila/RNAseq/raw/main/topGO_functions.R"
    )

    analyse_go(funcs, func_focus, project_folder, combos, paths)
  }
}

add_gene_funcs <- function(combos, funcs) {
  sapply(seqalong(combos[1, ]), function(x) {

  })
}
