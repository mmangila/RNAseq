lfc.suffixes <- data.frame(
  Level = c(0,1.5,2),
  Suffix = c("_detags.csv",
             "_detags_1point5FC.csv",
             "_detags_2FC.csv"
  )
)

find_all_de <- function () {
  combos <- eval(
    parse(
      text = paste0(
        "combn(as.data.frame(keyfile %>% distinct(",
        group,
        "))[,1],2)"
      )
    )
  )

  de_genes_summary <- t(combos)
  rownames(de_genes_summary) <- paste(
    de_genes_summary[,1],
    "vs",
    de_genes_summary[,1],
    sep = "."
  )
}

find_de_intersects <- function (combos) {

  lfc.suffixes <- data.frame(
    Level = c(0,1.5,2),
    Suffix = c("_detags.csv",
               "_detags_1point5FC.csv",
               "_detags_2FC.csv"
    )
  )

  de_summary <- data.frame()


  sapply(1:3, function (lfc) {

    lfc_summary <- eval(
      parse(
        text = paste0(
          "data.frame(",
          paste0("down_", lfc.suffixes[lfc, 1]),
          " = integer(), ",
          paste0("up_", lfc.suffixes[lfc, 1]),
          " = integer())"
        )
      )
    )



    sapply(1:length(combos[1,]), function (test) {
      test.name <- paste0(
        combos[1,test],
        ".vs.",
        combos[2,test]
      )

      summary_test <- find_de_intersect(test.name, lfc.suffixes, 0, hisat_paths)

      eval(
        parse(
          text = paste0(
            "rbind(lfc_summary, ",
            test.name,
            " = summary_test)"
          )
        )
      )
    })

  })
}

find_de_intersect <- function (test.name, lfc.suffixes, y, paths) {
  edger_genes <- read.csv(paste0(paths[3],
                                 "/edgeR/DE_tables/",
                                 test.name,
                                 "/",
                                 test.name,
                                 lfc.suffixes[y,2]))$X

  deseq_genes <- read.csv(paste0(paths[3],
                                 "/DESEQ2/DE_tables/",
                                 test.name,
                                 "/",
                                 test.name,
                                 lfc.suffixes[y,2]))$X

  comb_genes <- union(edger_genes,deseq_genes)

  edger_table <- edger.deset[which(edger.deset$X %in% comb_genes),]
  deseq_table <- deseq.deset[which(deseq.deset$X %in% comb_genes),]

  final_table <- merge(data.table(edger_table,
                                  key = names(edger_table)),
                       data.table(deseq_table,
                                  key = names(deseq_table)))

  # write_csv(
  #   final_table,
  #   file = paste0(combined.folder,
  #                 "/DE_tables/",
  #                 test.name,
  #                 "/",
  #                 test.name,
  #                 lfc.suffixes[y,2]))

  return(as.vector(table(sign(final_table$edger_logFC))))
}

find_de_combined <- function (combos, lfc.suffixes, combined.folder, funcs, func_focus, paths) {

  de_genes_summary <- t(combos)
  rownames(de_genes_summary) <- paste(
    de_genes_summary$V1,
    "vs",
    de_genes_summary$V2,
    sep = "."
  )

  sapply(1:length(combos[1,]), function (x) {
    test.name <- paste0(
      combos[1,x],
      ".vs.",
      combos[2,x]
    )

    dir.create(paste0(combined.folder,"/DE_tables/",test.name), showWarnings = FALSE)

    edger.deset <- read.csv(paste0(
      paths[3],
      "/edgeR/DE_tables/",
      test.name,
      "/",
      test.name,
      "_alltags.csv"
    ))


    colnames(edger.deset) <- paste("edger",
                                   colnames(edger.deset),
                                   sep = "_")
    colnames(edger.deset)[which(colnames(edger.deset) == "edger_X")] = "X"

    deseq.deset <- read.csv(paste0(
      paths[3],
      "/DESEQ2/DE_tables/",
      test.name,
      "/",
      test.name,
      "_alltags.csv"
    ))
    colnames(deseq.deset) <- paste("deseq",
                                   colnames(deseq.deset),
                                   sep = "_")
    colnames(deseq.deset)[which(colnames(deseq.deset) == "deseq_X")] = "X"

    sapply(1:3, function (y) {


        print(paste0(
          "Finding differentially expressed genes in ",
          test.name,
          " with a fold change greater than ",
          lfc.suffixes[y,1], "."
        ))

        edger_genes <- read.csv(paste0(paths[3],
                                       "/edgeR/DE_tables/",
                                       test.name,
                                       "/",
                                       test.name,
                                       lfc.suffixes[y,2]))$X

        deseq_genes <- read.csv(paste0(paths[3],
                                       "/DESEQ2/DE_tables/",
                                       test.name,
                                       "/",
                                       test.name,
                                       lfc.suffixes[y,2]))$X

        comb_genes <- union(edger_genes,deseq_genes)

        edger_table <- edger.deset[which(edger.deset$X %in% comb_genes),]
        deseq_table <- deseq.deset[which(deseq.deset$X %in% comb_genes),]

        gene_table  <- merge(data.table(edger_table,
                                        key = names(edger_table)),
                             data.table(deseq_table,
                                        key = names(deseq_table)))
        colnames(gene_table)[which(colnames(gene_table) == "X")] <- func_focus

        gene_funcs  <- funcs[which(funcs$X %in% comb_genes), ]
        final_table <- merge(data.table(gene_funcs, key = names(gene_funcs)),
                             data.table(gene_table, key = names(gene_table)))

        write_csv(
          final_table,
          file = paste0(combined.folder,
                        "/DE_tables/",
                        test.name,
                        "/",
                        test.name,
                        lfc.suffixes[y,2]))



        print(table(sign(final_table$edger_logFC)))
      })

    })
  }


find_combined_de <- function(keyfile, group, lfc.suffixes, func_path, func_focus, paths) {

  combined.folder <- paste0(
    paths[3],
    "/Combined"
  )
  dir.create(
    combined.folder,
    showWarnings = FALSE
  )
  dir.create(
    paste0(combined.folder,"/DE_tables/"),
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

  funcs <- read.csv(file = func_path)

  sink(file = paste0(combined.folder, "/DE_tables/de_genes_summary.txt"))
  find_de_combined(combos, lfc.suffixes, combined.folder, funcs, func_focus, paths)
  sink()
}

add_gene_funcs <- function (combos, funcs) {
  sapply(1:length(combos[1, ]), function (x) {

  })
}
