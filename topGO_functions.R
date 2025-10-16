analyse_topgo <- function(ontology_to_test,
                          de_direction,
                          test_name,
                          de_table_file,
                          out_dir,
                          func_focus,
                          geneid_to_go) {

  ### args
  # ontology_to_test = "BP"
  # test_name = "Gall_base.vs.Ungalled_leaf"
  # test_number = 3
  # de_table_file = paste0("FibrosaAlignment/FibrosaAlignment_batch1/DE_analysis_2022-05-23-1743/Combined/DE_tables/Surrounding_tissue.vs.Ungalled_leaf/Surrounding_tissue.vs.Ungalled_leaf_alltags.csv")
  # # need a column called "model", "adj.P.Val" and "logFC"
  # de_direction = "down"
  # out_dir = paste0("FibrosaAlignment/FibrosaAlignment_batch1/DE_analysis_2022-05-23-1743/Combined/GO_tests/")
  # #####

  dir.create(out_dir, showWarnings = FALSE)
  de_table <- read_csv(de_table_file)
  gene_universe <- de_table %>% pull(func_focus)

  if (de_direction == "up") {
    de_locus <- as.data.frame(
      de_table[which(de_table$edger_logFC > 1), 1]
    )[, 1]
  } else if (de_direction == "down") {
    de_locus <- as.data.frame(
      de_table[which(de_table$edger_logFC < 1), 1]
    )[, 1]
  } else {
    print("de direction confusion")
  }

  gene_list <- factor(as.integer(gene_universe %in% de_locus))
  names(gene_list) <- gene_universe
  str(gene_list)

  ### build go data

  go_data <- new(
    "topGOdata",
    description = test_name,
    ontology    = ontology_to_test,
    allGenes    = gene_list,
    annot       = annFUN.gene2GO,
    nodeSize    = 10,
    gene2GO     = geneid_to_go
  )

  result_fisher <- runTest(go_data,
                           algorithm = "classic",
                           statistic = "fisher")

  top_res       <- GenTable(go_data,
                            classicFisher = result_fisher,
                            topNodes = "10",
                            numChar = 1000,
                            orderBy = "classicFisher")

  # to plot the result

  all_go        <- usedGO(object = go_data)

  all_res       <- as.tibble(GenTable(go_data,
                                      Fisher_p_value = result_fisher,
                                      topNodes = length(all_go),
                                      numChar = 1000,
                                      orderBy = "classicFisher"))

  all_res_fdr   <- all_res %>%
    mutate(fdr = p.adjust(Fisher_p_value, method = "BH"),
           Aspect = ontology_to_test)

  sig_res_fdr   <- all_res_fdr %>% filter(fdr < 0.05)

  write_csv(sig_res_fdr,
            paste0(out_dir,
                   "/",
                   test_name,
                   "_",
                   de_direction,
                   "_",
                   ontology_to_test,
                   ".csv"))

  #### Make graphs
  result_ks      <- runTest(go_data, algorithm = "classic", statistic = "ks")
  result_ks_elim <- runTest(go_data, algorithm = "elim",    statistic = "ks")

  node_path <- paste0(out_dir,
                      "/",
                      test_name,
                      "_",
                      de_direction,
                      "_",
                      ontology_to_test,
                      "_nodegraphs.pdf")

  pdf(node_path)
  par(cex = 0.5)
  showSigOfNodes(go_data,
                 score(result_ks_elim),
                 firstSigNodes = 5,
                 useInfo       = "all")
  printGraph(go_data,
             result_fisher,
             firstSigNodes = 5,
             fn.prefix     = "tGO",
             useInfo       = "all",
             pdfSW         = TRUE)
  dev.off()

}



go_plot_comparison <- function(ontology_to_test,
                               de_direction,
                               test_name,
                               big_data,
                               data_folder) {

  #Arguments
  # plot_comparison = "Gall_base.vs.Surrounding_tissue_down_BP"

  plot_comparison <- paste0(test_name,
                            "_",
                            de_direction,
                            "_",
                            ontology_to_test)

  big_data_selection <- big_data %>%
    filter(sample %in% plot_comparison) %>%
    mutate(
      enrichment_ratio = Significant/Expected,
      GO_term = paste0(GO.ID, " ", Term)
    ) %>%
    arrange(enrichment_ratio) %>%
    mutate(GO_term = factor(GO_term, levels = GO_term))

  g <- ggplot(
    big_data_selection,aes(
      y = GO_term,
      x = enrichment_ratio,
      size = Significant,colour = fdr
    )
  ) +
    geom_point() +
    theme_minimal() +
    scale_colour_viridis(
      option = "A",
      direction = -1,
      begin = 0.2,
      end = 0.8
    ) +
    labs(
      size = "Number of genes",
      colour='FDR',
      y = NULL,
      x = "Relative enrichment (Observed/Expected)"
    )

  pdf(paste0(data_folder, "/", plot_comparison, "_GO_plot.pdf"))
  print(g)
  dev.off()
}

analyse_go <- function (funcs, func_focus, project_folder, combos, project_paths) {

  geneDescription_GO <- funcs %>%
    dplyr::select(func_focus, GO) %>%
    dplyr::rename(DB_Object_ID = func_focus, GOid_list = GO)
  geneDescription_GO <- geneDescription_GO %>% filter(!is.na(GOid_list))

  write_tsv(geneDescription_GO,
            paste0(project_folder, "/", project_folder, "_readMappings.tsv"),
            col_names = F)

  geneid_to_go <- readMappings(file = paste0(project_folder,
                                          "/",
                                          project_folder,
                                          "_readMappings.tsv"))

  go_opts <- expand.grid(c("BP", "MF", "CC"), c("up", "down"))

  sapply(1:length(combos[1, ]), function (x) {
    sapply(1:length(go_opts[, 1]), function (y) {
      analyse_topgo(as.character(go_opts[y, 1]),
                    as.character(go_opts[y, 2]),
                    paste0(combos[1,x], ".vs.", combos[2,x]),
                    paste0(project_paths[3],
                           "/Combined/DE_tables/",
                           paste0(combos[1,x], ".vs.", combos[2,x]),
                           "/",
                           paste0(combos[1,x], ".vs.", combos[2,x]),
                           "_detags_1point5FC.csv"),
                    paste0(project_paths[3], "/Combined/GO_tests/"),
                    func_focus,
                    geneid_to_go)
    })
    # map2(
    #   rep(c("BP", "MF", "CC"), 2), rep(c("up", "down"), each = 3),
    #   analyse_topgo,
    #   test_name     = paste0(combos[1,x], ".vs.", combos[2,x]),
    #   de_table_file = paste0(project_paths[3],
    #                          "/Combined/DE_tables/",
    #                          paste0(combos[1,x], ".vs.", combos[2,x]),
    #                          "/",
    #                          paste0(combos[1,x], ".vs.", combos[2,x]),
    #                          "_1point5FC.csv"),
    #   # need a column called "locusName", "adj.P.Val" and "logFC",
    #   out_dir = paste0(project_paths[3],"/Combined/GO_tests/")
    # )
  })

  data_folder <- file.path(paste0(project_paths[3], "/Combined/GO_tests/"))   # path to the data

  # make list of files
  file_names <- list.files(path = data_folder, pattern = "*.csv")

  # make list of file paths
  sample_paths <- file.path(data_folder, file_names)

  ############ split file names to match samples to keyfile

  Sample_names <- file_names %>% str_replace(".csv", "")

  data <- data_frame(sample = Sample_names, paths = sample_paths) %>%
    mutate(
      file_contents = map(sample_paths, ~ read_csv(., col_names = TRUE))
    ) %>%
    dplyr::select(-paths)

  big_data <- unnest(data[which(!isEmpty(data$file_contents)),],
                     cols = c(file_contents))

  sapply(1:combos[1, ], function (x) {
    map2(
      rep(c("BP", "MF", "CC"), 2), rep(c("up", "down"), each = 3),
      analyse_topgo,
      test_name  = paste0(combos[1,x], ".vs.", combos[2,x]),
      big_data   = big_data,
      data_folder = data_folder
    )
  })
}
